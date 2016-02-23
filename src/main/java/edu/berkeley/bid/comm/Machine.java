package edu.berkeley.bid.comm;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.net.InetSocketAddress;
import java.net.ServerSocket;
import java.net.Socket;
import java.net.SocketException;
import java.net.SocketTimeoutException;
import java.net.ConnectException;
import java.nio.ByteBuffer;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.Hashtable;

public class Machine {
	/* Machine Configuration Variables */	
	public final int D;                                               // Depth of the network
	public final int M;                                               // Number of Machines
	public final int imachine;                                        // My identity
	public final Groups groups;                                       // group membership indices
	public final int replicate;                                       // replication factor
	public final String [] machineIP;                                 // String IP names
	public final boolean doSim;                                       // Simulation on one machine: send messages directly without sockets
	public int sockBase = 50000;                                      // Socket base address
	public int sendTimeout = 1000;                                    // in msec
	public int trace = 0;                                             // 0: no trace, 1: high-level, 2: everything
	public int configTimeout = 1000;
	public int reduceTimeout = 1000;

	public Vec tov;                                                   // Snapshot of the reduced matrix at each layer going down
	public Vec fromv;                                                 // Snapshot of the result matrix coming up
	public Layer [] layers;                                           // All the layers
	public ByteBuffer [] sendbuf;                                     // buffers, one for each destination in a group
	public ByteBuffer [] recbuf;
	public IVec finalMap;                                             // Map from to --> up at layer D-1
	public Msg [][] messages;                                         // Message queue 
	public boolean [][] msgrecvd;                                     // Receiver status
	public boolean [][] amsending;                                    // Sender status
	public ExecutorService executor;
	public ExecutorService sockExecutor;
	public Future<?> listener;
	public Hashtable<Machine.SockWriter, Future<?>> writers;
	public Hashtable<Machine.SockReader, Future<?>> readers;
	Network network;

	public Machine(Network p0, Groups groups0, int imachine0, int M0, int bufsize, boolean doSim0, int trace0, 
			int replicate0, String [] machineIP0) {
		network = p0;
		M = M0;
		doSim = doSim0;
		imachine = imachine0;
		groups = groups0;
		replicate = replicate0;
		writers = new Hashtable<Machine.SockWriter, Future<?>>();
		readers = new Hashtable<Machine.SockReader, Future<?>>();
		if (machineIP0 == null) {
			machineIP = new String[M*replicate];
			for (int i = 0; i < M*replicate; i++) machineIP[i] = "localhost";
		} else {
			machineIP = machineIP0;
		}
		D = groups.depth();
		trace = trace0;
		layers = new Layer[D];
		int cumk = 1;
		int maxk = 1;
		int cumPos = 0;
		for (int level = 0; level < D; level++) {
			int k = groups.nodesInGroup(imachine, level).length;
			maxk = Math.max(maxk, k);
		}
		executor = Executors.newFixedThreadPool(maxk); // set to 1 for sequential messaging. 
		sendbuf = new ByteBuffer[maxk];
		recbuf = new ByteBuffer[maxk];
		for (int i = 0; i < maxk; i++) {
			sendbuf[i] = ByteBuffer.wrap(new byte[4*bufsize]);
			recbuf[i] = ByteBuffer.wrap(new byte[4*bufsize]);
		}
		messages = new Msg[M*replicate][];
		msgrecvd = new boolean[M*replicate][];
		amsending = new boolean[M*replicate][];
		for (int i = 0; i < M*replicate; i++) {
			messages[i] = new Msg[3*D];
			msgrecvd[i] = new boolean[3*D];
			amsending[i] = new boolean[3*D];
		}
		for (int level = 0; level < D; level++) {
			int k = groups.nodesInGroup(imachine, level).length;
			int pos = groups.posInGroup(imachine, level);
			layers[level] = new Layer(this, k, cumk, cumPos, pos, level);
			cumPos = cumPos * k + pos;
			cumk *= k;
		}
		if (!doSim) {
			sockExecutor = Executors.newFixedThreadPool(1+4*maxk); 
			listener = sockExecutor.submit(new Listener());
		}
	}

	public void stop() {
		executor.shutdownNow();
		if (sockExecutor != null) {
			listener.cancel(true);
			sockExecutor.shutdownNow();
		}
	}

	public void config(IVec toi, IVec upi) {
		IVec [] outputs = new IVec[2];
		for (int i = 0; i < D; i++) {
			layers[i].config(toi, upi, outputs);
			toi = outputs[0];
			upi = outputs[1];
		}
		finalMap = IVec.mapInds(upi, toi);
	}
	
	public void config(int [] toi, int [] upi) {
		config(new IVec(toi), new IVec(upi));
	}
	

	public Vec reduce(Vec tov0, int stride) {
		tov = tov0;
		for (int d = 0; d < D; d++) {
			tov = layers[d].reduceDown(tov, stride);
		}
		fromv = tov.mapFrom(finalMap, stride);
		for (int d = D-1; d >= 0; d--) {
			fromv = layers[d].reduceUp(fromv, stride);
		}
		if (trace > 0) {
			synchronized (network) {
				System.out.format("Reduce machine %d result nnz %d out of %d\n", imachine, fromv.nnz(), fromv.size());
			}
		}
		return fromv;
	}
	
	public float [] reduce(float [] tov, int stride) {
		return reduce(new Vec(tov), stride).data;
	}
	
	public Vec configReduce(IVec toi, IVec fromi, Vec tov0, int stride) {
		tov = tov0;
		IVec [] outputs = new IVec[2];
		for (int d = 0; d < D; d++) {
			tov = layers[d].configReduce(toi, fromi, outputs, tov, stride);
			toi = outputs[0];
			fromi = outputs[1];
		}
		finalMap = IVec.mapInds(fromi, toi);
		fromv = tov.mapFrom(finalMap, stride);
		for (int d = D-1; d >= 0; d--) {
			fromv = layers[d].reduceUp(fromv, stride);
		}
		if (trace > 0) {
			synchronized (network) {
				System.out.format("ConfigReduce machine %d result nnz %d out of %d\n", imachine, fromv.nnz(), fromv.size());
			}
		}
		return fromv;
	}
	
	public float [] configReduce(int [] toi, int [] fromi, float [] tov, int stride) {
		return configReduce(new IVec(toi), new IVec(fromi), new Vec(tov), stride).data;
	}



	public class SockWriter implements Runnable {
		int dest;
		Msg msg;

		public SockWriter(int dest0, Msg msg0) {
			msg = msg0;
			dest = dest0;
		}

		public void run() {
			Socket socket = null;
			try {
				socket = new Socket();
				socket.connect(new InetSocketAddress(machineIP[dest], sockBase + dest), sendTimeout);
				if (socket.isConnected()) {
					amsending[dest][msg.tag] = true;
					DataOutputStream ostr = new DataOutputStream(socket.getOutputStream());
					ostr.writeInt(msg.size);
					ostr.writeInt(msg.sender);
					ostr.writeInt(msg.tag);
					ostr.write(msg.buf, 0, msg.size*4);		
				}
			}	catch (SocketTimeoutException e) {
				// No need to do anything
			} catch (ConnectException e) {
				// Server may have been killed - OK
			}	catch (SocketException e) {
				// Server may have been killed - OK
			}	catch (Exception e) {
				throw new RuntimeException("Problem writing socket "+e);
			} finally {
				try { if (socket != null) socket.close(); } catch (Exception e) {}
				amsending[dest][msg.tag] = false;
				synchronized (writers) {
					if (writers.containsKey(this)) writers.remove(this);
				}
			}
		}
	}

	public class SockReader implements Runnable {
		Socket socket = null;

		public SockReader(Socket sock) {
			socket = sock;
		}

		public void run() {
			try {
				DataInputStream istr = new DataInputStream(socket.getInputStream());
				int len = istr.readInt();
				int src = istr.readInt();
				int tag = istr.readInt();
				if (!msgrecvd[src][tag]) {
					Msg msg = new Msg(len, src, imachine, tag);
					istr.readFully(msg.buf, 0, len*4);
					if (!msgrecvd[src][tag]) {
						messages[src][tag] = msg;		
					}
				}
			} catch (Exception e) {
				throw new RuntimeException("Problem reading socket "+e);
			} finally {
				try {socket.close();} catch (IOException e) {}
				synchronized (readers) {
					if (readers.containsKey(this)) readers.remove(this);
				}
			}
		}
	}

	public class Listener implements Runnable {
		boolean stop = false;
		ServerSocket ss = null;

		public Listener() {
			try {
				ss = new ServerSocket(sockBase + imachine);
			} catch (Exception e) {
				throw new RuntimeException("Couldnt start socket listener "+e);
			}			
		}

		public void run() {
			while (!stop) {
				try {
					SockReader scs = new SockReader(ss.accept());
					Future<?> fut = sockExecutor.submit(scs);
					synchronized (readers) {
						readers.put(scs, fut);
					}
				} catch (SocketException e) {
					// This is probably due to the server shutting to. Don't do anything.
				}
				catch (Exception e) {
					throw new RuntimeException("Socket listener had a problem "+e);
				}
			}
		}

		public boolean stillSending() {
			boolean sending = false;
			for (int i = 0; i < amsending.length; i++) {
				boolean [] sendrow = amsending[i];
				for (int j = 0; j < sendrow.length; j++) {
					if (amsending[i][j]) sending = true;
				}
			}			
			return sending;
		}

		public void stop() {
			while (stillSending()) {
				try { Thread.sleep(1); } catch (InterruptedException e) {}
			}
			try {
				stop = true;
				ss.close();
			} catch (Exception e) {
				throw new RuntimeException("Trouble closing listener");
			}			
		}
	}
	
	public boolean sendrecv(int igroup, ByteBuffer [] sbuf, int sendn, int outi, ByteBuffer [] rbuf, int recn, int ini, int tag) {
		sbuf[igroup].rewind();
		Msg msg = new Msg(sbuf[igroup].array(), sendn, imachine, outi, tag);
		if (imachine == outi) {
			rbuf[igroup].clear();
			rbuf[igroup].put(msg.buf, 0, 4*sendn);
			rbuf[igroup].rewind();
			return true;				
		} else { 
			if (doSim) {					
				for (int i = 0; i < replicate; i++) { 
					network.machines[outi + i*M].messages[imachine][tag] = msg;
				}
			} else {
				for (int i = 0; i < replicate; i++) { 
					SockWriter w = new SockWriter(outi + i*M, msg);
					Future<?> fut = sockExecutor.submit(w);
					synchronized (writers) {
						writers.put(w, fut);
					}
				}
			}
			boolean gotit = false;
			while (!gotit) {
				for (int i = 0; i < replicate; i++) {
					if (messages[ini + i*M][tag] != null) {
						Msg rmsg = messages[ini + i*M][tag];
						rbuf[igroup].clear();
						rbuf[igroup].put(rmsg.buf, 0, 4*rmsg.size);
						rbuf[igroup].rewind();
						gotit = true;
						break;
					}
					try {
						Thread.sleep(1);
					} catch (InterruptedException e) {}
				}
			}
			for (int i = 0; i < replicate; i++) {
				messages[ini + i*M][tag] = null;
				msgrecvd[ini + i*M][tag] = true;
			}
			return true;
		}
	}	
}