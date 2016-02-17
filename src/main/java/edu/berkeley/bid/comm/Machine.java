package edu.berkeley.bid.comm;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.net.ConnectException;
import java.net.InetSocketAddress;
import java.net.ServerSocket;
import java.net.Socket;
import java.net.SocketException;
import java.net.SocketTimeoutException;
import java.nio.ByteBuffer;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class Machine {
	/* Machine Configuration Variables */	
	final int N;                                               // Number of features
	final int D;                                               // Depth of the network
	final int M;                                               // Number of Machines
	final int imachine;                                        // My identity
	final Groups groups;                                       // group membership indices
	final int replicate;                                       // replication factor
	final String [] machineIP;                                 // String IP names
	final boolean doSim;                                       // Simulation on one machine: send messages directly without sockets
	int sockBase = 50000;                                      // Socket base address
	int sendTimeout = 1000;                                    // in msec
	int trace = 0;                                             // 0: no trace, 1: high-level, 2: everything

	Layer [] layers;                                           // All the layers
	ByteBuffer [] sendbuf;                                     // buffers, one for each destination in a group
	ByteBuffer [] recbuf;
	IVec finalMap;                                             // Map from down --> up at layer D-1
	Msg [][] messages;                                         // Message queue 
	boolean [][] msgrecvd;                                     // Receiver status
	boolean [][] amsending;                                    // Sender status
	ExecutorService executor;
	ExecutorService sockExecutor;
	Listener listener;
	AllReduceY parent;

	public Machine(AllReduceY p0, int N0, Groups groups0, int imachine0, int M0, int bufsize, boolean doSim0, int trace0, 
			int replicate0, String [] machineIP0) {
		parent = p0;
		N = N0;
		M = M0;
		doSim = doSim0;
		imachine = imachine0;
		groups = groups0;
		replicate = replicate0;
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
		for (int i = 0; i < D; i++) {
			int k = groups.nodesInGroup(i, imachine).length;
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
			cumPos = cumPos * k + pos;
			layers[level] = new Layer(this, k, cumk, cumPos, pos, level);
			cumk *= k;
		}
		if (!doSim) {
			sockExecutor = Executors.newFixedThreadPool(1+4*maxk); 
			listener = new Listener();
			sockExecutor.execute(listener);
		}
	}

	public void stop() {
		if (listener != null) {
			listener.stop();
		}
	}

	public void config(IVec downi, IVec upi) {
		IVec [] outputs = new IVec[2];
		for (int i = 0; i < D; i++) {
			layers[i].config(downi, upi, outputs);
			downi = outputs[0];
			upi = outputs[1];
		}
		finalMap = IVec.mapInds(upi, downi);
	}

	public Vec reduce(Vec downv, int stride) {
		for (int d = 0; d < D; d++) {
			downv = layers[d].reduceDown(downv, stride);
		}
		Vec upv = downv.mapFrom(finalMap);
		for (int d = D-1; d >= 0; d--) {
			upv = layers[d].reduceUp(upv, stride);
		}
		if (trace > 0) {
			synchronized (parent) {
				System.out.format("machine %d reduce result nnz %d out of %d\n", imachine, upv.nnz(), upv.size());
			}
		}
		return upv;
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
					Socket cs = ss.accept();
					sockExecutor.execute(new SockReader(cs));
				} catch (SocketException e) {
					// This is probably due to the server shutting down. Don't do anything.
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
					parent.simNetwork[outi + i*M].messages[imachine][tag] = msg;
				}
			} else {
				for (int i = 0; i < replicate; i++) { 
					sockExecutor.execute(new SockWriter(outi + i*M, msg));
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