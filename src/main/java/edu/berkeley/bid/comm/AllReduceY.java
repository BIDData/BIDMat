
package edu.berkeley.bid.comm;

// Third version of sparse allreduce
// Includes support for matrices, feature range limits, long feature indices
//
// 

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.Arrays;
import java.io.*;
import java.nio.*;
import java.net.*;
//import mpi.*;


public class AllReduceY {

	public class Machine {
		/* Machine Configuration Variables */	
		final int N;                                               // Number of features
		final int D;                                               // Depth of the network
		final int M;                                               // Number of Machines
		final int imachine;                                        // My identity
		final int [][] groups;                                     // group membership indices
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

		// return the size of the group at level "depth" containing this node, and the position of this node in that group.

		public int [] inMyGroup(int imachine, int depth) {
			int [] grouplev = groups[depth];
			int count = 0;
			int pos = 0;
			int mygroup = grouplev[imachine];
			for (int i = 0; i < N; i++) {
				if (imachine == i) pos = count;
				if (grouplev[i] == mygroup) count++;
			}
			return new int[]{count, pos};
		}

		// return the machines in the group at level "depth" containing this node, and the position of this node in that group.

		public int [] nodesInMyGroup(int imachine, int depth) {
			int count = 0;
			int [] grouplev = groups[depth];
			int mygroup = grouplev[imachine];
			int [] img = inMyGroup(imachine, depth);
			int [] machines = new int[img[0]];
			for (int i = 0; i < N; i++) {
				if (grouplev[i] == mygroup) {
					machines[count] = i;
					count++;
				}
			}
			return machines;
		}

		public Machine(int N0, int [][] groups0, int imachine0, int M0, int bufsize, boolean doSim0, int trace0, 
				int replicate0, String [] machineIP0) {
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
			D = groups.length;
			trace = trace0;
			layers = new Layer[D];
			int cumk = 1;
			int maxk = 1;
			int cumPos = 0;
			for (int i = 0; i < D; i++) {
				int [] inGroup = inMyGroup(imachine, i);
				int k = inGroup[0];
				int pos = inGroup[1];
				cumPos = cumPos * k + pos;
				layers[i] = new Layer(k, cumk, cumPos, pos, imachine, i);
				cumk *= k;
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
				synchronized (AllReduceY.this) {
					System.out.format("machine %d reduce result nnz %d out of %d\n", imachine, upv.nnz(), upv.size());
				}
			}
			return upv;
		}

	

		class Layer {

			/* Layer Configuration Variables */

			int k;        																					 // Size of this group
			int cumk;                                                // cumulative product of sizes
			int cumPos;                                              // cumulative position in groups
			int depth;
			int posInMyGroup;                                        // Position in this machines group
			int [] outNbr;                                           // Machines we talk to 
			int [] inNbr;                                            // Machines we listen to
			IVec [] downMaps;                                        // Maps to indices below for down indices
			IVec [] upMaps;                                          // Maps to indices below for up indices
			IVec [] interleave;                                      // Indices to undo partitioning for up data
			IVec dpart;                                              // partition number for each element
			IVec upart;
			IVec [] downpartinds;
			IVec [] uppartinds;
			Vec [] upparts;
			int downn;                                               // Size of the down master list
			int upn;                                                 // Size of the up vector
			Partitioner partitioner;

			public Layer(int k0, int cumk0, int cumPos0, int posInMyGroup0, int imachine, int depth0) {
				k = k0;
				cumk = cumk0;
				cumPos = cumPos0;
				posInMyGroup = posInMyGroup0;
				depth = depth0;
				inNbr = nodesInMyGroup(imachine, depth);
				outNbr = nodesInMyGroup(imachine, depth);		
				downMaps = new IVec[k];
				upMaps = new IVec[outNbr.length];
			}		

			class ConfigThread implements Runnable {
				IVec [] downp;
				IVec [] upp;
				IVec [] dtree;
				IVec [] utree;
				int i;
				int repno;
				CountDownLatch latch;

				public ConfigThread(IVec [] downp0, IVec [] upp0,	IVec [] dtree0, IVec [] utree0, int i0, CountDownLatch latch0) {
					downp = downp0;
					upp = upp0;					
					dtree = dtree0;
					utree = utree0;
					i = i0;
					latch = latch0;
				}

				public void run() {
					sendbuf[i].clear();
					recbuf[i].clear();
					IntBuffer sbuf = sendbuf[i].asIntBuffer();
					IntBuffer rbuf = recbuf[i].asIntBuffer();	
					int seg1 = downp[i].size();
					int seg2 = seg1 + upp[i].size();
					sbuf.put(seg1);
					sbuf.put(seg2);
					sbuf.put(downp[i].data, 0, seg1);
					sbuf.put(upp[i].data, 0, seg2-seg1);

					if (trace > 1) {
						synchronized (AllReduceY.this) {
							System.out.format("config layer %d machine %d sent msg to %d, from %d, sizes %d %d\n", depth, imachine, outNbr[i],  inNbr[i],  sbuf.get(0), sbuf.get(1));
						}
					}
					sendrecv(i, sendbuf, seg2+2, outNbr[i], recbuf, rbuf.capacity(), inNbr[i], depth*3);
					seg1 = rbuf.get();
					seg2 = rbuf.get();
					if (trace > 1) {
						synchronized (AllReduceY.this) {
							System.out.format("config layer %d machine %d got msg from %d, sizes %d %d\n", depth, imachine, inNbr[i], seg1, seg2);
						}
					}

					IVec downout = new IVec(seg1);
					IVec upout =   new IVec(seg2-seg1);
					rbuf.get(downout.data, 0, seg1);
					rbuf.get(upout.data, 0, seg2-seg1);
					IVec.checkTree(dtree, downout, i, k);
					IVec.checkTree(utree, upout, i, k);	
					downp[i] = downout;
					upp[i] = upout;

					latch.countDown();
				}
			}

			public void config(IVec downi, IVec upi, IVec [] outputs) {
				dpart = partitioner.part(downi);
				upart = partitioner.part(upi);
				downpartinds = partitioner.partition(downi, dpart, null);
				uppartinds = partitioner.partition(upi, upart, interleave);
				IVec [] dtree = new IVec[2*k-1];
				IVec [] utree = new IVec[2*k-1];
				if (trace > 0) {
					synchronized (AllReduceY.this) {
						System.out.format("machine %d layer %d, dparts (%d", imachine, depth, downpartinds[0].size());
					}
				}

				CountDownLatch latch = new CountDownLatch(k);
				for (int i = 0; i < k; i++) {
					int ix = (i + posInMyGroup) % k;                               // Try to stagger the traffic
					executor.execute(new ConfigThread(downpartinds, uppartinds, dtree, utree, ix, latch));
				}
				try {	latch.await(); } catch (InterruptedException e) {}
				IVec dmaster = dtree[0];
				Arrays.fill(dtree, null);
				downn = dmaster.size();
				IVec umaster = utree[0];
				Arrays.fill(utree, null);
				upn = upi.size();
				for (int i = 0; i < k; i++) {
					downMaps[i] = IVec.mapInds(downpartinds[i], dmaster);
					upMaps[i] = IVec.mapInds(uppartinds[i], umaster);
					if (trace > 0) {
						synchronized (AllReduceY.this) { 
							System.out.format("machine %d dmap(%d) size %d\n", imachine, i, downMaps[i].size());
							System.out.format("machine %d umap(%d) size %d\n", imachine, i, upMaps[i].size());
						}
					}
				}
				outputs[0] = dmaster;
				outputs[1] = umaster;
			}

			public class ReduceDownThread implements Runnable {
				Vec newv;
				Vec downv;
				int i;
				CountDownLatch latch;
				
				// Note, these have to be partitions of the data now. 

				public ReduceDownThread(Vec newv0, Vec downv0, int i0, CountDownLatch latch0) {
					newv = newv0;
					downv = downv0;
					i = i0;
					latch = latch0;
				}

				public void run() {
					sendbuf[i].clear();
					recbuf[i].clear();
					IntBuffer isbuf = sendbuf[i].asIntBuffer();
					IntBuffer irbuf = recbuf[i].asIntBuffer();
					FloatBuffer sbuf = sendbuf[i].asFloatBuffer();
					FloatBuffer rbuf = recbuf[i].asFloatBuffer();
					int msize = downv.size();
					isbuf.put(msize);
					sbuf.position(1);
					sbuf.put(downv.data, 0, msize);

					if (trace > 1) {
						synchronized (AllReduceY.this) {
							System.out.format("reduce layer %d machine %d sent msg to %d, from %d, size %d\n", depth, imachine, outNbr[i],  inNbr[i],  msize);
						}
					}
					sendrecv(i, sendbuf, msize+1, outNbr[i], recbuf, rbuf.capacity(), inNbr[i], depth*3+1);
					msize = irbuf.get();
					if (trace > 1) {
						synchronized (AllReduceY.this) {
							System.out.format("reduce layer %d machine %d got msg from %d, size %d\n", depth, imachine, inNbr[i], msize);
						}
					}

					Vec res = new Vec(msize);
					rbuf.position(1);
					rbuf.get(res.data, 0, msize);
					synchronized (newv) {
						res.addTo(newv, downMaps[i]);	
					}
					latch.countDown();
				}
			}

			public Vec reduceDown(Vec downv, int stride) {
				CountDownLatch latch = new CountDownLatch(k);
				Vec newv = new Vec(downn);
				Vec [] vparts = partitioner.partition(downv, dpart, downpartinds, stride);
				for (int i = 0; i < k; i++) {
					int ix = (i + posInMyGroup) % k;   // Try to stagger the traffic
					executor.execute(new ReduceDownThread(newv, vparts[ix], ix, latch));
				}
				try { latch.await(); } catch (InterruptedException e) {}
				return newv;
			}

			public class ReduceUpThread implements Runnable {
				Vec newv;
				Vec upv;
				int i;
				CountDownLatch latch;

				public ReduceUpThread(Vec newv0, Vec upv0, int i0, CountDownLatch latch0) {
					newv = newv0;
					upv = upv0;
					i = i0;
					latch = latch0;
				}

				public void run () {
					sendbuf[i].clear();
					recbuf[i].clear();
					IntBuffer isbuf = sendbuf[i].asIntBuffer();
					IntBuffer irbuf = recbuf[i].asIntBuffer();
					FloatBuffer sbuf = sendbuf[i].asFloatBuffer();
					FloatBuffer rbuf = recbuf[i].asFloatBuffer();
					Vec up = upv.mapFrom(upMaps[i]);
					int msize = up.size();
					isbuf.put(msize);
					sbuf.position(1);
					sbuf.put(up.data, 0, msize);

					if (trace > 1) {
						synchronized (AllReduceY.this) {
							System.out.format("reduce up layer %d machine %d sent msg to %d, from %d, size %d\n", depth, imachine, outNbr[i],  inNbr[i],  msize);
						}
					}
					sendrecv(i, sendbuf, msize+1, inNbr[i], recbuf, irbuf.capacity(), outNbr[i], depth*3+2);
					msize = irbuf.get();
					if (trace > 1) {
						synchronized (AllReduceY.this) {
							System.out.format("reduce up layer %d machine %d got msg from %d, size %d\n", depth, imachine, inNbr[i], msize);
						}
					}

					int psize = upMaps[i].size();
					if (msize != psize) throw new RuntimeException("ReduceUp size mismatch "+msize+" "+psize);
					rbuf.position(1);
					rbuf.get(upparts[i].data, 0, msize);
					latch.countDown();
				}
			}

			public Vec reduceUp(Vec upv, int stride) {
				Vec newv = new Vec(upn);
				CountDownLatch latch = new CountDownLatch(k);
				for (int i = 0; i < k; i++) {
					int ix = (i + posInMyGroup) % k;                               // Try to stagger the traffic
					executor.execute(new ReduceUpThread(newv, upv, ix, latch));
				}
				try { latch.await(); } catch (InterruptedException e) {}
				partitioner.merge(newv, stride, upparts, interleave);
				return newv;
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
						simNetwork[outi + i*M].messages[imachine][tag] = msg;
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
	}

	public class Msg {
		byte [] buf;
		int size;
		int sender;
		int receiver;
		int tag;

		public Msg(int size0, int sender0, int receiver0, int tag0) {
			buf = new byte[4*size0];
			size = size0;
			sender = sender0;
			receiver = receiver0;
			tag = tag0;
		}

		public Msg(byte [] inbuf, int size0, int sender0, int receiver0, int tag0) {
			buf = new byte[4*size0];
			System.arraycopy(inbuf, 0, buf, 0, 4*size0);
			size = size0;
			sender = sender0;
			receiver = receiver0;
			tag = tag0;
		}
	}

	public Machine [] simNetwork = null;

	public AllReduceY(int M) {
		simNetwork = new Machine[M];
	}

	public void stop() {
		if (simNetwork != null) {
			for (int i = 0; i < simNetwork.length; i++) simNetwork[i].stop();
		}
	}

}
