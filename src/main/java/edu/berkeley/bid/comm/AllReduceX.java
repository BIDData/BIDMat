
package edu.berkeley.bid.comm;

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.Arrays;
import java.io.*;
import java.nio.*;
import java.net.*;
//import mpi.*;

public class AllReduceX {
	
	public class Machine {
		/* Machine Configuration Variables */	
		int N;                                                     // Number of features
		int D;                                                     // Depth of the network
		int M;                                                     // Number of Machines
		int imachine;                                              // My identity
		int [] allks;                                              // k values
		int replicate = 1;                                         // replication factor
		int sockBase = 50000;                                      // Socket base address
		int sendTimeout = 1000;                                    // in msec
		String [] machineIP;                                       // String IP names
		
		Layer [] layers;                                           // All the layers
		ByteBuffer [] sendbuf;                                     // buffers, one for each destination in a group
		ByteBuffer [] recbuf;
		IVec finalMap;                                             // Map to down from down --> up at layer D-1
		Msg [][] messages;                                         // Message queue for the simulation
		boolean [][] msgrecvd;
		boolean doSim = true;
		ExecutorService executor;
		ExecutorService sockExecutor;
		Listener listener;
		int trace = 0;                                             // 0: no trace, 1: high-level, 2: everything

		public Machine(int N0, int [] allks0, int imachine0, int M0, int bufsize, boolean doSim0, int trace0, 
				           int replicate0, String [] machineIP0) {
			N = N0;
			M = M0;
			doSim = doSim0;
			replicate = replicate0;
			machineIP = machineIP0;
			if (machineIP == null) {
				machineIP = new String[M*replicate];
				for (int i = 0; i < M*replicate; i++) machineIP[i] = "localhost";
			}
			imachine = imachine0;
			allks = allks0;
			D = allks.length;
			trace = trace0;
			layers = new Layer[D];
			int left = 0;
			int right = N;
			int cumk = 1;
			int maxk = 1;
			for (int i = 0; i < D; i++) {
				int k = allks[i];
				layers[i] = new Layer(k, cumk, left, right, imachine, i);
				int pimg = layers[i].posInMyGroup;
				left = layers[i].left;
				if (pimg > 0) left = layers[i].partBoundaries.data[pimg-1];
				right = layers[i].partBoundaries.data[pimg];
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
			for (int i = 0; i < M*replicate; i++) {
				messages[i] = new Msg[3*D];
				msgrecvd[i] = new boolean[3*D];
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
		
		public Vec reduce(Vec downv) {
			for (int d = 0; d < D; d++) {
				downv = layers[d].reduceDown(downv);
			}
			Vec upv = downv.mapFrom(finalMap);
			for (int d = D-1; d >= 0; d--) {
				upv = layers[d].reduceUp(upv);
			}
			if (trace > 0) {
				synchronized (AllReduceX.this) {
					System.out.format("machine %d reduce result nnz %d out of %d\n", imachine, upv.nnz(), upv.size());
				}
			}
			return upv;
		}
	
		class Layer {
			
			/* Layer Configuration Variables */
			
			int k;        																					 // Size of this group
			int left;                                                // Left boundary of its indices
			int right;                                               // Right boundary of its indices
			int depth;
			int posInMyGroup;                                        // Position in this machines group
			int [] outNbr;                                           // Machines we talk to 
			int [] inNbr;                                            // Machines we listen to
			IVec partBoundaries;                                     // Partition boundaries
			IVec [] downMaps;                                        // Maps to indices below for down indices
			IVec [] upMaps;                                          // Maps to indices below for up indices
			int downn;                                               // Size of the down master list
			int upn;                                                 // Size of the up vector
			int [] dPartInds;
			int [] uPartInds;

			public Layer(int k0, int cumk, int left0, int right0, int imachine, int depth0) {
				k = k0;
				int i;
				left = left0;
				right = right0;
				depth = depth0;
				partBoundaries = new IVec(k);
				inNbr = new int [k];
				outNbr = new int [k];
				dPartInds = new int[k+1];
				uPartInds = new int[k+1];
				int ckk = cumk * k;
				posInMyGroup = (imachine % ckk) / cumk;
				int ibase = (imachine % M) - posInMyGroup * cumk;
				for (i = 0; i < k; i++) {
					partBoundaries.data[i] = left + (int)(((long)(right - left)) * (i+1) / k);
					outNbr[i] = ibase + i * cumk;
					int toMe = (k + 2*posInMyGroup - i) % k;
					inNbr[i] = ibase + toMe * cumk;
				}		
				downMaps = new IVec[k];
				upMaps = new IVec[k];
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
						synchronized (AllReduceX.this) {
							System.out.format("config layer %d machine %d sent msg to %d, from %d, sizes %d %d\n", depth, imachine, outNbr[i],  inNbr[i],  sbuf.get(0), sbuf.get(1));
						}
					}
					sendrecv(i, sendbuf, seg2+2, outNbr[i], recbuf, rbuf.capacity(), inNbr[i], depth*3);
					seg1 = rbuf.get();
					seg2 = rbuf.get();
					if (trace > 1) {
						synchronized (AllReduceX.this) {
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
				IVec [] downp = IVec.partition(downi, partBoundaries);
				IVec [] upp = IVec.partition(upi, partBoundaries);
				IVec [] dtree = new IVec[2*k-1];
				IVec [] utree = new IVec[2*k-1];
				if (trace > 0) {
					synchronized (AllReduceX.this) {
						System.out.format("machine %d layer %d, dparts (%d", imachine, depth, downp[0].size());
						for (int i = 1; i < downp.length; i++) System.out.format(", %d", downp[i].size());
						System.out.format(") from %d, bounds %d %d\n", downi.size(), partBoundaries.data[0], partBoundaries.data[partBoundaries.size()-1]);
						System.out.format("machine %d layer %d, uparts (%d", imachine, depth, upp[0].size());
						for (int i = 1; i < upp.length; i++) System.out.format(", %d", upp[i].size());
						System.out.format(") from %d, bounds %d %d\n", upi.size(), partBoundaries.data[0], partBoundaries.data[partBoundaries.size()-1]);
					}
				}
				dPartInds[0] = 0;
				uPartInds[0] = 0;
				for (int i = 0; i < k; i++) {
					dPartInds[i+1] = dPartInds[i] + downp[i].size();
					uPartInds[i+1] = uPartInds[i] + upp[i].size();
				}
				CountDownLatch latch = new CountDownLatch(k);
				for (int i = 0; i < k; i++) {
					int ix = (i + posInMyGroup) % k;                               // Try to stagger the traffic
					executor.execute(new ConfigThread(downp, upp, dtree, utree, ix, latch));
				}
				try {	latch.await(); } catch (InterruptedException e) {}
				IVec dmaster = dtree[0];
				Arrays.fill(dtree, null);
				downn = dmaster.size();
				IVec umaster = utree[0];
				Arrays.fill(utree, null);
				upn = upi.size();
				for (int i = 0; i < k; i++) {
					downMaps[i] = IVec.mapInds(downp[i], dmaster);
					upMaps[i] = IVec.mapInds(upp[i], umaster);
  				if (trace > 0) {
  					synchronized (AllReduceX.this) { 
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
					int msize = dPartInds[i+1] - dPartInds[i];
					isbuf.put(msize);
					sbuf.position(1);
					sbuf.put(downv.data, dPartInds[i], msize);
					
					if (trace > 1) {
						synchronized (AllReduceX.this) {
							System.out.format("reduce layer %d machine %d sent msg to %d, from %d, size %d\n", depth, imachine, outNbr[i],  inNbr[i],  msize);
						}
					}
					sendrecv(i, sendbuf, msize+1, outNbr[i], recbuf, rbuf.capacity(), inNbr[i], depth*3+1);
					msize = irbuf.get();
					if (trace > 1) {
						synchronized (AllReduceX.this) {
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
			
			public Vec reduceDown(Vec downv) {
				Vec newv = new Vec(downn);
				CountDownLatch latch = new CountDownLatch(k);
				for (int i = 0; i < k; i++) {
					int ix = (i + posInMyGroup) % k;                               // Try to stagger the traffic
					executor.execute(new ReduceDownThread(newv, downv, ix, latch));
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
						synchronized (AllReduceX.this) {
							System.out.format("reduce up layer %d machine %d sent msg to %d, from %d, size %d\n", depth, imachine, outNbr[i],  inNbr[i],  msize);
						}
					}
					sendrecv(i, sendbuf, msize+1, inNbr[i], recbuf, irbuf.capacity(), outNbr[i], depth*3+2);
					msize = irbuf.get();
					if (trace > 1) {
						synchronized (AllReduceX.this) {
							System.out.format("reduce up layer %d machine %d got msg from %d, size %d\n", depth, imachine, inNbr[i], msize);
						}
					}
					
					int psize = uPartInds[i+1] - uPartInds[i];
					if (uPartInds[i+1] > newv.size()) throw new RuntimeException("ReduceUp index out of range "+uPartInds[i+1]+" "+newv.size());
					if (msize != psize) throw new RuntimeException("ReduceUp size mismatch "+msize+" "+psize);
					rbuf.position(1);
					rbuf.get(newv.data, uPartInds[i], msize);
					latch.countDown();
				}
			}
			
			public Vec reduceUp(Vec upv) {
				Vec newv = new Vec(upn);
				CountDownLatch latch = new CountDownLatch(k);
				for (int i = 0; i < k; i++) {
					int ix = (i + posInMyGroup) % k;                               // Try to stagger the traffic
					executor.execute(new ReduceUpThread(newv, upv, ix, latch));
				}
				try { latch.await(); } catch (InterruptedException e) {}
				return newv;
			}
		}

		public boolean sendrecv(int igroup, ByteBuffer [] sbuf, int sendn, int outi, ByteBuffer [] rbuf, int recn, int ini, int tag) {
			sbuf[igroup].rewind();
			Msg msg = new Msg(sbuf[igroup].array(), sendn, imachine, outi, tag);
			if (imachine == outi) {
				rbuf[igroup].clear();
				rbuf[igroup].put(msg.buf, 0, 4*sendn);
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
					
/*					try {
            sbuf.rewind();
            rbuf.clear();
						MPI.COMM_WORLD.sendRecv(sbuf, 4*sendn, MPI.BYTE, outi, tag, rbuf, 4*recn, MPI.BYTE, ini, tag); 
						sbuf.rewind();
						rbuf.rewind();
					} catch (MPIException e) {
						throw new RuntimeException("Exception in sendrecv "+e);
					} */
					
					// JFC: Use this code
				/*	try {	 
            Request [] sreq = new Request[replicate];
            Request [] rreq = new Request[replicate];
            boolean sdone = false;
            boolean rdone = false;
            for (int i = 0; i < replicate; i++) {
              sbuf[igroup + i*k].rewind();
              rbuf[igroup + i*k].clear();  
            	if (i > 0) {
            		sbuf[igroup + i*k].put(sbuf[igroup].array(), 0, sendn);
            	}
            	sreq[i] = MPI.COMM_WORLD.iSend(sbuf[igroup + i*k].array(), 4*sendn, MPI.BYTE, outi + i*M, tag);
            	rreq[i] = MPI.COMM_WORLD.iRecv(rbuf[igroup + i*k].array(), 4*recn, MPI.BYTE, ini + i*M, tag);
            }
            // Wait until timeout or when one send and one receive are done, then cancel others
						long timeout = 2000;   // Wait this many msecs
						long then = System.currentTimeMillis();
						while ((!sdone || !rdone) && System.currentTimeMillis() - then < timeout) {
							if (!rdone) {
								for (int i = 0; i < replicate; i++) {
									if (rreq[i].testStatus() != null) {
										if (i > 0) {
											int msize = rbuf[igroup + i*k].asIntBuffer().get(0);
											rbuf[igroup].put(rbuf[igroup + i*k].array(), 0, msize);
										}
										rreq[i] = null;
										rdone = true;
										break;
									}
								}
							}
							if (!sdone) {
								for (int i = 0; i < replicate; i++) {
									if (sreq[i].testStatus() != null) {
										sreq[i] = null;
										sdone = true;
										break;
									}
								}
							}
							Thread.sleep(1);
            }  
						for (int i = 0; i < replicate; i++) {
							if (sreq[i] != null && sreq[i].testStatus() == null) sreq[i].cancel();
							if (rreq[i] != null && rreq[i].testStatus() == null) rreq[i].cancel();
						}

						if (!rdone || !sdone) {
							return false;
						} 
						sbuf[igroup].rewind();
						rbuf[igroup].rewind();
					} catch (Exception e) {
						throw new RuntimeException("Exception in sendrecv "+e);
					} */
		}	
		
		public class SockWriter implements Runnable {
			int dest;
			Msg msg;

			public SockWriter(int dest0, Msg msg0) {
				msg = msg0;
				dest = dest0;
			}

			public void run() {
				try {
					Socket socket = new Socket();
					socket.connect(new InetSocketAddress(machineIP[dest], sockBase + dest), sendTimeout);
					if (socket.isConnected()) {
						DataOutputStream ostr = new DataOutputStream(socket.getOutputStream());
						ostr.writeInt(msg.size);
						ostr.writeInt(msg.sender);
						ostr.writeInt(msg.tag);
						ostr.write(msg.buf, 0, msg.size*4);		
						socket.close();
					}
				}	catch (SocketTimeoutException e) {
					// No need to do anything
				} catch (ConnectException e) {
					// Server may have been killed - OK
				}	catch (Exception e) {
					throw new RuntimeException("Problem writing socket "+e);
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
			
			public void stop() {
				try {
					stop = true;
					ss.close();
				} catch (Exception e) {}			
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
	
	public AllReduceX(int M) {
		simNetwork = new Machine[M];
	}
	
	public void stop() {
		if (simNetwork != null) {
			for (int i = 0; i < simNetwork.length; i++) simNetwork[i].stop();
		}
	}

}
