
package edu.berkeley.bid.comm;

import java.util.LinkedList;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.Arrays;
import java.nio.*;
//import mpi.*;

public class AllReduce {
	
	public class Machine {
		/* Machine Configuration Variables */	
		int N;                                                     // Number of features
		int D;                                                     // Depth of the network
		int M;                                                     // Number of Machines
		int imachine;                                              // My identity
		int [] allks;                                              // k values
		
		Layer [] layers;                                           // All the layers
		ByteBuffer [] sendbuf;                                     // buffers, one for each destination in a group
		ByteBuffer [] recbuf;
		IVec finalMap;                                             // Map to down from down --> up at layer D-1
		LinkedList<Msg> [][] messages;                             // Message queue for the simulation
		boolean doSim = true;
		ExecutorService executor;
		int trace = 0;                                             // 0: no trace, 1: high-level, 2: everything

		public Machine(int N0, int [] allks0, int imachine0, int M0, int bufsize, boolean doSim0, int trace0) {
			N = N0;
			M = M0;
			doSim = doSim0;
			if (doSim) {
				imachine = imachine0;
			} else {
/*				try {
					String [] args = {""};
					MPI.InitThread(args, MPI.THREAD_MULTIPLE);
					imachine = MPI.COMM_WORLD.getRank();
				} catch (MPIException e) {
					throw new RuntimeException("Couldnt init MPI "+e);
				} */
			}
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
				sendbuf[i] = ByteBuffer.wrap(new byte[bufsize*4]);
				recbuf[i] = ByteBuffer.wrap(new byte[bufsize*4]);
			}
			if (doSim) {
				messages = new LinkedList[M][];
				for (int i = 0; i < M; i++) {
					messages[i] = new LinkedList[3*D];
					for (int j = 0; j < 3*D; j++) {
						messages[i][j] = new LinkedList<Msg>();
					}
				}
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
				synchronized (AllReduce.this) {
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
				int ibase = imachine - posInMyGroup * cumk;
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
						synchronized (AllReduce.this) {
							System.out.format("config layer %d machine %d sent msg to %d, from %d, sizes %d %d\n", depth, imachine, outNbr[i],  inNbr[i],  sbuf.get(0), sbuf.get(1));
						}
					}
					sendrecv(sbuf, seg2+2, outNbr[i], rbuf, rbuf.capacity(), inNbr[i], depth*3);
					seg1 = rbuf.get();
					seg2 = rbuf.get();
					if (trace > 1) {
						synchronized (AllReduce.this) {
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
					synchronized (AllReduce.this) {
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
  					synchronized (AllReduce.this) { 
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
						synchronized (AllReduce.this) {
							System.out.format("reduce layer %d machine %d sent msg to %d, from %d, size %d\n", depth, imachine, outNbr[i],  inNbr[i],  msize);
						}
					}
					sendrecv(isbuf, msize+1, outNbr[i], irbuf, rbuf.capacity(), inNbr[i], depth*3+1);
					msize = irbuf.get();
					if (trace > 1) {
						synchronized (AllReduce.this) {
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
						synchronized (AllReduce.this) {
							System.out.format("reduce up layer %d machine %d sent msg to %d, from %d, size %d\n", depth, imachine, outNbr[i],  inNbr[i],  msize);
						}
					}
					sendrecv(isbuf, msize+1, inNbr[i], irbuf, irbuf.capacity(), outNbr[i], depth*3+2);
					msize = irbuf.get();
					if (trace > 1) {
						synchronized (AllReduce.this) {
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

		public boolean sendrecv(IntBuffer sbuf, int sendn, int outi, IntBuffer rbuf, int recn, int ini, int tag) {
			if (imachine == outi) {
				Msg a = new Msg(sbuf, sendn, imachine, outi);
				rbuf.clear();
				rbuf.put(a.buf, 0, sendn);
				rbuf.rewind();
				return true;				
			} else {
				if (doSim) {
					synchronized (simNetwork[outi].messages[imachine][tag]) {
						simNetwork[outi].messages[imachine][tag].add(new Msg(sbuf, sendn, imachine, outi));
						simNetwork[outi].messages[imachine][tag].notify();
					}
					synchronized (messages[ini][tag]) {
						while (messages[ini][tag].size() == 0) {
							try {
								messages[ini][tag].wait();
							} catch (InterruptedException e) {}
						}
						Msg msg = messages[ini][tag].removeFirst();
						rbuf.clear();
						rbuf.put(msg.buf, 0, msg.size);
						rbuf.rewind();
					}
					return true;
				} else {
/*					try {
            sbuf.rewind();
            rbuf.clear();
						MPI.COMM_WORLD.sendRecv(sbuf, sendn, MPI.INT, outi, tag, rbuf, recn, MPI.INT, ini, tag); 
						sbuf.rewind();
						rbuf.rewind();
					} catch (MPIException e) {
						throw new RuntimeException("Exception in sendrecv "+e);
					} */
/*					try {	
            sbuf.rewind();
            rbuf.clear();	
						Request sreq = MPI.COMM_WORLD.iSend(sbuf, sendn, MPI.INT, outi, tag);
						Request rreq = MPI.COMM_WORLD.iRecv(rbuf, recn, MPI.INT, ini, tag);
						Status rdone = null;
						Status sdone = null;
						long timeout = 1000;   // Wait this many msecs
						long then = System.currentTimeMillis();
						while ((sdone == null || rdone == null) && System.currentTimeMillis() - then < timeout) {
							if (rdone == null) rdone = rreq.testStatus();
							if (sdone == null) sdone = sreq.testStatus();
							Thread.sleep(1);
						}
						if (rdone == null) rreq.cancel();
						if (sdone == null) sreq.cancel();
						if (rdone == null || sdone == null) {
							return false;
						} 
						sbuf.rewind();
						rbuf.rewind();
					} catch (Exception e) {
						throw new RuntimeException("Exception in sendrecv "+e);
					} */
				}
				return true;
			}		
		}
	}
	
	public class Msg {
		int [] buf;
		int size;
		int sender;
		int receiver;
		
		public Msg(IntBuffer inbuf, int size0, int sender0, int receiver0) {
			buf = new int[size0];
			inbuf.rewind();
			inbuf.get(buf, 0, size0);
			inbuf.rewind();
			size = size0;
			sender = sender0;
			receiver = receiver0;
		}
	}
	
	public Machine [] simNetwork;
	
	public AllReduce(int M) {
		simNetwork = new Machine[M];
	}

}
