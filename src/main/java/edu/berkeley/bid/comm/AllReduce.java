
package edu.berkeley.bid.comm;

import java.util.List;
import java.util.LinkedList;
import edu.berkeley.bid.UTILS;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class AllReduce {
	
	public class Machine {
		/* Machine Configuration Variables */	
		int N;                                                   // Number of features
		int D;                                                   // Depth of the network
		int M;                                                   // Number of Machines
		int imachine;                                            // My identity
		int [] allks;                                            // k values
		
		Layer [] layers;                                         // All the layers
		int [][] sendbuf;                                        // buffers, one for each destination in a group
		int [][] recbuf;
		IVec finalMap;                                           // Map to down from down --> up at layer D-1
		LinkedList<Msg> [][] messages;                           // Message queue for the simulation
		boolean doSim = true;
		ExecutorService executor;
		boolean trace = false;

		public Machine(int N0, int [] allks0, int imachine0, int M0, int bufsize, boolean doSim0) {
			N = N0;
			M = M0;
			imachine = imachine0;
			allks = allks0;
			D = allks.length;
			doSim = doSim0;
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
			executor = Executors.newFixedThreadPool(maxk);
			sendbuf = new int[maxk][];
			recbuf = new int[maxk][];
			for (int i = 0; i < maxk; i++) {
				sendbuf[i] = new int[bufsize];
				recbuf[i] = new int[bufsize];
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
				int ioff = imachine % cumk;
				int ibase = imachine/(cumk*k);
				posInMyGroup = ioff / cumk;
				for (i = 0; i < k; i++) {
					partBoundaries.data[i] = left + (int)(((long)(right - left)) * (i+1) / k);
					outNbr[i] = ibase + (ioff + i * cumk) % (cumk * k);
					inNbr[i] = outNbr[i];
//					inNbr[i] = ibase + (ioff + (k - i) * cumk) % (cumk * k);
				}		
				downMaps = new IVec[k];
				upMaps = new IVec[k];
			}		
			
			class ConfigThread implements Runnable {
				IVec [] downp;
				IVec [] upp;
				IVec [] newdownp;
				IVec [] newupp; 
				IVec [] dtree;
				IVec [] utree;
				int i;
				CountDownLatch latch;
				
				public ConfigThread(IVec [] downp0, IVec [] upp0, IVec [] newdownp0, IVec [] newupp0, 
						IVec [] dtree0, IVec [] utree0, int i0, CountDownLatch latch0) {
					downp = downp0;
					upp = upp0;					
					newdownp = newdownp0;
					newupp = newupp0;
					dtree = dtree0;
					utree = utree0;
					i = i0;
					latch = latch0;
				}

				public void run() {
					int [] segments = new int[2];
					int [] sbuf = sendbuf[i];
					int [] rbuf = recbuf[i];				
					System.arraycopy(downp[i].data, 0, sbuf, 2, downp[i].size());
					segments[0] = downp[i].size();
					System.arraycopy(upp[i].data, 0, sbuf, segments[0]+2, upp[i].size());
					segments[1] = segments[0] + upp[i].size();
					sbuf[0] = segments[0];
					sbuf[1] = segments[1];

					if (trace) System.out.format("config layer %d machine %d sent msg to %d, from %d, sizes %d %d\n", depth, imachine, outNbr[i],  inNbr[i],  sbuf[0], sbuf[1]);
					sendrecv(sbuf, segments[1]+2, outNbr[i], rbuf, rbuf.length, inNbr[i], depth*3);
					if (trace) System.out.format("config layer %d machine %d got msg from %d, sizes %d %d\n", depth, imachine, inNbr[i], rbuf[0], rbuf[1]);

					IVec downout = new IVec(rbuf[0]);
					IVec upout =   new IVec(rbuf[1]-rbuf[0]);
					System.arraycopy(rbuf, 2, downout.data, 0, rbuf[0]);
					System.arraycopy(rbuf, 2 + rbuf[0], upout.data, 0, rbuf[1]-rbuf[0]);
					IVec.checkTree(dtree, downout, i, k);
					IVec.checkTree(utree, upout, i, k);	
					newdownp[i] = downout;
					newupp[i] = upout;
					latch.countDown();
				}
			}

			public void config(IVec downi, IVec upi, IVec [] outputs) {
				IVec [] downp = IVec.partition(downi, partBoundaries);
				IVec [] upp = IVec.partition(upi, partBoundaries);
				IVec [] newdownp = new IVec[downp.length];
				IVec [] newupp = new IVec[upp.length];
				IVec [] dtree = new IVec[2*k-1];
				IVec [] utree = new IVec[2*k-1];
//				System.out.format("machine %d layer %d, dparts %d %d\n", imachine, depth, downp[0].size(), downp[1].size());
				if (trace) System.out.format("machine %d layer %d, uparts %d %d\n", imachine, depth, upp[0].size(), upp[1].size());
				dPartInds[0] = 0;
				uPartInds[0] = 0;
				
				CountDownLatch latch = new CountDownLatch(k);
				for (int i = 0; i < k; i++) {
					dPartInds[i+1] = dPartInds[i] + downp[i].size();
					uPartInds[i+1] = uPartInds[i] + upp[i].size();
					
					executor.execute(new ConfigThread(downp, upp, newdownp, newupp, dtree, utree, i, latch));
				}
				try {	latch.await(); } catch (InterruptedException e) {}
				IVec dmaster = dtree[0];
				downn = dmaster.size();
				IVec umaster = utree[0];
				upn = upi.size();
				for (int i = 0; i < k; i++) {
					downMaps[i] = IVec.mapInds(newdownp[i], dmaster);
					upMaps[i] = IVec.mapInds(newupp[i], umaster);
//					System.out.format("machine %d dmap(%d) size %d\n", imachine, i, downMaps[i].size());
					if (trace) System.out.format("machine %d umap(%d) size %d\n", imachine, i, upMaps[i].size());
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
					int [] sbuf = sendbuf[i];
					int [] rbuf = recbuf[i];
					int msize = dPartInds[i+1] - dPartInds[i];
					sbuf[0] = msize;
					UTILS.memcpyfi(msize, downv.data, dPartInds[i], sbuf, 1);
					
					if (trace) System.out.format("reduce layer %d machine %d sent msg to %d, from %d, size %d\n", depth, imachine, outNbr[i],  inNbr[i],  sbuf[0]);
					sendrecv(sbuf, msize+1, outNbr[i], rbuf, rbuf.length, inNbr[i], depth*3+1);
					if (trace) System.out.format("reduce layer %d machine %d got msg from %d, size %d\n", depth, imachine, inNbr[i], rbuf[0]);
					
					Vec res = new Vec(rbuf[0]);
					UTILS.memcpyif(rbuf[0], rbuf, 1, res.data, 0);
					res.addTo(newv, downMaps[i]);	
					latch.countDown();
				}
			}
			
			public Vec reduceDown(Vec downv) {
				Vec newv = new Vec(downn);
				CountDownLatch latch = new CountDownLatch(k);
				for (int i = 0; i < k; i++) {
					executor.execute(new ReduceDownThread(newv, downv, i, latch));
				}
				try { latch.await(); } catch (InterruptedException e) {}
				System.out.format("layer %d machine %d reduce down finished\n", depth, imachine);
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
					int [] sbuf = sendbuf[i];
					int [] rbuf = recbuf[i];
					Vec up = upv.mapFrom(upMaps[i]);
					sbuf[0] = up.size();
					UTILS.memcpyfi(up.size(), up.data, 0, sbuf, 1);
					
					if (trace) System.out.format("reduce up layer %d machine %d sent msg to %d, from %d, size %d\n", depth, imachine, outNbr[i],  inNbr[i],  sbuf[0]);
					sendrecv(sbuf, up.size()+1, inNbr[i], rbuf, rbuf.length, outNbr[i], depth*3+2);
					if (trace) System.out.format("reduce up layer %d machine %d got msg from %d, size %d\n", depth, imachine, inNbr[i], rbuf[0]);
					
					int msize = uPartInds[i+1] - uPartInds[i];
					if (uPartInds[i+1] > newv.size()) throw new RuntimeException("ReduceUp index out of range "+uPartInds[i+1]+" "+newv.size());
					UTILS.memcpyif(msize, rbuf, 1, newv.data, uPartInds[i]);	
					latch.countDown();
				}
			}
			
			public Vec reduceUp(Vec upv) {
				Vec newv = new Vec(upn);
				CountDownLatch latch = new CountDownLatch(k);
				for (int i = 0; i < k; i++) {
					executor.execute(new ReduceUpThread(newv, upv, i, latch));
				}
				try { latch.await(); } catch (InterruptedException e) {}
				System.out.format("layer %d machine %d reduce up finished\n", depth, imachine);
				return newv;
			}
		}

		public boolean sendrecv(int [] sbuf, int sendn, int outi, int [] rbuf, int recn, int ini, int tag) {
			if (doSim) {
				synchronized (simNetwork[outi].messages[imachine][tag]) {
//					if (trace) System.out.format("Message sent %d %d %d\n", imachine, sendn, outi);
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
//					if (trace) System.out.format("Message recv %d %d %d\n", imachine, msg.sender, msg.receiver, msg.size);
					System.arraycopy(msg.buf, 0, rbuf, 0, msg.size);
				}
				return true;
			} else {
/*				MPI.COMM_WORLD.Sendrecv(sbuf, 0, sendn, MPI.INT, outi, tag, rbuf, 0, recn, MPI.INT, ini, tag);
			  Request sreq = MPI.COMM_WORLD.ISend(sbuf, 0, sendn, MPI.INT, outi, tag)
				Request rreq = MPI.COMM_WORLD.IRecv(rbuf, 0, recn, MPI.INT, ini, tag)
				Status rdone = null;
				Status sdone = null;
				long timeout = 1000;   // Wait this many msecs
				long then = System.currentTimeMillis();
				while ((sdone == null || rdone == null) && System.currentTimeMillis() - then < timeout) {
					if (rdone == null) rdone = rreq.Test();
					if (sdone == null) sdone = sreq.Test();
					sleep(1);
				} 
				if (rdone == null) rreq.Cancel();
				if (sdone == null) sreq.Cancel();
				if (rdone == null || sdone == null) {
				  return false;
				} else {
				  return true;
				}
				*/
				return true;			
			}		
		}
	}
	
	public class Msg {
		int [] buf;
		int size;
		int sender;
		int receiver;
		
		public Msg(int [] inbuf0, int size0, int sender0, int receiver0) {
			buf = new int[size0];
			System.arraycopy(inbuf0, 0, buf, 0, size0);
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
