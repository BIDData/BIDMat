
package edu.berkeley.bid.comm;

import java.util.List;
import java.util.LinkedList;
import edu.berkeley.bid.UTILS;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

class AllReduce {	
	
	class Machine {
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
		LinkedList<Msg> [] messages;                             // Message queue for the simulation
		boolean doSim = false;
		ExecutorService executor;

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
				layers[i] = new Layer(k, cumk, left, right, imachine);
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
				messages = new LinkedList[M];
				for (int i = 0; i < M; i++) {
					messages[i] = new LinkedList<Msg>();
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
			int posInMyGroup;                                        // Position in this machines group
			int [] outNbr;                                           // Machines we talk to 
			int [] inNbr;                                            // Machines we listen to
			IVec partBoundaries;                                     // Partition boundaries
			IVec [] downMaps;                                        // Maps to indices below for down indices
			IVec [] upMaps;                                          // Maps to indices below for up indices
			int downn;                                               // Size of the down master list
			int upn;                                                 // Size of the up master list
			int [] dPartInds;
			int [] uPartInds;

			public Layer(int k0, int cumk, int left0, int right0, int imachine) {
				k = k0;
				int i;
				left = left0;
				right = right0;
				partBoundaries = new IVec(k);
				inNbr = new int [k];
				outNbr = new int [k];
				dPartInds = new int[k+1];
				uPartInds = new int[k+1];
				int ioff = imachine % (cumk * k);
				int ibase = imachine - ioff;
				posInMyGroup = ioff / cumk;
				for (i = 0; i < k; i++) {
					partBoundaries.data[i] = left + (int)(((long)(right - left)) * (i+1) / k);
					outNbr[i] = ibase + (ioff + i * cumk) % (cumk * k);
					inNbr[i] = ibase + (ioff + (k - i) * cumk) % (cumk * k);
				}		
			}		
			
			class ConfigThread implements Runnable {
				IVec [] downp;
				IVec [] upp;
				IVec [] dtree;
				IVec [] utree;
				int i;
				CountDownLatch latch;
				
				public ConfigThread(IVec [] downp0, IVec [] upp0, IVec [] dtree0, IVec [] utree0, int i0, CountDownLatch latch0) {
					downp = downp0;
					upp = upp0;
					dtree = dtree0;
					utree = utree0;
					i = i0;
					latch = latch0;
				}

				public void run() {
					int [] segments = new int[2];
					int [] sbuf = sendbuf[i];
					int [] rbuf = recbuf[i];
					dPartInds[i+1] = dPartInds[i] + downp[i].data.length;
					uPartInds[i+1] = uPartInds[i] + upp[i].data.length;					
					System.arraycopy(downp[i].data, 0, sbuf, 2, downp[i].size());
					segments[0] = downp[i].size();
					System.arraycopy(upp[i].data, 0, sbuf, segments[0]+2, upp[i].size());
					segments[1] = segments[0] + upp[i].size();
					sbuf[0] = segments[0];
					sbuf[1] = segments[1];

					sendrecv(sbuf, segments[1]+2, outNbr[i], rbuf, rbuf.length, inNbr[i]);

					IVec downout = new IVec(rbuf[0]);
					IVec upout =   new IVec(rbuf[1]-rbuf[0]);
					System.arraycopy(rbuf, 2, downout.data, 0, rbuf[0]);
					System.arraycopy(rbuf, 2 + rbuf[0], upout.data, 0, rbuf[1]-rbuf[0]);
					IVec.checkTree(dtree, downout, i, k);
					IVec.checkTree(utree, upout, i, k);	
					latch.countDown();
				}
			}

			public void config(IVec downi, IVec upi, IVec [] outputs) {
				IVec [] downp = IVec.partition(downi, partBoundaries);
				IVec [] upp = IVec.partition(upi, partBoundaries);
				IVec [] dtree = new IVec[2*k-1];
				IVec [] utree = new IVec[2*k-1];
				dPartInds[0] = 0;
				uPartInds[0] = 0;
				
				CountDownLatch latch = new CountDownLatch(k);
				for (int i = 0; i < k; i++) {
					executor.execute(new ConfigThread(downp, upp, dtree, utree, i, latch));
				}
				try {	latch.await(); } catch (InterruptedException e) {}
				IVec dmaster = dtree[0];
				downn = dmaster.size();
				IVec umaster = utree[0];
				upn = umaster.size();
				for (int i = 0; i < k; i++) {
					downMaps[i] = IVec.mapInds(downp[i], dmaster);
					upMaps[i] = IVec.mapInds(upp[i], umaster);
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
					
					sendrecv(sbuf, msize+1, outNbr[i], rbuf, rbuf.length, inNbr[i]);
					
					Vec res = new Vec(rbuf[0]);
					UTILS.memcpyif(res.size(), rbuf, 1, res.data, 0);
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
					
					sendrecv(sbuf, up.size()+1, outNbr[i], rbuf, rbuf.length, inNbr[i]);
					
					int msize = uPartInds[i+1] - uPartInds[i];
					UTILS.memcpyif(msize, rbuf, 1, newv.data, uPartInds[i]);						
				}
			}
			
			public Vec reduceUp(Vec upv) {
				Vec newv = new Vec(upn);
				CountDownLatch latch = new CountDownLatch(k);
				for (int i = 0; i < k; i++) {
					executor.execute(new ReduceUpThread(newv, upv, i, latch));
				}
				try { latch.await(); } catch (InterruptedException e) {}
				return newv;
			}
		}

		public void sendrecv(int [] sbuf, int sendn, int outi, int [] rbuf, int recn, int ini) {
			if (doSim) {
				synchronized (simNetwork[outi].messages[imachine]) {
					simNetwork[outi].messages[imachine].add(new Msg(sbuf, sendn, imachine, outi));
					simNetwork[outi].messages[imachine].notify();
				}
				synchronized (messages[ini]) {
					while (messages[ini].size() == 0) {
						try {
							messages[ini].wait();
						} catch (InterruptedException e) {					
						}
					}
					Msg msg = messages[ini].removeFirst();
					System.arraycopy(msg.inbuf, 0, rbuf, 0, sendn);
				}
			} else {
//				MPI.COMM_WORLD.Sendrecv(sbuf, 0, sendn, MPI.INT, outi, 0, rbuf, 0, recn, MPI.INT, ini, 0);
			}
			
		}
	}
	
	public class Msg {
		int [] inbuf;
		int insize;
		int sender;
		int receiver;
		
		public Msg(int [] inbuf0, int insize0, int sender0, int receiver0) {
			inbuf = inbuf0;
			insize = insize0;
			sender = sender0;
			receiver = receiver0;
		}
	}
	
	Machine [] simNetwork;
	
	public static void runSim() {
		
	}

}