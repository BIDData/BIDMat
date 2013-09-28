
package edu.berkeley.bid.comm;

import java.util.List;
import java.util.LinkedList;

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
		IVec finalMap;
		LinkedList<Msg> [] messages;
		boolean doSim = false;

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
			int k;        																					 // size of this group
			int left;                                                // left boundary of its indices
			int right;                                               // right boundary of its indices
			int posInMyGroup;                                        // position in this machines group
			int [] outNbr;                                           // machines we talk to 
			int [] inNbr;                                            // machines we listen to
			IVec partBoundaries;                                     // partition boundaries
			IVec [] downMaps;                                        // maps to indices below for down indices
			IVec [] upMaps;                                          // maps to indices below for up indices
			int [] dPartInds;
			int [] uPartInds;

			public Layer(int k0, int cumk, int left0, int right0, int imachine) {
				k = k0;
				int i;
				left = left0;
				right = right0;
				partBoundaries = new IVec(k);
				for (i = 0; i < k; i++) {
					partBoundaries.data[i] = left + (int)(((long)(right - left)) * (i+1) / k);				
				}
				inNbr = new int [k];
				outNbr = new int [k];
				dPartInds = new int[k+1];
				uPartInds = new int[k+1];
				int ioff = imachine % (cumk * k);
				int ibase = imachine - ioff;
				posInMyGroup = ioff / cumk;
				for (i = 0; i < k; i++) {
					outNbr[i] = ibase + (ioff + i * cumk) % (cumk * k);
					inNbr[i] = ibase + (ioff + (k - i) * cumk) % (cumk * k);
				}		
			}		

			public void config(IVec downi, IVec upi, IVec [] outputs) {
				IVec [] downp = IVec.partition(downi, partBoundaries);
				IVec [] upp = IVec.partition(upi, partBoundaries);
				IVec [] dtree = new IVec[2*k-1];
				IVec [] utree = new IVec[2*k-1];
				dPartInds[0] = 0;
				uPartInds[0] = 0;
				for (int i = 0; i < k; i++) {
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
				}
				IVec dmaster = dtree[0];
				IVec umaster = utree[0];
				for (int i = 0; i < k; i++) {
					downMaps[i] = IVec.mapInds(downp[i], dmaster);
					upMaps[i] = IVec.mapInds(upp[i], umaster);
				}
				outputs[0] = dmaster;
				outputs[1] = umaster;
			}
			
			public Vec reduceDown(Vec downv) {
				return downv;
			}
			
			public Vec reduceUp(Vec upv) {
				return upv;
			}
		}

		public void sendrecv(int [] sbuf, int sendn, int outi, int [] rbuf, int recn, int ini) {
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

}