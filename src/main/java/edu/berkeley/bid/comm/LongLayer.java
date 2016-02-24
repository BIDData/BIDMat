package edu.berkeley.bid.comm;


import java.nio.FloatBuffer;
import java.nio.LongBuffer;
import java.util.Arrays;
import java.util.concurrent.Future;

/**
 * 
 * @author jfc
 * A Layer with Long indices
 *
 */
public class LongLayer extends Layer {

	/* Layer Configuration Variables */

	public LVec [] toindsparts;                                          // Partition of to indices 
	public LVec [] fromindsparts;                                        // Partition of from indices

	public LongLayer(Machine machine0, int k0, int cumk0, int cumPos0, int posInMyGroup0, int depth0) {
		super(machine0, k0, cumk0, cumPos0, posInMyGroup0, depth0);
	}		

	class ConfigThread implements Runnable {
		LVec [] top;
		LVec [] fromp;
		LVec [] totree;
		LVec [] fromtree;
		int i;
		int repno;

		public ConfigThread(LVec [] top0, LVec [] fromp0,	LVec [] totree0, LVec [] fromtree0, int i0) {
			top = top0;
			fromp = fromp0;					
			totree = totree0;
			fromtree = fromtree0;
			i = i0;
		}

		public void run() {
			sendbuf[i].clear();
			recbuf[i].clear();
			LongBuffer sbuf = sendbuf[i].asLongBuffer();
			LongBuffer rbuf = recbuf[i].asLongBuffer();	
			int seg1 = top[i].size();
			int seg2 = fromp[i].size();
			sbuf.put(seg1);
			sbuf.put(seg2);
			sbuf.put(top[i].data, 0, seg1);
			sbuf.put(fromp[i].data, 0, seg2);

			if (trace > 1 && network != null) {
				synchronized (network) {
					System.out.format("Config layer %d machine %d sent msg to %d, from %d, sizes %d %d\n", depth, imachine, outNbr[i],  inNbr[i],  sbuf.get(0), sbuf.get(1));
				}
			}
			machine.sendrecv(i, sendbuf, 2*(seg1+seg2+2), outNbr[i], recbuf, rbuf.capacity(), inNbr[i], depth*3);
			seg1 = (int)rbuf.get();
			seg2 = (int)rbuf.get();
			if (trace > 1 && network != null) {
				synchronized (network) {
					System.out.format("Config layer %d machine %d got msg from %d, sizes %d %d\n", depth, imachine, inNbr[i], seg1, seg2);
				}
			}

			LVec toout = new LVec(seg1);
			LVec upout =   new LVec(seg2);
			rbuf.get(toout.data, 0, seg1);
			rbuf.get(upout.data, 0, seg2);
			LVec.checkTree(totree, toout, i, k);
			LVec.checkTree(fromtree, upout, i, k);	
			top[i] = toout;
			fromp[i] = upout;
		}
	}

	public void config(LVec toi, LVec fromi, LVec [] outputs) {
		LVec [] totree = new LVec[2*k-1];
		LVec [] fromtree = new LVec[2*k-1];
		config_pre(toi, fromi, "Config ");

		Future<?> [] futures = new Future<?>[k];
		Future<?> timeoutf = executor.submit(new TimeoutThread(machine.configTimeout, futures), null);
		for (int i = 0; i < k; i++) {
			int ix = (i + posInMyGroup) % k;                               // Try to stagger the traffic
			futures[i] = executor.submit(new ConfigThread(toindsparts, fromindsparts, totree, fromtree, ix), null);
		};
		for (int i = 0; i < k; i++) {
			try {
				futures[i].get();
			} catch (Exception e) {
				if (trace > 0 && network != null) {
					synchronized (network) {
						e.printStackTrace();
					}
				}
			}
		}	
		timeoutf.cancel(true);
		
		config_post(toi, fromi, outputs, totree, fromtree, "Config ");
	}
	
	public void config_pre(LVec toi, LVec fromi, String prefix) {
		topart = partitioner.part(toi);                                   // Compute partition indices
		frompart = partitioner.part(fromi);
		interleave = new IVec[k];
		toindsparts = partitioner.partition(toi, topart, null);          // Actually partition the indices into k LVec parts
		fromindsparts = partitioner.partition(fromi, frompart, interleave);
		topartsizes = new IVec(k);
		frompartsizes = new IVec(k);
		fromparts = new Vec[k];
		for (int i = 0; i < k; i++) {
			topartsizes.data[i] = toindsparts[i].length;
			frompartsizes.data[i] = fromindsparts[i].length;
		}
		if (trace > 0 && network != null) {
			synchronized (network) {
				System.out.format(prefix + "machine %d layer %d, toparts (%d", imachine, depth, toindsparts[0].size());
			}
		}
	}
	
	public void config_post(LVec toi, LVec fromi, LVec [] outputs, LVec [] totree, LVec [] fromtree, String prefix) {
		LVec dmaster = totree[0];
		Arrays.fill(totree, null);
		ton = dmaster.size();
		LVec umaster = fromtree[0];
		Arrays.fill(fromtree, null);
		fromn = fromi.size();
		for (int i = 0; i < k; i++) {
			toMaps[i] = LVec.mapInds(toindsparts[i], dmaster);
			fromMaps[i] = LVec.mapInds(fromindsparts[i], umaster);
			if (trace > 0 && network != null) {
				synchronized (network) { 
					System.out.format(prefix + "machine %d layer %d dmap(%d) size %d\n", imachine, depth, i, toMaps[i].size());
					System.out.format(prefix + "machine %d layer %d umap(%d) size %d\n", imachine, depth, i, fromMaps[i].size());
				}
			}
		}
		outputs[0] = dmaster;
		outputs[1] = umaster;
	}
	
	class ConfigReduceThread implements Runnable {
		LVec [] top;
		LVec [] fromp;
		LVec [] totree;
		LVec [] fromtree;
		Vec downv;
		Vec [] res;
		int stride;
		int i;
		int repno;

		public ConfigReduceThread(LVec [] top0, LVec [] fromp0,	LVec [] totree0, LVec [] fromtree0, Vec downv0, Vec [] res0, int i0, int stride0) {
			top = top0;
			fromp = fromp0;					
			totree = totree0;
			fromtree = fromtree0;
			downv = downv0;
			res = res0;
			i = i0;
			stride = stride0;
		}

		public void run() {
			sendbuf[i].clear();
			recbuf[i].clear();
			LongBuffer isbuf = sendbuf[i].asLongBuffer();
			LongBuffer irbuf = recbuf[i].asLongBuffer();	
			FloatBuffer sbuf = sendbuf[i].asFloatBuffer();
			FloatBuffer rbuf = recbuf[i].asFloatBuffer();
			int seg1 = top[i].size();
			int seg2 = fromp[i].size();
			int msize = downv.size();
			isbuf.put(seg1);
			isbuf.put(seg2);
			isbuf.put(msize);
			isbuf.put(top[i].data, 0, seg1);
			isbuf.put(fromp[i].data, 0, seg2);
			sbuf.position(2*(3+seg1+seg2));
			sbuf.put(downv.data, 0, msize);

			if (trace > 1 && network != null) {
				synchronized (network) {
					System.out.format("ConfigReduce layer %d machine %d sent msg to %d, from %d, sizes %d %d %d\n", depth, imachine, outNbr[i],  inNbr[i],  isbuf.get(0), isbuf.get(1), isbuf.get(2));
				}
			}
			machine.sendrecv(i, sendbuf, 2*(seg1+seg2+3)+msize, outNbr[i], recbuf, rbuf.capacity(), inNbr[i], depth*3);
			seg1 = (int)irbuf.get();
			seg2 = (int)irbuf.get();
			msize = (int)irbuf.get();
			if (trace > 1 && network != null) {
				synchronized (network) {
					System.out.format("ConfigReduce layer %d machine %d got msg from %d, sizes %d %d %d\n", depth, imachine, inNbr[i], seg1, seg2, msize);
				}
			}

			LVec toout = new LVec(seg1);
			LVec fromout = new LVec(seg2);
			irbuf.get(toout.data, 0, seg1);
			irbuf.get(fromout.data, 0, seg2);
			LVec.checkTree(totree, toout, i, k);
			LVec.checkTree(fromtree, fromout, i, k);	
			top[i] = toout;
			fromp[i] = fromout;
			
			res[i] = new Vec(msize);
			rbuf.position(2*(3+seg1+seg2));
			rbuf.get(res[i].data, 0, msize);
		}
	}
	
	public Vec configReduce(LVec toi, LVec fromi, LVec [] outputs, Vec downv, int stride) {
		
		config_pre(toi, fromi, "ConfigReduce ");
		
		LVec [] totree = new LVec[2*k-1];
		LVec [] fromtree = new LVec[2*k-1];
		Vec [] vparts = partitioner.partition(downv, topart, topartsizes, stride);
		Vec [] res = new Vec[k];

		Future<?> [] futures = new Future<?>[k];
		Future<?> timeoutf = executor.submit(new TimeoutThread(machine.configTimeout, futures), null);
		for (int i = 0; i < k; i++) {
			int ix = (i + posInMyGroup) % k;                               // Try to stagger the traffic
			futures[i] = executor.submit(new ConfigReduceThread(toindsparts, fromindsparts, totree, fromtree, vparts[ix], res, ix, stride), null);
		};
		for (int i = 0; i < k; i++) {
			try {
				futures[i].get();
			} catch (Exception e) {
				if (trace > 0 && network != null) {
					synchronized (network) {
						e.printStackTrace();
					}
				}
			}
		}	
		timeoutf.cancel(true);
		
		config_post(toi, fromi, outputs, totree, fromtree, "ConfigReduce ");
		
		Vec newv = new Vec(ton*stride);
		for (int i = 0; i < k; i++) {
			res[i].addTo(newv, toMaps[i], stride);
		}
		return newv;
	}
}