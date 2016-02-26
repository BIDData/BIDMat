package edu.berkeley.bid.comm;


import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.FloatBuffer;
import java.nio.LongBuffer;
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
		int ix;
		int iy;
		int repno;
		int tag;

		public ConfigThread(LVec [] top0, LVec [] fromp0,	 int ix0, int iy0, int tag0) {
			top = top0;
			fromp = fromp0;					
			ix = ix0;
			iy = iy0;
			tag = tag0;
		}

		public void run() {
			sendbuf[ix].clear();
			recbuf[ix].clear();
			LongBuffer lsbuf = sendbuf[ix].asLongBuffer();
			LongBuffer lrbuf = recbuf[ix].asLongBuffer();	
			int seg1 = top[ix].size();
			int seg2 = fromp[ix].size();
			lsbuf.put(seg1);
			lsbuf.put(seg2);
			lsbuf.put(top[ix].data, 0, seg1);
			lsbuf.put(fromp[ix].data, 0, seg2);

			if (trace > 1) log(String.format("Config layer %d machine %d sent msg to %d, from %d, sizes %d %d\n", depth, imachine, outNbr[ix],  outNbr[iy],  lsbuf.get(0), lsbuf.get(1)));
			machine.sendrecv(ix, sendbuf[ix], 2*(seg1+seg2+2), outNbr[ix], recbuf[ix], lrbuf.capacity(), outNbr[iy], depth*3 + tag);
			seg1 = (int)lrbuf.get();
			seg2 = (int)lrbuf.get();
			if (trace > 1) log(String.format("Config layer %d machine %d got msg from %d, sizes %d %d\n", depth, imachine, inNbr[iy], seg1, seg2));

			LVec toout = new LVec(seg1);
			LVec upout =   new LVec(seg2);
			lrbuf.get(toout.data, 0, seg1);
			lrbuf.get(upout.data, 0, seg2);	
			top[ix] = toout;
			fromp[ix] = upout;
		}
	}

	public void config(LVec toi, LVec fromi, LVec [] outputs, int round) {
		config_pre(toi, fromi, "Config ");

		Future<?> [] futures = new Future<?>[k];
		Future<?> timeoutf = executor.submit(new TimeoutThread(machine, machine.configTimeout, futures, "config"));
		for (int i = 0; i < k; i++) {
			int ix = (i + posInMyGroup) % k;                   
			int iy = (posInMyGroup - i + k) % k;                     
			futures[i] = executor.submit(new ConfigThread(toindsparts, fromindsparts, ix, iy, 3*D*round), null);
		};
		for (int i = 0; i < k; i++) {
			try {
				futures[i].get();
			} catch (Exception e) {
				if (trace > 0) {
					ByteArrayOutputStream baos = new ByteArrayOutputStream();
					PrintStream ps = new PrintStream(baos);
					e.printStackTrace(ps);
					log(baos.toString());
					ps.close();
				}
			}
		}	
		timeoutf.cancel(true);
		
		config_post(futures, toi, fromi, outputs, "Config ");
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
		if (trace > 1) log(String.format(prefix + "machine %d layer %d, toparts %d\n", imachine, depth, toindsparts[0].size()));
	}
	
	public void config_post(Future<?> [] futures, LVec toi, LVec fromi, LVec [] outputs, String prefix) {
		LVec [] totree = new LVec[k];
		LVec [] fromtree = new LVec[k];
		for (int i = 0 ; i < k; i++) {
			if (futures[i].isDone()) {
				LVec.treeAdd(toindsparts[i], totree);
				LVec.treeAdd(fromindsparts[i], fromtree);
			}
		}
		LVec tomaster = LVec.treeFlush(totree);
		LVec frommaster = LVec.treeFlush(fromtree);
		ton = tomaster.size();
		fromn = fromi.size();
		for (int i = 0; i < k; i++) {
			if (futures[i].isDone()) {
				toMaps[i] = LVec.mapInds(toindsparts[i], tomaster);
				fromMaps[i] = LVec.mapInds(fromindsparts[i], frommaster);
				if (trace > 1) {
					log(String.format(prefix + "machine %d layer %d dmap(%d) size %d\n", imachine, depth, i, toMaps[i].size()));
					log(String.format(prefix + "machine %d layer %d umap(%d) size %d\n", imachine, depth, i, fromMaps[i].size()));
				}
			} else {
				toMaps[i] = null;
				fromMaps[i] = null;
			}
		}
		outputs[0] = tomaster;
		outputs[1] = frommaster;
	}
	
	class ConfigReduceThread implements Runnable {
		LVec [] top;
		LVec [] fromp;
		Vec downv;
		Vec [] res;
		int stride;
		int ix;
		int iy;
		int repno;
		int tag;

		public ConfigReduceThread(LVec [] top0, LVec [] fromp0,	Vec downv0, Vec [] res0, int ix0, int iy0, int stride0, int tag0) {
			top = top0;
			fromp = fromp0;					
			downv = downv0;
			res = res0;
			ix = ix0;
			iy = iy0;
			stride = stride0;
			tag = tag0;
		}

		public void run() {
			sendbuf[ix].clear();
			recbuf[ix].clear();
			LongBuffer isbuf = sendbuf[ix].asLongBuffer();
			LongBuffer irbuf = recbuf[ix].asLongBuffer();	
			FloatBuffer fsbuf = sendbuf[ix].asFloatBuffer();
			FloatBuffer frbuf = recbuf[ix].asFloatBuffer();
			int seg1 = top[ix].size();
			int seg2 = fromp[ix].size();
			int msize = downv.size();
			isbuf.put(seg1);
			isbuf.put(seg2);
			isbuf.put(msize);
			isbuf.put(top[ix].data, 0, seg1);
			isbuf.put(fromp[ix].data, 0, seg2);
			fsbuf.position(2*(3+seg1+seg2));
			fsbuf.put(downv.data, 0, msize);

			if (trace > 1) log(String.format("ConfigReduce layer %d machine %d sent msg to %d, from %d, sizes %d %d %d\n", depth, imachine, outNbr[ix],  outNbr[iy],  isbuf.get(0), isbuf.get(1), isbuf.get(2)));
			machine.sendrecv(ix, sendbuf[ix], 2*(seg1+seg2+3)+msize, outNbr[ix], recbuf[ix], frbuf.capacity(), outNbr[iy], depth*3 + tag);
			seg1 = (int)irbuf.get();
			seg2 = (int)irbuf.get();
			msize = (int)irbuf.get();
			if (trace > 1) log(String.format("ConfigReduce layer %d machine %d got msg from %d, sizes %d %d %d\n", depth, imachine, outNbr[iy], seg1, seg2, msize));

			LVec toout = new LVec(seg1);
			LVec fromout = new LVec(seg2);
			irbuf.get(toout.data, 0, seg1);
			irbuf.get(fromout.data, 0, seg2);
			top[ix] = toout;
			fromp[ix] = fromout;
			
			res[ix] = new Vec(msize);
			frbuf.position(2*(3+seg1+seg2));
			frbuf.get(res[ix].data, 0, msize);
		}
	}
	
	public Vec configReduce(LVec toi, LVec fromi, LVec [] outputs, Vec downv, int stride, int round) {
		
		config_pre(toi, fromi, "ConfigReduce ");
		
		Vec [] vparts = partitioner.partition(downv, topart, topartsizes, stride);
		Vec [] res = new Vec[k];

		Future<?> [] futures = new Future<?>[k];
		Future<?> timeoutf = executor.submit(new TimeoutThread(machine, machine.configTimeout, futures, "configReduce"));
		for (int i = 0; i < k; i++) {
			int ix = (i + posInMyGroup) % k;
			int iy = (posInMyGroup - i + k) % k;
			futures[i] = executor.submit(new ConfigReduceThread(toindsparts, fromindsparts, vparts[ix], res, ix, iy, stride, round), null);
		};
		for (int i = 0; i < k; i++) {
			try {
				futures[i].get();
			} catch (Exception e) {
				if (trace > 0) {
					ByteArrayOutputStream baos = new ByteArrayOutputStream();
					PrintStream ps = new PrintStream(baos);
					e.printStackTrace(ps);
					log(baos.toString());
					ps.close();
				}
			}
		}	
		timeoutf.cancel(true);
		
		config_post(futures, toi, fromi, outputs, "ConfigReduce ");
		
		Vec newv = new Vec(ton*stride);
		for (int i = 0; i < k; i++) {
			res[i].addTo(newv, toMaps[i], stride);
		}
		return newv;
	}
	
	private void log(String msg) {
		machine.log(msg);
	}
}