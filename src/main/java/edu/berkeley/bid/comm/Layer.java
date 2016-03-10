package edu.berkeley.bid.comm;

import java.nio.ByteBuffer;
import java.nio.FloatBuffer;
import java.nio.IntBuffer;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

public class Layer {

	/* Layer Configuration Variables */

	public final int imachine;                                           // This machines global id number
	public final int k;        																					 // Size of this group
	public final int cumk;                                               // cumulative product of sizes of layers above
	public final int cumPos;                                             // cumulative position in groups
	public final int depth;                                              // How deep is this layer from the top
	public final int posInMyGroup;                                       // Position in this machines group
	public int [] outNbr;                                                // Machines we talk to 
	public int [] inNbr;                                                 // Machines we listen to
	public IVec topart;                                                  // partition number for each index
	public IVec frompart;
	public IVec [] toindsparts;                                          // Partition of to indices 
	public IVec [] fromindsparts;                                        // Partition of from indices
	public IVec topartsizes;                                             // Sizes of the above
	public IVec frompartsizes;
	public IVec [] toMaps;                                               // Maps to indices in layer below for "TO" indices
	public IVec [] fromMaps;                                             // Maps to indices in layer below for "FROM" indices
	public IVec [] interleave;                                           // Maps to indices above for "FROM" indices
	public Vec [] fromparts;
	public int ton;                                                      // Size of the "TO" union list
	public int fromn;                                                    // Size of the "FROM" union list
	public int trace;
	public Partitioner partitioner;
	public Machine machine;
	public Groups groups;
	public int D;
	ExecutorService executor;
	Network network;
	public ByteBuffer [] sendbuf;                                        // buffers, one for each destination in a group
  public ByteBuffer [] recbuf;


	public Layer(Machine machine0, int k0, int cumk0, int cumPos0, int posInMyGroup0, int depth0) {
		machine = machine0;
		k = k0;
		cumk = cumk0;
		cumPos = cumPos0;
		posInMyGroup = posInMyGroup0;
		depth = depth0;
		groups = machine.groups;
		sendbuf = machine.sendbuf;
		recbuf = machine.recbuf;
		network = machine.network;
		trace = machine.trace;
		imachine = machine.imachine;
		executor = machine.executor;
		D = machine.D;
		partitioner = new Partitioner(cumk, cumPos, k);
		inNbr = groups.nodesInGroup(imachine, depth);
		outNbr = groups.nodesInGroup(imachine, depth);		
		toMaps = new IVec[k];
		fromMaps = new IVec[outNbr.length];
	}		

	class ConfigThread implements Runnable {
		IVec [] top;
		IVec [] fromp;
		int ix;
		int iy;
		int repno;
		int tag;

		public ConfigThread(IVec [] top0, IVec [] fromp0,	int ix0, int iy0, int tag0) {
			top = top0;
			fromp = fromp0;					
			ix = ix0;
			iy = iy0;
			tag = tag0;
		}

		public void run() {
			sendbuf[ix].clear();
			recbuf[ix].clear();
			IntBuffer isbuf = sendbuf[ix].asIntBuffer();
			IntBuffer irbuf = recbuf[ix].asIntBuffer();	
			int seg1 = top[ix].size();
			int seg2 = fromp[ix].size();
			isbuf.put(seg1);
			isbuf.put(seg2);
			isbuf.put(top[ix].data, 0, seg1);
			isbuf.put(fromp[ix].data, 0, seg2);

			if (trace > 1) log(String.format("Config layer %d machine %d sent msg to %d, from %d, sizes %d %d\n", depth, imachine, outNbr[ix], outNbr[iy],  isbuf.get(0), isbuf.get(1)));
			machine.sendrecv(ix, sendbuf[ix], seg1+seg2+2, outNbr[ix], recbuf[ix], irbuf.capacity(), outNbr[iy], depth*3 + tag);
			seg1 = irbuf.get();
			seg2 = irbuf.get();
			if (trace > 1) log(String.format("Config layer %d machine %d got msg from %d, sizes %d %d\n", depth, imachine, outNbr[iy], seg1, seg2));

			IVec toout = new IVec(seg1);
			IVec fromout = new IVec(seg2);
			irbuf.get(toout.data, 0, seg1);
			irbuf.get(fromout.data, 0, seg2);	
			top[ix] = toout;
			fromp[ix] = fromout;
		}
	}

	public void config(IVec toi, IVec fromi, IVec [] outputs, int round) {
		config_pre(toi, fromi, "Config ");

		Future<?> [] futures = new Future<?>[k];
		Future<?> timeoutf = executor.submit(new TimeoutThread(machine, machine.configTimeout, futures, "config "+round));
		for (int i = 0; i < k; i++) {
			int ix = (i + posInMyGroup) % k;                               // Try to stagger the traffic
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
	
	public void config_pre(IVec toi, IVec fromi, String prefix) {
		topart = partitioner.part(toi);                                   // Compute partition indices
		frompart = partitioner.part(fromi);
		interleave = new IVec[k];
		toindsparts = partitioner.partition(toi, topart, null);          // Actually partition the indices into k IVec parts
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
	
	public void config_post(Future<?> [] futures, IVec toi, IVec fromi, IVec [] outputs, String prefix) {
		IVec [] totree = new IVec[k];
		IVec [] fromtree = new IVec[k];
		for (int i = 0; i < k; i++) {
			if (futures[i].isDone()) {
				IVec.treeAdd(toindsparts[i], totree);
				IVec.treeAdd(fromindsparts[i], fromtree);
			}
		}
		IVec tomaster = IVec.treeFlush(totree);
		IVec frommaster = IVec.treeFlush(fromtree);
		ton = tomaster.size();
		fromn = fromi.size();
		for (int i = 0; i < k; i++) {
			if (futures[i].isDone()) {
				toMaps[i] = IVec.mapInds(toindsparts[i], tomaster);
				fromMaps[i] = IVec.mapInds(fromindsparts[i], frommaster);
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

	public class ReduceDownThread implements Runnable {
		Vec downv;
		int ix;
		int iy;
		int stride; 
		int tag;
		Vec [] res;

		public ReduceDownThread(Vec downv0, Vec [] res0, int ix0, int iy0, int stride0, int tag0) {
			downv = downv0;
			ix = ix0;
			iy = iy0;
			stride = stride0;
			tag = tag0;
			res = res0;
		}

		public void run() {
			sendbuf[ix].clear();
			recbuf[ix].clear();
			IntBuffer isbuf = sendbuf[ix].asIntBuffer();
			IntBuffer irbuf = recbuf[ix].asIntBuffer();
			FloatBuffer fsbuf = sendbuf[ix].asFloatBuffer();
			FloatBuffer frbuf = recbuf[ix].asFloatBuffer();
			int msize = downv.size();
			isbuf.put(msize);
			fsbuf.position(1);
			fsbuf.put(downv.data, 0, msize);

			if (trace > 1) log(String.format("Reduce down layer %d machine %d sent msg to %d, from %d, size %d\n", depth, imachine, outNbr[ix],  outNbr[iy],  msize));
			machine.sendrecv(ix, sendbuf[ix], msize+1, outNbr[ix], recbuf[ix], frbuf.capacity(), outNbr[iy], depth*3+1 + tag);
			msize = irbuf.get();
			if (trace > 1) log(String.format("Reduce down layer %d machine %d got msg from %d, size %d\n", depth, imachine, outNbr[iy], msize));

			res[ix] = new Vec(msize);
			frbuf.position(1);
			frbuf.get(res[ix].data, 0, msize);
			if (msize != toMaps[ix].size() * stride) {
				throw new RuntimeException(String.format("Exception in ReduceDownThread,  mismatched sizes %d %d", msize, toMaps[ix].size()*stride));
			}
		}
	}

	public Vec reduceDown(Vec downv, int stride, int round) {
		Vec [] vparts = partitioner.partition(downv, topart, topartsizes, stride);
		Vec [] res = new Vec[k];
		Future<?> [] futures = new Future<?>[k];
		Future<?> timeoutf = executor.submit(new TimeoutThread(machine, machine.reduceTimeout, futures, "reduceDown "+round));
		for (int i = 0; i < k; i++) {
			int ix = (i + posInMyGroup) % k;  
			int iy = (posInMyGroup - i + k) % k;  
			futures[i] = executor.submit(new ReduceDownThread(vparts[ix], res, ix, iy, stride, 3*D*round));
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
		Vec newv = new Vec(ton*stride);
		for (int i = 0; i < k; i++) {
			res[i].addTo(newv, toMaps[i], stride);
		}
		return newv;
	}

	public class ReduceUpThread implements Runnable {
		Vec newv;
		Vec upv;
		int ix;
		int iy;
		int stride;
		int tag;

		public ReduceUpThread(Vec newv0, Vec upv0, int ix0, int iy0, int stride0, int tag0) {
			newv = newv0;
			upv = upv0;
			ix = ix0;
			iy = iy0;
			stride = stride0;
			tag = tag0;
		}

		public void run () {
			sendbuf[ix].clear();
			recbuf[ix].clear();
			IntBuffer isbuf = sendbuf[ix].asIntBuffer();
			IntBuffer irbuf = recbuf[ix].asIntBuffer();
			FloatBuffer fsbuf = sendbuf[ix].asFloatBuffer();
			FloatBuffer frbuf = recbuf[ix].asFloatBuffer();
			Vec up = upv.mapFrom(fromMaps[ix], stride);
			int msize = up.size();
			isbuf.put(msize);
			fsbuf.position(1);
			fsbuf.put(up.data, 0, msize);

			if (trace > 1) log(String.format("Reduce up layer %d machine %d sent msg to %d, from %d, size %d\n", depth, imachine, outNbr[ix],  outNbr[iy],  msize));
			machine.sendrecv(ix, sendbuf[ix], msize+1, outNbr[iy], recbuf[ix], irbuf.capacity(), outNbr[ix], depth*3+2 + tag);
			msize = irbuf.get();
			if (trace > 1) log(String.format("Reduce up layer %d machine %d got msg from %d, size %d\n", depth, imachine, outNbr[iy], msize));

			int psize = interleave[ix].size();
			if (msize != psize*stride) throw new RuntimeException("ReduceUp size mismatch "+msize+" "+(psize*stride));
			frbuf.position(1);
			if (fromparts[ix] == null || fromparts[ix].size() != msize) fromparts[ix] = new Vec(msize);
			frbuf.get(fromparts[ix].data, 0, msize);
		}
	}

	public Vec reduceUp(Vec upv, int stride, int round) {
		Vec newv = new Vec(fromn*stride);
		Future<?> [] futures = new Future<?>[k];
		Future<?> timeoutf = executor.submit(new TimeoutThread(machine, machine.reduceTimeout, futures, "reduceUp " + round));
		for (int i = 0; i < k; i++) {
			int ix = (i + posInMyGroup) % k;
			int iy = (posInMyGroup - i + k) % k;
			futures[i] = executor.submit(new ReduceUpThread(newv, upv, ix, iy, stride, 3*D*round));
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
		partitioner.merge(newv, stride, fromparts, interleave);
		return newv;
	}
	
	class ConfigReduceThread implements Runnable {
		IVec [] top;
		IVec [] fromp;
		Vec downv;
		Vec [] res;
		int stride;
		int ix;
		int iy;
		int repno;
		int tag;

		public ConfigReduceThread(IVec [] top0, IVec [] fromp0,	Vec downv0, Vec [] res0, int ix0, int iy0, int stride0, int tag0) {
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
			IntBuffer isbuf = sendbuf[ix].asIntBuffer();
			IntBuffer irbuf = recbuf[ix].asIntBuffer();	
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
			fsbuf.position(3+seg1+seg2);
			fsbuf.put(downv.data, 0, msize);

			if (trace > 1) log(String.format("ConfigReduce layer %d machine %d sent msg to %d, from %d, sizes %d %d %d\n", depth, imachine, outNbr[ix],  outNbr[iy],  isbuf.get(0), isbuf.get(1), isbuf.get(2)));
			machine.sendrecv(ix, sendbuf[ix], seg1+seg2+msize+3, outNbr[ix], recbuf[ix], frbuf.capacity(), outNbr[iy], depth*3 + tag);
			seg1 = irbuf.get();
			seg2 = irbuf.get();
			msize = irbuf.get();
			if (trace > 1) log(String.format("ConfigReduce layer %d machine %d got msg from %d, sizes %d %d %d\n", depth, imachine, outNbr[iy], seg1, seg2, msize));

			IVec toout = new IVec(seg1);
			IVec fromout = new IVec(seg2);
			irbuf.get(toout.data, 0, seg1);
			irbuf.get(fromout.data, 0, seg2);	
			top[ix] = toout;
			fromp[ix] = fromout;
			
			res[ix] = new Vec(msize);
			frbuf.position(3+seg1+seg2);
			frbuf.get(res[ix].data, 0, msize);
		}
	}
	
	public Vec configReduce(IVec toi, IVec fromi, IVec [] outputs, Vec downv, int stride, int round) {
		
		config_pre(toi, fromi, "ConfigReduce ");
		
		Vec [] vparts = partitioner.partition(downv, topart, topartsizes, stride);
		Vec [] res = new Vec[k];

		Future<?> [] futures = new Future<?>[k];
		Future<?> timeoutf = executor.submit(new TimeoutThread(machine, machine.configTimeout, futures, "configReduce " + round));
		for (int i = 0; i < k; i++) {
			int ix = (i + posInMyGroup) % k;
			int iy = (posInMyGroup - i + k) % k;  
			futures[i] = executor.submit(new ConfigReduceThread(toindsparts, fromindsparts, vparts[ix], res, ix, iy, stride, 3*D*round));
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