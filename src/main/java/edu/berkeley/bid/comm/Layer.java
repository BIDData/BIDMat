package edu.berkeley.bid.comm;

import java.nio.ByteBuffer;
import java.nio.FloatBuffer;
import java.nio.IntBuffer;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

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
		int i;
		int repno;
		int tag;

		public ConfigThread(IVec [] top0, IVec [] fromp0,	int i0, int tag0) {
			top = top0;
			fromp = fromp0;					
			i = i0;
			tag = tag0;
		}

		public void run() {
			sendbuf[i].clear();
			recbuf[i].clear();
			IntBuffer sbuf = sendbuf[i].asIntBuffer();
			IntBuffer rbuf = recbuf[i].asIntBuffer();	
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
			machine.sendrecv(i, sendbuf, seg1+seg2+2, outNbr[i], recbuf, rbuf.capacity(), inNbr[i], depth*3 + tag);
			seg1 = rbuf.get();
			seg2 = rbuf.get();
			if (trace > 1 && network != null) {
				synchronized (network) {
					System.out.format("Config layer %d machine %d got msg from %d, sizes %d %d\n", depth, imachine, inNbr[i], seg1, seg2);
				}
			}

			IVec toout = new IVec(seg1);
			IVec fromout = new IVec(seg2);
			rbuf.get(toout.data, 0, seg1);
			rbuf.get(fromout.data, 0, seg2);	
			top[i] = toout;
			fromp[i] = fromout;
		}
	}

	public void config(IVec toi, IVec fromi, IVec [] outputs, int round) {
		config_pre(toi, fromi, "Config ");

		Future<?> [] futures = new Future<?>[k];
		Future<?> timeoutf = executor.submit(new TimeoutThread(machine.configTimeout, futures, "config "+round));
		for (int i = 0; i < k; i++) {
			int ix = (i + posInMyGroup) % k;                               // Try to stagger the traffic
			futures[i] = executor.submit(new ConfigThread(toindsparts, fromindsparts, ix, 3*D*round), null);
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
		if (trace > 1 && network != null) {
			synchronized (network) {
				System.out.format(prefix + "machine %d layer %d, toparts %d\n", imachine, depth, toindsparts[0].size());
			}
		}
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
				if (trace > 1 && network != null) {
					synchronized (network) { 
						System.out.format(prefix + "machine %d layer %d dmap(%d) size %d\n", imachine, depth, i, toMaps[i].size());
						System.out.format(prefix + "machine %d layer %d umap(%d) size %d\n", imachine, depth, i, fromMaps[i].size());
					}
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
		int stride; 
		int tag;
		Vec [] res;

		public ReduceDownThread(Vec downv0, Vec [] res0, int i0, int stride0, int tag0) {
			downv = downv0;
			ix = i0;
			stride = stride0;
			tag = tag0;
			res = res0;
		}

		public void run() {
			sendbuf[ix].clear();
			recbuf[ix].clear();
			IntBuffer isbuf = sendbuf[ix].asIntBuffer();
			IntBuffer irbuf = recbuf[ix].asIntBuffer();
			FloatBuffer sbuf = sendbuf[ix].asFloatBuffer();
			FloatBuffer rbuf = recbuf[ix].asFloatBuffer();
			int msize = downv.size();
			isbuf.put(msize);
			sbuf.position(1);
			sbuf.put(downv.data, 0, msize);

			if (trace > 1 && network != null) {
				synchronized (network) {
					System.out.format("Reduce down layer %d machine %d sent msg to %d, from %d, size %d\n", depth, imachine, outNbr[ix],  inNbr[ix],  msize);
				}
			}
			machine.sendrecv(ix, sendbuf, msize+1, outNbr[ix], recbuf, rbuf.capacity(), inNbr[ix], depth*3+1 + tag);
			msize = irbuf.get();
			if (trace > 1 && network != null) {
				synchronized (network) {
					System.out.format("Reduce down layer %d machine %d got msg from %d, size %d\n", depth, imachine, inNbr[ix], msize);
				}
			}

			res[ix] = new Vec(msize);
			rbuf.position(1);
			rbuf.get(res[ix].data, 0, msize);
			if (msize != toMaps[ix].size() * stride) {
				throw new RuntimeException(String.format("Exception in ReduceDownThread,  mismatched sizes %d %d", msize, toMaps[ix].size()*stride));
			}
		}
	}

	public Vec reduceDown(Vec downv, int stride, int round) {
		Vec [] vparts = partitioner.partition(downv, topart, topartsizes, stride);
		Vec [] res = new Vec[k];
		Future<?> [] futures = new Future<?>[k];
		Future<?> timeoutf = executor.submit(new TimeoutThread(machine.reduceTimeout, futures, "reduceDown "+round));
		for (int i = 0; i < k; i++) {
			int ix = (i + posInMyGroup) % k;                               // Try to stagger the traffic
			futures[i] = executor.submit(new ReduceDownThread(vparts[ix], res, ix, stride, 3*D*round));
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
		Vec newv = new Vec(ton*stride);
		for (int i = 0; i < k; i++) {
			res[i].addTo(newv, toMaps[i], stride);
		}
		return newv;
	}

	public class ReduceUpThread implements Runnable {
		Vec newv;
		Vec upv;
		int i;
		int stride;
		int tag;

		public ReduceUpThread(Vec newv0, Vec upv0, int i0, int stride0, int tag0) {
			newv = newv0;
			upv = upv0;
			i = i0;
			stride = stride0;
			tag = tag0;
		}

		public void run () {
			sendbuf[i].clear();
			recbuf[i].clear();
			IntBuffer isbuf = sendbuf[i].asIntBuffer();
			IntBuffer irbuf = recbuf[i].asIntBuffer();
			FloatBuffer sbuf = sendbuf[i].asFloatBuffer();
			FloatBuffer rbuf = recbuf[i].asFloatBuffer();
			Vec up = upv.mapFrom(fromMaps[i], stride);
			int msize = up.size();
			isbuf.put(msize);
			sbuf.position(1);
			sbuf.put(up.data, 0, msize);

			if (trace > 1 && network != null) {
				synchronized (network) {
					System.out.format("Reduce up layer %d machine %d sent msg to %d, from %d, size %d\n", depth, imachine, outNbr[i],  inNbr[i],  msize);
				}
			}
			machine.sendrecv(i, sendbuf, msize+1, inNbr[i], recbuf, irbuf.capacity(), outNbr[i], depth*3+2 + tag);
			msize = irbuf.get();
			if (trace > 1 && network != null) {
				synchronized (network) {
					System.out.format("Reduce up layer %d machine %d got msg from %d, size %d\n", depth, imachine, inNbr[i], msize);
				}
			}

			int psize = interleave[i].size();
			if (msize != psize*stride) throw new RuntimeException("ReduceUp size mismatch "+msize+" "+(psize*stride));
			rbuf.position(1);
			if (fromparts[i] == null || fromparts[i].size() != msize) fromparts[i] = new Vec(msize);
			rbuf.get(fromparts[i].data, 0, msize);
		}
	}

	public Vec reduceUp(Vec upv, int stride, int round) {
		Vec newv = new Vec(fromn*stride);
		Future<?> [] futures = new Future<?>[k];
		Future<?> timeoutf = executor.submit(new TimeoutThread(machine.reduceTimeout, futures, "reduceUp " + round));
		for (int i = 0; i < k; i++) {
			int ix = (i + posInMyGroup) % k;                               // Try to stagger the traffic
			futures[i] = executor.submit(new ReduceUpThread(newv, upv, ix, stride, 3*D*round));
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
		partitioner.merge(newv, stride, fromparts, interleave);
		return newv;
	}
	
	class ConfigReduceThread implements Runnable {
		IVec [] top;
		IVec [] fromp;
		Vec downv;
		Vec [] res;
		int stride;
		int i;
		int repno;
		int tag;

		public ConfigReduceThread(IVec [] top0, IVec [] fromp0,	Vec downv0, Vec [] res0, int i0, int stride0, int tag0) {
			top = top0;
			fromp = fromp0;					
			downv = downv0;
			res = res0;
			i = i0;
			stride = stride0;
			tag = tag0;
		}

		public void run() {
			sendbuf[i].clear();
			recbuf[i].clear();
			IntBuffer isbuf = sendbuf[i].asIntBuffer();
			IntBuffer irbuf = recbuf[i].asIntBuffer();	
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
			sbuf.position(3+seg1+seg2);
			sbuf.put(downv.data, 0, msize);

			if (trace > 1 && network != null) {
				synchronized (network) {
					System.out.format("ConfigReduce layer %d machine %d sent msg to %d, from %d, sizes %d %d %d\n", depth, imachine, outNbr[i],  inNbr[i],  isbuf.get(0), isbuf.get(1), isbuf.get(2));
				}
			}
			machine.sendrecv(i, sendbuf, seg1+seg2+msize+3, outNbr[i], recbuf, rbuf.capacity(), inNbr[i], depth*3 + tag);
			seg1 = irbuf.get();
			seg2 = irbuf.get();
			msize = irbuf.get();
			if (trace > 1 && network != null) {
				synchronized (network) {
					System.out.format("ConfigReduce layer %d machine %d got msg from %d, sizes %d %d %d\n", depth, imachine, inNbr[i], seg1, seg2, msize);
				}
			}

			IVec toout = new IVec(seg1);
			IVec fromout = new IVec(seg2);
			irbuf.get(toout.data, 0, seg1);
			irbuf.get(fromout.data, 0, seg2);	
			top[i] = toout;
			fromp[i] = fromout;
			
			res[i] = new Vec(msize);
			rbuf.position(3+seg1+seg2);
			rbuf.get(res[i].data, 0, msize);
		}
	}
	
	public Vec configReduce(IVec toi, IVec fromi, IVec [] outputs, Vec downv, int stride, int round) {
		
		config_pre(toi, fromi, "ConfigReduce ");
		
		Vec [] vparts = partitioner.partition(downv, topart, topartsizes, stride);
		Vec [] res = new Vec[k];

		Future<?> [] futures = new Future<?>[k];
		Future<?> timeoutf = executor.submit(new TimeoutThread(machine.configTimeout, futures, "configReduce " + round));
		for (int i = 0; i < k; i++) {
			int ix = (i + posInMyGroup) % k;                               // Try to stagger the traffic
			futures[i] = executor.submit(new ConfigReduceThread(toindsparts, fromindsparts, vparts[ix], res, ix, stride, 3*D*round));
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
		
		config_post(futures, toi, fromi, outputs, "ConfigReduce ");
		
		Vec newv = new Vec(ton*stride);
		for (int i = 0; i < k; i++) {
			res[i].addTo(newv, toMaps[i], stride);
		}
		return newv;
	}
	
	public class TimeoutThread implements Runnable {
		Future <?> [] futures;
		int mtime;
		String canceller;

		public TimeoutThread(int mtime0, Future <?> [] futures0, String canceller0) {		
			futures = futures0;
			mtime = mtime0;
			canceller = canceller0;
		}

		public void run() {
			try {
				Thread.sleep(mtime);

				for (int i = 0; i < futures.length; i++) {
					if (futures[i] != null) {
						if (trace > 1 && network != null) {
							synchronized (network) {
								System.out.format("Machine %d cancelling %s future %d\n", imachine, canceller, i);
							}
						}
						futures[i].cancel(true);
					}
				}
			} catch (InterruptedException e) {
			}
		}
	}
}