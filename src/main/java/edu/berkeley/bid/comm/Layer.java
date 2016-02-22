package edu.berkeley.bid.comm;

import java.nio.ByteBuffer;
import java.nio.FloatBuffer;
import java.nio.IntBuffer;
import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

public class Layer {

	/* Layer Configuration Variables */

	public final int k;        																					 // Size of this group
	public final int cumk;                                               // cumulative product of sizes
	public final int cumPos;                                             // cumulative position in groups
	public final int depth;
	public final int posInMyGroup;                                       // Position in this machines group
	public int [] outNbr;                                                // Machines we talk to 
	public int [] inNbr;                                                 // Machines we listen to
	public IVec [] downpartinds;                                         // partition indices for down indices
	public IVec [] uppartinds;
	public IVec downpartsizes;
	public IVec uppartsizes;
	public IVec [] downMaps;                                             // Maps to indices below for down indices
	public IVec [] upMaps;                                               // Maps to indices below for up indices
	public IVec [] interleave;                                           // Indices to undo partitioning for up data
	public IVec dpart;                                                   // partition number for each element
	public IVec upart;
	public Vec [] upparts;
	public int trace;
	public final int imachine;
	public int downn;                                                    // Size of the down master list
	public int upn;
	public int timeout;
	public Partitioner partitioner;
	public Machine machine;
	public Groups groups;
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
		timeout = machine.timeout;
		executor = machine.executor;
		partitioner = new Partitioner(cumk, cumPos, k);
		inNbr = groups.nodesInGroup(imachine, depth);
		outNbr = groups.nodesInGroup(imachine, depth);		
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

		public ConfigThread(IVec [] downp0, IVec [] upp0,	IVec [] dtree0, IVec [] utree0, int i0) {
			downp = downp0;
			upp = upp0;					
			dtree = dtree0;
			utree = utree0;
			i = i0;
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
				synchronized (network) {
					System.out.format("config layer %d machine %d sent msg to %d, from %d, sizes %d %d\n", depth, imachine, outNbr[i],  inNbr[i],  sbuf.get(0), sbuf.get(1));
				}
			}
			machine.sendrecv(i, sendbuf, seg2+2, outNbr[i], recbuf, rbuf.capacity(), inNbr[i], depth*3);
			seg1 = rbuf.get();
			seg2 = rbuf.get();
			if (trace > 1) {
				synchronized (network) {
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
		}
	}

	public void config(IVec downi, IVec upi, IVec [] outputs) {
		dpart = partitioner.part(downi);                                   // Compute partition indices
		upart = partitioner.part(upi);
		interleave = new IVec[k];
		downpartinds = partitioner.partition(downi, dpart, null);          // Actually partition the indices into k IVec parts
		uppartinds = partitioner.partition(upi, upart, interleave);
		downpartsizes = new IVec(k);
		uppartsizes = new IVec(k);
		upparts = new Vec[k];
		for (int i = 0; i < k; i++) {
			downpartsizes.data[i] = downpartinds[i].length;
			uppartsizes.data[i] = uppartinds[i].length;
		}
		IVec [] dtree = new IVec[2*k-1];
		IVec [] utree = new IVec[2*k-1];
		if (trace > 0) {
			synchronized (network) {
				System.out.format("machine %d layer %d, dparts (%d", imachine, depth, downpartinds[0].size());
			}
		}

		Future<?> [] futures = new Future<?>[k];
		Future<?> timeoutf = executor.submit(new TimeoutThread(timeout, futures), null);
		for (int i = 0; i < k; i++) {
			int ix = (i + posInMyGroup) % k;                               // Try to stagger the traffic
			futures[i] = executor.submit(new ConfigThread(downpartinds, uppartinds, dtree, utree, ix), null);
		};
		for (int i = 0; i < k; i++) {
			try {
				futures[i].get();
			} catch (Exception e) {
				synchronized (network) {
					e.printStackTrace();
				}
			}
		}	
		timeoutf.cancel(true);

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
				synchronized (network) { 
					System.out.format("config machine %d dmap(%d) size %d\n", imachine, i, downMaps[i].size());
					System.out.format("config machine %d umap(%d) size %d\n", imachine, i, upMaps[i].size());
				}
			}
		}
		outputs[0] = dmaster;
		outputs[1] = umaster;
	}

	public class ReduceDownThread implements Runnable {
		Vec newv;
		Vec downv;
		int ix;
		int stride;

		// Note, these have to be partitions of the data now. 

		public ReduceDownThread(Vec newv0, Vec downv0, int i0, int stride0) {
			newv = newv0;
			downv = downv0;
			ix = i0;
			stride = stride0;
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

			if (trace > 1) {
				synchronized (network) {
					System.out.format("reduce down layer %d machine %d sent msg to %d, from %d, size %d\n", depth, imachine, outNbr[ix],  inNbr[ix],  msize);
				}
			}
			machine.sendrecv(ix, sendbuf, msize+1, outNbr[ix], recbuf, rbuf.capacity(), inNbr[ix], depth*3+1);
			msize = irbuf.get();
			if (trace > 1) {
				synchronized (network) {
					System.out.format("reduce down layer %d machine %d got msg from %d, size %d\n", depth, imachine, inNbr[ix], msize);
				}
			}

			Vec res = new Vec(msize);
			rbuf.position(1);
			rbuf.get(res.data, 0, msize);
			synchronized (newv) {
				res.addTo(newv, downMaps[ix], stride);	
			}
		}
	}

	public Vec reduceDown(Vec downv, int stride) {
		Vec newv = new Vec(downn*stride);
		Vec [] vparts = partitioner.partition(downv, dpart, downpartsizes, stride);
		Future<?> [] futures = new Future<?>[k];
		Future<?> timeoutf = executor.submit(new TimeoutThread(timeout, futures));
		for (int i = 0; i < k; i++) {
			int ix = (i + posInMyGroup) % k;                               // Try to stagger the traffic
			futures[i] = executor.submit(new ReduceDownThread(newv, vparts[ix], ix, stride));
		};
		for (int i = 0; i < k; i++) {
			try {
				futures[i].get();
			} catch (Exception e) {
				synchronized (network) {
					e.printStackTrace();
				}
			}
		}
		timeoutf.cancel(true);
		return newv;
	}

	public class ReduceUpThread implements Runnable {
		Vec newv;
		Vec upv;
		int i;
		int stride;

		public ReduceUpThread(Vec newv0, Vec upv0, int i0, int stride0) {
			newv = newv0;
			upv = upv0;
			i = i0;
			stride = stride0;
		}

		public void run () {
			sendbuf[i].clear();
			recbuf[i].clear();
			IntBuffer isbuf = sendbuf[i].asIntBuffer();
			IntBuffer irbuf = recbuf[i].asIntBuffer();
			FloatBuffer sbuf = sendbuf[i].asFloatBuffer();
			FloatBuffer rbuf = recbuf[i].asFloatBuffer();
			Vec up = upv.mapFrom(upMaps[i], stride);
			int msize = up.size();
			isbuf.put(msize);
			sbuf.position(1);
			sbuf.put(up.data, 0, msize);

			if (trace > 1) {
				synchronized (network) {
					System.out.format("reduce up layer %d machine %d sent msg to %d, from %d, size %d\n", depth, imachine, outNbr[i],  inNbr[i],  msize);
				}
			}
			machine.sendrecv(i, sendbuf, msize+1, inNbr[i], recbuf, irbuf.capacity(), outNbr[i], depth*3+2);
			msize = irbuf.get();
			if (trace > 1) {
				synchronized (network) {
					System.out.format("reduce up layer %d machine %d got msg from %d, size %d\n", depth, imachine, inNbr[i], msize);
				}
			}

			int psize = interleave[i].size();
			if (msize != psize*stride) throw new RuntimeException("ReduceUp size mismatch "+msize+" "+(psize*stride));
			rbuf.position(1);
			if (upparts[i] == null || upparts[i].size() != msize) upparts[i] = new Vec(msize);
			rbuf.get(upparts[i].data, 0, msize);
		}
	}

	public Vec reduceUp(Vec upv, int stride) {
		Vec newv = new Vec(upn*stride);
		Future<?> [] futures = new Future<?>[k];
		executor.submit(new TimeoutThread(timeout, futures), null);
		for (int i = 0; i < k; i++) {
			int ix = (i + posInMyGroup) % k;                               // Try to stagger the traffic
			futures[i] = executor.submit(new ReduceUpThread(newv, upv, ix, stride), null);
		};
		for (int i = 0; i < k; i++) {
			try {
				futures[i].get();
			} catch (Exception e) {
				synchronized (network) {
					e.printStackTrace();
				}
			}
		}
		partitioner.merge(newv, stride, upparts, interleave);
		return newv;
	}
	
	public class TimeoutThread implements Runnable {
		Future <?> [] futures;
		int mtime;

		public TimeoutThread(int mtime0, Future <?> [] futures0) {		
			futures = futures0;
			mtime = mtime0;
		}

		public void run() {
			try {
				Thread.sleep(mtime);

				for (int i = 0; i < futures.length; i++) {
					if (futures[i] != null) {
						futures[i].cancel(true);
					}
				}
			} catch (InterruptedException e) {
			}
		}
	}
}