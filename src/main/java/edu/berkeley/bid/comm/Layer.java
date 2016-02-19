package edu.berkeley.bid.comm;

import java.nio.ByteBuffer;
import java.nio.FloatBuffer;
import java.nio.IntBuffer;
import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.FutureTask;

class Layer {

	/* Layer Configuration Variables */

	int k;        																					 // Size of this group
	int cumk;                                                // cumulative product of sizes
	int cumPos;                                              // cumulative position in groups
	int depth;
	int posInMyGroup;                                        // Position in this machines group
	int [] outNbr;                                           // Machines we talk to 
	int [] inNbr;                                            // Machines we listen to
	IVec [] downMaps;                                        // Maps to indices below for down indices
	IVec [] upMaps;                                          // Maps to indices below for up indices
	IVec [] interleave;                                      // Indices to undo partitioning for up data
	IVec dpart;                                              // partition number for each element
	IVec upart;
	IVec [] downpartinds;
	IVec [] uppartinds;
	Vec [] upparts;
	int trace;
	int imachine;
	int downn;                                               // Size of the down master list
	int upn;
	int timeout;
	Partitioner partitioner;
	Machine machine;
	Groups groups;
	ExecutorService executor;
	AllReduceY allreducer;
	ByteBuffer [] sendbuf;                                     // buffers, one for each destination in a group
	ByteBuffer [] recbuf;


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
		allreducer = machine.allreducer;
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
				synchronized (allreducer) {
					System.out.format("config layer %d machine %d sent msg to %d, from %d, sizes %d %d\n", depth, imachine, outNbr[i],  inNbr[i],  sbuf.get(0), sbuf.get(1));
				}
			}
			machine.sendrecv(i, sendbuf, seg2+2, outNbr[i], recbuf, rbuf.capacity(), inNbr[i], depth*3);
			seg1 = rbuf.get();
			seg2 = rbuf.get();
			if (trace > 1) {
				synchronized (allreducer) {
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
		dpart = partitioner.part(downi);
		upart = partitioner.part(upi);
		downpartinds = partitioner.partition(downi, dpart, null);
		uppartinds = partitioner.partition(upi, upart, interleave);
		IVec [] dtree = new IVec[2*k-1];
		IVec [] utree = new IVec[2*k-1];
		if (trace > 0) {
			synchronized (allreducer) {
				System.out.format("machine %d layer %d, dparts (%d", imachine, depth, downpartinds[0].size());
			}
		}
		
		Future<?> [] futures = new Future<?>[k];
		Future<?> timeoutf = executor.submit(new FutureTask<Void>(new TimeoutThread(timeout, futures), null));
		for (int i = 0; i < k; i++) {
			int ix = (i + posInMyGroup) % k;                               // Try to stagger the traffic
			futures[i] = executor.submit(new FutureTask<Void>(new ConfigThread(downpartinds, uppartinds, dtree, utree, ix), null));
		};
		for (int i = 0; i < k; i++) {
			try {
				futures[i].get();
			} catch (Exception e) {}
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
				synchronized (allreducer) { 
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

		// Note, these have to be partitions of the data now. 

		public ReduceDownThread(Vec newv0, Vec downv0, int i0) {
			newv = newv0;
			downv = downv0;
			i = i0;
		}

		public void run() {
			sendbuf[i].clear();
			recbuf[i].clear();
			IntBuffer isbuf = sendbuf[i].asIntBuffer();
			IntBuffer irbuf = recbuf[i].asIntBuffer();
			FloatBuffer sbuf = sendbuf[i].asFloatBuffer();
			FloatBuffer rbuf = recbuf[i].asFloatBuffer();
			int msize = downv.size();
			isbuf.put(msize);
			sbuf.position(1);
			sbuf.put(downv.data, 0, msize);

			if (trace > 1) {
				synchronized (allreducer) {
					System.out.format("reduce layer %d machine %d sent msg to %d, from %d, size %d\n", depth, imachine, outNbr[i],  inNbr[i],  msize);
				}
			}
			machine.sendrecv(i, sendbuf, msize+1, outNbr[i], recbuf, rbuf.capacity(), inNbr[i], depth*3+1);
			msize = irbuf.get();
			if (trace > 1) {
				synchronized (allreducer) {
					System.out.format("reduce layer %d machine %d got msg from %d, size %d\n", depth, imachine, inNbr[i], msize);
				}
			}

			Vec res = new Vec(msize);
			rbuf.position(1);
			rbuf.get(res.data, 0, msize);
			synchronized (newv) {
				res.addTo(newv, downMaps[i]);	
			}
		}
	}

	public Vec reduceDown(Vec downv, int stride) {
		Vec newv = new Vec(downn);
		Vec [] vparts = partitioner.partition(downv, dpart, downpartinds, stride);
		Future<?> [] futures = new Future<?>[k];
		Future<?> timeoutf = executor.submit(new FutureTask<Void>(new TimeoutThread(timeout, futures), null));
		for (int i = 0; i < k; i++) {
			int ix = (i + posInMyGroup) % k;                               // Try to stagger the traffic
			futures[i] = executor.submit(new FutureTask<Void>(new ReduceDownThread(newv, vparts[ix], ix), null));
		};
		for (int i = 0; i < k; i++) {
			try {
				futures[i].get();
			} catch (Exception e) {}
		}
		timeoutf.cancel(true);
		return newv;
	}

	public class ReduceUpThread implements Runnable {
		Vec newv;
		Vec upv;
		int i;

		public ReduceUpThread(Vec newv0, Vec upv0, int i0) {
			newv = newv0;
			upv = upv0;
			i = i0;
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
				synchronized (allreducer) {
					System.out.format("reduce up layer %d machine %d sent msg to %d, from %d, size %d\n", depth, imachine, outNbr[i],  inNbr[i],  msize);
				}
			}
			machine.sendrecv(i, sendbuf, msize+1, inNbr[i], recbuf, irbuf.capacity(), outNbr[i], depth*3+2);
			msize = irbuf.get();
			if (trace > 1) {
				synchronized (allreducer) {
					System.out.format("reduce up layer %d machine %d got msg from %d, size %d\n", depth, imachine, inNbr[i], msize);
				}
			}

			int psize = upMaps[i].size();
			if (msize != psize) throw new RuntimeException("ReduceUp size mismatch "+msize+" "+psize);
			rbuf.position(1);
			rbuf.get(upparts[i].data, 0, msize);
		}
	}

	public Vec reduceUp(Vec upv, int stride) {
		Vec newv = new Vec(upn);
		Future<?> [] futures = new Future<?>[k];
		executor.submit(new FutureTask<Void>(new TimeoutThread(timeout, futures), null));
		for (int i = 0; i < k; i++) {
			int ix = (i + posInMyGroup) % k;                               // Try to stagger the traffic
			futures[i] = executor.submit(new FutureTask<Void>(new ReduceUpThread(newv, upv, ix), null));
		};
		for (int i = 0; i < k; i++) {
			try {
				futures[i].get();
			} catch (Exception e) {}
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