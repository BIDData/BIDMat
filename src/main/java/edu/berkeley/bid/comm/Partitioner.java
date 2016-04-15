package edu.berkeley.bid.comm;

public class Partitioner {
	public final int cumk;
	public final int ibase;
	public final int k;

	public Partitioner(int cumk0, int ibase0, int k0) {
		cumk = cumk0;
		k = k0;
		ibase = ibase0;
	}

	public int part(int v) {
		int r = (v - ibase) / cumk;
		return r % k;
	}
	
	public IVec part(IVec v) {
		IVec parts = new IVec(v.size());
		for (int i = 0; i < v.size(); i++) {
			int r = ((v.data[i] - ibase) / cumk);
			parts.data[i] = r % k;
		}
		return parts;
	}
	
	// Take the input indices and map them to values in 0..(k-1) using the partition map
	public IVec part(LVec v) {
		IVec parts = new IVec(v.size());
		for (int i = 0; i < v.size(); i++) {
			int r = (int)((v.data[i] - ibase) / cumk);
			parts.data[i] = r % k;
		}
		return parts;
	}
	
	// Take the input indices, a partition map defined by part() above, and return k IVec's containing the partitioned indices. 
	// Optionally return k IVecs in "mapback" that specify the indices in the original vector for partition indices.
	public IVec [] partition(IVec vv, IVec part, IVec [] mapback) {
		if (part.size() != vv.size()) {
			throw new RuntimeException(String.format("matrix partition: mismatched lengths %d %d", part.size(), vv.size()));
		}
		int n = vv.size();
		int [] lens = new int[k];
		IVec [] parts = new IVec[k];
		for (int i = 0; i < k; i++) {
			lens[i] = 0;
		}
		for (int i = 0; i < n; i++) {
			int r = part.data[i];
			lens[r] ++;
		}
		for (int i = 0; i < k; i++) {
			parts[i] = new IVec(lens[i]);
			if (mapback != null) mapback[i] = new IVec(lens[i]);
			lens[i] = 0;
		}
		if (mapback != null) {
			for (int i = 0; i < n; i++) {
				int v = vv.data[i];
				int r = part.data[i];
				parts[r].data[lens[r]] = v;
				mapback[r].data[lens[r]] = i;
				lens[r] ++;
			}
		} else {
			for (int i = 0; i < n; i++) {
				int v = vv.data[i];
				int r = part.data[i];
				parts[r].data[lens[r]] = v;
				lens[r] ++;
			}
		}
		return parts;
	}
	
	public LVec [] partition(LVec vv, IVec part, IVec [] mapback) {
		if (part.size() != vv.size()) {
			throw new RuntimeException(String.format("matrix partition: mismatched lengths %d %d", part.size(), vv.size()));
		}
		int [] lens = new int[k];
		LVec [] parts = new LVec[k];
		int n = vv.size();
		for (int i = 0; i < k; i++) {
			lens[i] = 0;
		}
		for (int i = 0; i < n; i++) {
			int r = part.data[i];
			lens[r] ++;
		}
		for (int i = 0; i < k; i++) {
			parts[i] = new LVec(lens[i]);
			if (mapback != null) mapback[i] = new IVec(lens[i]);
			lens[i] = 0;
		}
		if (mapback != null) {
			for (int i = 0; i < n; i++) {
				long v = vv.data[i];
				int r = part.data[i];
				parts[r].data[lens[r]] = v;
				mapback[r].data[lens[r]] = i;
				lens[r] ++;
			}
		} else {
			for (int i = 0; i < n; i++) {
				long v = vv.data[i];
				int r = part.data[i];
				parts[r].data[lens[r]] = v;
				lens[r] ++;
			}
		}
		return parts;
	}
	
	public Vec [] partition(Vec vv, IVec part, IVec partsizes, int stride) {
		if (part.size()*stride != vv.size()) {
			throw new RuntimeException(String.format("matrix partition: mismatched lengths %d %d",part.size()*stride, vv.size()));
		}
		int n = part.size();
		int [] nlens = new int[k];
		Vec [] parts = new Vec[k];
		for (int i = 0; i < k; i++) {
			parts[i] = new Vec(partsizes.data[i] * stride);
			nlens[i] = 0;
		}
		for (int i = 0; i < n; i++) {
			int r = part.data[i];
			int ir = nlens[r] * stride;
			int ii = i * stride;
			for (int j = 0; j < stride; j++) {
				parts[r].data[ir + j] = vv.data[j + ii];
			}
			nlens[r] ++;
		}
		return parts;
	}
	
	public void merge(Vec vv, int stride, Vec [] parts, IVec [] mapback) {
		for (int i = 0; i < vv.length; i++) {
			vv.data[i] = 0;
		}
		for (int i = 0; i < k; i++) {
			Vec part = parts[i];
			IVec mapb = mapback[i];
			if (part != null) {
				for (int j = 0; j < mapb.size(); j++) {
					int jj = j * stride;
					int rj = mapb.data[j] * stride;
					for (int kk = 0; kk < stride; kk++) {
						vv.data[kk + rj] = part.data[kk + jj];
					}
				}
			}
		}
	}
	
}