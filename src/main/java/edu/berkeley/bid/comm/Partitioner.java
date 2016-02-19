package edu.berkeley.bid.comm;

public class Partitioner {
	int cumk;
	int ibase;
	int k;

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
	
	public IVec part(LVec v) {
		IVec parts = new IVec(v.size());
		for (int i = 0; i < v.size(); i++) {
			int r = (int)((v.data[i] - ibase) / cumk);
			parts.data[i] = r % k;
		}
		return parts;
	}
	
	public IVec [] partition(IVec vv, IVec part, IVec [] mapback) {
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
			mapback[i] = new IVec(lens[i]);
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
	
	public Vec [] partition(Vec vv, IVec part, IVec [] iparts, int stride) {
		int n = vv.size();
		int [] nlens = new int[k];
		Vec [] parts = new Vec[k];
		for (int i = 0; i < k; i++) {
			parts[i] = new Vec(iparts[i].size() * stride);
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
		for (int i = 0; i < k; i++) {
			Vec part = parts[i];
			IVec mapb = mapback[i];
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