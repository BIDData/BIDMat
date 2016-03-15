package edu.berkeley.bid.comm;

import java.util.Random;
import java.lang.Math;

public class Groups {
	public final int N;
	public int [][] groupIds;
	public int [][] posInGroups;
	public int [][][] nbrs;
	public int [][] groupSizes;
	public int [] machineCodes;
	public int [] perm;
  public int [] invperm;
	public int [] gmods;
	public int gprod;
	public int D;
	public float initg;
	public int trace;
	
	public Groups(int N0, float initg0) {
		N = N0;
		initg = initg0;
		trace = 0;
		getGroupSizes(N);			
	}
	
	public Groups(int N0, float initg0, int [] machineCodes0, int seed) {
		N = N0;
		initg = initg0;
		machineCodes = machineCodes0;
		trace = 0;
		getGroupSizes(N);	
		createPerm(seed);
		assignGroups();
	}
	
	public Groups(int N0, int [] gmods0, int [] machineCodes0, int seed) {
		N = N0;
		machineCodes = machineCodes0;
		trace = 0;
		gmods = gmods0;
		D = gmods0.length;
		gprod = 1;
		for (int i = 0; i < D; i++) gprod *= gmods[i];
		createPerm(seed);
		assignGroups();
	}

	
	static private void permuteArray(int [] arr, int n, int seed) {
		Random gen = new Random(seed);
		for (int i = 0; i < n; i++) {
			arr[i] = i;
		}
		for (int i = 0; i < (n-1); i++) {
			int iswap = Math.min(n - 1, i + (int)(gen.nextDouble() * (n-i)));
			int tmp = arr[i];
			arr[i] = arr[iswap];
			arr[iswap] = tmp;
		}
	}
	
	private void createPerm(int seed) {
		perm = new int[N];
		invperm = new int[N];
		permute(seed);
	}
			
	public void permute(int seed) {
		permuteArray(perm, N, seed);
		for (int i = 0; i < N; i++) {
			invperm[perm[i]] = i;
		}
	}
	// Get the target layer sizes
	// Layer sizes form a geometric series, starting with exp(vD);
	//
	// Solve for Sum vi = log N
	// where vi = v0 + i alpha
	// alpha is the log of the layer growth factor
	public void getGroupSizes(int N) {
		double [] gsizes;
		double alpha = 0.3;
		double vD = Math.log(initg);
		double a = alpha / 2;
		double b = vD - alpha / 2;
		double c = -Math.log(N);
		double rad = Math.sqrt(b*b - 4*a*c);
		double soln = (-b + rad)/(2*a);
		D = Math.max(1, (int)Math.round(soln));
		gsizes = new double[D];
		gmods = new int[D];
		vD = -c / D - (D-1) * alpha / 2;
		gprod = 1;
		for (int i = 0; i < D-1; i++) {
			gsizes[i] = Math.exp(vD + (D-i-1) * alpha);
			gmods[i] = (int)Math.round(gsizes[i]);
			gprod *= gmods[i];
		}
		gmods[D-1] = 1 + (N-1) / gprod;                    // Make sure gprod >= N 
		gprod = gprod * gmods[D-1];  
	}
	
	public void assignGroups() {
		groupIds = new int[N][];
		posInGroups = new int[D][];
		nbrs = new int[D][][];
		groupSizes = new int[D][];
		for (int i = 0; i < N; i++) {
			groupIds[i] = new int[D];
		}
		int pprod = 1;

		for (int d = 0; d < D; d++) {
			int numgroups = gprod / gmods[d];
			posInGroups[d] = new int[N];
			nbrs[d] = new int[numgroups][];
			groupSizes[d] = new int[numgroups];
			int []groupPos = new int[numgroups];

			for (int i = 0; i < N; i++) {
				int ii = machineCodes[i];
				int left = ii / (pprod * gmods[d]);
				int right = ii % pprod;
				int gnum = right + pprod * left;
				groupIds[i][d] = gnum;
				groupSizes[d][gnum]++;
			}

			for (int i = 0; i < numgroups; i++) {
				nbrs[d][i] = new int[groupSizes[d][i]];
			}
			
			for (int i = 0; i < N; i++) {
				int gnum = groupIds[i][d];
				int pos = groupPos[gnum];
				nbrs[d][gnum][pos] = i;
				posInGroups[d][i] = pos;
				groupPos[gnum]++;
			}
			
			pprod *= gmods[d];
		}
	}
	
	public int [] minSizes() {
		int [] msizes = new int[D];
		for (int d = 0; d < D; d++) {
			int mm = 1000000;
			for (int j = 0; j < groupSizes[d].length; j++) {
				if (groupSizes[d][j] < mm) {
					mm = groupSizes[d][j];
				}
			}
			msizes[d] = mm;
		}
		return msizes;
	}
	
	public int [] maxSizes() {
		int [] msizes = new int[D];
		for (int d = 0; d < D; d++) {
			int mm = 0;
			for (int j = 0; j < groupSizes[d].length; j++) {
				if (groupSizes[d][j] > mm) {
					mm = groupSizes[d][j];
				}
			}
			msizes[d] = mm;
		}
		return msizes;
	}
	
	public int compare(Groups g) {
		int [] mymin = minSizes();
		int [] thatmin = g.minSizes();
		int minc = 1;
		int maxc = -1;
		int d = 0;
		while (maxc - minc < 2 && d < D) {
			int cmp = (int)Math.signum(mymin[d] - thatmin[d]);
			if (cmp > maxc) maxc = cmp;
			if (cmp < minc) minc = cmp;
			d++;
		}			
		return maxc + minc;
	}
	
	public void printArray(int [] arr) {
		for (int i = 0; i < arr.length-1; i++) {
			System.out.format("%d,", arr[i]);
		}
		System.out.format("%d", arr[arr.length-1]);
	}
	
	static final double mapfn(double v) {
		return Math.sqrt(v);
	}
	//
	// Minimize the sum of the squares of the group sizes
	//
	
	public int [] optimize(int howlong, double prob) {
		getGroupSizes(N);
		machineCodes = new int[gprod];
		groupSizes = new int[D][];
		groupIds = new int[gprod][];
		permuteArray(machineCodes, gprod, 1000);
		double sumv = 0;
		for (int d = 0; d < D; d++) {
			groupSizes[d] = new int[gprod/gmods[d]];
		}
		for (int i = 0; i < gprod; i++) {
			groupIds[i] = new int[D];
			int ii = machineCodes[i];
			int pprod = 1;
			for (int d = 0; d < D; d++) {
				int left = ii / (pprod * gmods[d]);
				int right = ii % pprod;
				int gnum = right + pprod * left;
				groupIds[i][d] = gnum;
				if (i < N) {
					groupSizes[d][gnum]++;
				}
				pprod *= gmods[d];
			}
		}
		for (int d = 0; d < D; d++) {
			for (int i = 0; i < groupSizes[d].length; i++) {
				sumv += mapfn(groupSizes[d][i]);
			}
		}
		Random gen = new Random(1231231);
		for (int iter = 0; iter < howlong; iter++) {
			int i = (int)Math.min(N-1, (gen.nextDouble() * N));
			int j = (int)Math.min(gprod-1, N + (gen.nextDouble() * (gprod - N)));
			// Remove i's counts and add j's
			double testsumv = sumv;
			for (int d = 0; d < D; d++) {
				int gi = groupIds[i][d];
				int gj = groupIds[j][d];
				int sizei = groupSizes[d][gi];
				int sizej = groupSizes[d][gj];
				if (gi != gj) {
					testsumv += mapfn(sizej+1) - mapfn(sizej) + mapfn(sizei-1) - mapfn(sizei);
				}
			}
			if (testsumv > sumv || gen.nextGaussian()*prob > sumv - testsumv) { 
				sumv = testsumv;
				
				for (int d = 0; d < D; d++) {
					int gi = groupIds[i][d];
					int gj = groupIds[j][d];
					groupSizes[d][gi]--;
					groupSizes[d][gj]++;
				}
				
				int mm = machineCodes[i];
				machineCodes[i] = machineCodes[j];
				machineCodes[j] = mm;
				
				int [] vv = groupIds[i];
				groupIds[i] = groupIds[j];
				groupIds[j] = vv;	
				
				if (trace > 0) {		
					System.out.format("sumv = %4.3f, ",testsumv);
					System.out.print("min=[");
					printArray(minSizes());
					System.out.print("], ");
					System.out.print("max=[");
					printArray(maxSizes());
					System.out.print("]\n");
				}
			}
		}		
		
		return machineCodes;
	}

	public int groupId(int imachine, int level) {
		return groupIds[perm[imachine]][level];
	}

	public int posInGroup(int imachine, int level) {
		return posInGroups[level][perm[imachine]];
	}

	public int [] nodesInGroup(int imachine, int level) {
		int g = groupIds[perm[imachine]][level];
		int [] grp = nbrs[level][g];
		int [] out = new int[grp.length];
		for (int i = 0; i < grp.length; i++) {
			out[i] = invperm[grp[i]];
		}
		return out;
	}

	public int depth() {
		return D;
	}
}