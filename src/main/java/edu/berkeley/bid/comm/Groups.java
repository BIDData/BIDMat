package edu.berkeley.bid.comm;

import java.util.Random;
import java.lang.Math;

public class Groups {
	public final int N;
	public int [][] groupNums;
	public int [][] posInGroups;
	public int [][][] nbrs;
	public int [] perm;
  public int [] invperm;
	public double [] gsizes;
	public int [] gmods;
	public int D;

	private void permute(int seed) {
		perm = new int[N];
		invperm = new int[N];
		Random gen = new Random(seed);
		for (int i = 0; i < N; i++) {
			perm[i] = i;
		}
		for (int i = 0; i < (N-1); i++) {
			int iswap = Math.min(N - 1, i + (int)(gen.nextDouble() * (N-i)));
			int tmp = perm[i];
			perm[i] = perm[iswap];
			perm[iswap] = tmp;
		}
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
		double alpha = 0.3;
		double vD = Math.log(6);
		double a = alpha / 2;
		double b = vD - alpha / 2;
		double c = -Math.log(N);
		double rad = Math.sqrt(b*b - 4*a*c);
		double soln = (-b + rad)/(2*a);
		D = Math.max(1, (int)Math.round(soln));
		gsizes = new double[D];
		gmods = new int[D];
		vD = -c / D - (D-1) * alpha / 2;
		int prod = 1;
		for (int i = 0; i < D; i++) {
			gsizes[i] = Math.exp(vD + (D-i-1) * alpha);
			gmods[i] = (int)Math.round(N/gsizes[i]);
			prod *= gmods[i];
		}
		if (prod < N) gmods[D-1] = 1 + (N-1) / (prod/gmods[D-1]);
	}
	
	public void assignGroups() {
		groupNums = new int[D][];
		posInGroups = new int[D][];
		nbrs = new int[D][][];
		int [][] groupPos = new int[D][];
		for (int d = 0; d < D; d++) {
			groupNums[d] = new int[N];
			posInGroups[d] = new int[N];
			groupPos[d] = new int[gmods[d]];
			nbrs[d] = new int[gmods[d]][];
		}
		for (int i = 0; i < N; i++) {
			int res = i;
			for (int d = 0; d < D; d++) {
				int q = res / gmods[d];
				int v = res - q * gmods[d];
				groupNums[d][i] = v;
				posInGroups[d][i] = groupPos[d][v];
				groupPos[d][v] ++;
				res = q;
			}
		}
		for (int d = 0; d < D; d++) {
			for (int i = 0; i < gmods[d]; i++) {
				nbrs[d][i] = new int[groupPos[d][i]];
				groupPos[d][i] = 0;
			}
			for (int i = 0; i < N; i++) {
				int g = groupNums[d][i];
				int pos = groupPos[d][g];
				nbrs[d][g][pos] = i;
				groupPos[d][g]++;
			}
		}
	}

	public Groups(int N0, int seed) {
		N = N0;
		permute(seed);
		getGroupSizes(N);			
		assignGroups();
	}

	public int groupNum(int imachine, int level) {
		return groupNums[level][perm[imachine]];
	}

	public int posInGroup(int imachine, int level) {
		return posInGroups[level][perm[imachine]];
	}

	public int [] nodesInGroup(int imachine, int level) {
		int g = groupNums[level][perm[imachine]];
		int [] grp = nbrs[level][g];
		int [] out = new int[grp.length];
		for (int i = 0; i < grp.length; i++) {
			out[i] = invperm[grp[i]];
		}
		return out;
	}

	public int depth() {
		return groupNums.length;
	}
}