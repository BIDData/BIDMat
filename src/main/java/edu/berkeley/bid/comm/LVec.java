package edu.berkeley.bid.comm;

/* Class containing sparse, sorted long indices */
/* Supports merging and finding mappings from one set of indices to another */

import java.util.List;
import java.util.LinkedList;
import java.lang.StringBuffer;


public class LVec {
	public long[] data;
	public int length;
	
	public LVec(int n) {
		data = new long[n];
		length = n;
	}
	
	public LVec(long [] d) {
		data = d;
		length = d.length;
	} 
	
	static public LVec lrow(long... d) {
		LVec out = new LVec(d);
		return out;
	}
	
	public int size() {
		return data.length;
	}
	
	@Override public String toString() {
		StringBuffer out = new StringBuffer();
		out.append("LVec (");
		int lim = Math.min(20, length);
		for (int i = 0; i < lim; i++) {
			out.append(data[i]);
			if (i < lim-1) out.append(", ");
		}
		out.append(")");
		return out.toString();
	}
	
	/* Partition the input indices into ranges set by the elements of part */
	/* there are part.length ranges: (>-infty,...,<part(0)), (>=part(0) && <part(1)), ...(>=part(n-2) && <part(n-1) */
	
	static public LinkedList<LVec> partition(LVec inds, LVec part) {
		LinkedList<LVec> out = new LinkedList<LVec>();
		long [] id = inds.data;
		int here = 0;
		int newsize = 0;
		for (int i = 0; i < part.length; i++) {
			long end = part.data[i];
			while (here < id.length && id[here] < end) {
				newsize++;
				here++;
			}
			out.add(new LVec(newsize));
			newsize = 0;
		}
		here = 0;
		int newi = 0;
		int i = 0;
		for (LVec ovec : out) {
			long [] ovd = ovec.data;
			long end = part.data[i];
			while (here < id.length && id[here] < end) {
				ovd[newi] = id[here];
				newi++;
				here++;
			}
			newi = 0;
			i++;
		}		
		return out;
	}
	
	
	/* Input: list of k ranges of indices */
	/* Output: k+1 objects. The first LVec is the merged set of indices, call it "master" */
	/* the other k outputs (IVec) are the mappings from each input set to master */
	
	static public LinkedList<Object> mergeAndMap(List<LVec> indlist) {
		LinkedList<Object> out = new LinkedList<Object>();
		LVec master = merge(indlist);
		out.add(master);
		for (LVec iv : indlist) {
			out.add(mapInds(iv, master));
		}		
		return out;		
	}
	
	static public LVec merge(List<LVec> indlist) {
		int treesize = 1;
		int height = 1;
		while (treesize < indlist.size()) {
			height += 1;
			treesize *= 2;
		}
		LVec [] tree = new LVec[height];
		for (LVec iv : indlist) {
			treeAdd(iv, tree);
		}
		LVec master = treeFlush(tree);
		return master;		
	}
	
	// compute a map from sorted indices a to sorted indices b
	// so if c = mapInds(a,b) then a[i] = b[c[i]];
	
	static public IVec mapInds(LVec a, LVec b)  { 
		IVec out = new IVec(a.length);
		long [] ad = a.data;
		long [] bd = b.data;
		int [] od = out.data;
		for (int i = 0; i < out.length; i++) {od[i] = -1;}
		int i = 0;
		int ito = 0;
		while (i < a.length && ito < b.length) { 
		  long xx = ad[i] - bd[ito];
		  if (xx <= 0) {
		  	if (xx == 0) {
		  		od[i] = ito;
		  		ito += 1;
		  	}
		  	i += 1;
		  } else {
		  	ito += 1;
		  } 
		}
		return out;
	}
	
  static public void treeAdd(LVec x, LVec [] tree) {
    if (x != null) {
    	LVec dd = x;
    	int j = 0;
    	while (tree[j] != null) {
    		dd = merge2(tree[j], dd);
    		tree[j] = null;
    		j += 1;
    	}
    	tree[j] = dd;
    }
  }
  
  static public LVec treeFlush(LVec [] tree) {
    int j = 0;
    LVec dd = null;
    while (j < tree.length) {
    	if (tree[j] != null) {
    	  if (dd != null) {
    	  	dd = merge2(tree[j], dd);
    	  } else {
    	    dd = tree[j];
    	  }
    	  tree[j] = null;
    	}
    	j += 1;
    }
    return dd;
  }
	
	static public LVec merge2(LVec a, LVec b) {
    int i = 0;
    int j = 0;
    long [] ad = a.data;
    long [] bd = b.data;
    int nout = 0;
    while (i < a.length && j < b.length) {
    	long comp = ad[i] - bd[j];
      if (comp <= 0) {
        i += 1;
      } 
      if (comp >= 0) {
        j += 1;
      } 
      nout += 1;
    }
    if (i < a.length) {
      nout += a.length - i;
    }
    if (j < b.length) {
      nout += b.length - j;
    }
    LVec out = new LVec(nout);
    long [] od = out.data;
    i = 0;
    j = 0;
    nout = 0;
    while (i < a.length && j < b.length) {
    	long comp = ad[i] - bd[j];
      if (comp <= 0) {
      	od[nout] = ad[i];
        i += 1;
      }
      if (comp >= 0) {
      	if (comp > 0) {
      		od[nout] = bd[j];
      	}
      	j += 1;
      } 
      nout += 1;
    }
    while (i < a.length) {
    	od[nout] = ad[i];
      i += 1;
      nout += 1;
    }
    while (j < b.length) {
    	od[nout] = bd[j];
      j += 1;
      nout += 1;
    }
    return out;
  }	
	
}