package edu.berkeley.bid.comm;

/* Class containing sparse, sorted integer indices */
/* Supports merging and finding mappings from one set of indices to another */

import java.util.List;
import java.util.LinkedList;
import java.lang.StringBuffer;


public class IVec {
	public int[] data;
	public int length;
	
	public IVec(int n) {
		data = new int[n];
		length = n;
	}
	
	public IVec(int [] d) {
		data = d;
		length = d.length;
	} 
	
	static public IVec irow(int... d) {
		IVec out = new IVec(d.length);
		out.data = d;
		return out;
	}
	
	public int size() {
		return data.length;
	}
	
	@Override public String toString() {
		StringBuffer out = new StringBuffer();
		out.append("IVec (");
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
	
	static public LinkedList<IVec> partition(IVec inds, IVec part) {
		LinkedList<IVec> out = new LinkedList<IVec>();
		int [] id = inds.data;
		int here = 0;
		int newsize = 0;
		for (int i = 0; i < part.length; i++) {
			int end = part.data[i];
			while (here < id.length && id[here] < end) {
				newsize++;
				here++;
			}
			out.add(new IVec(newsize));
			newsize = 0;
		}
		here = 0;
		int newi = 0;
		int i = 0;
		for (IVec ovec : out) {
			int [] ovd = ovec.data;
			int end = part.data[i];
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
	/* Output: k+1 IVecs. The first is the merged set of indices, call it "master" */
	/* the other k outputs are the mappings from each input set to master */
	
	static public LinkedList<IVec> mergeAndMap(List<IVec> indlist) {
		LinkedList<IVec> out = new LinkedList<IVec>();
		IVec master = merge(indlist);
		out.add(master);
		for (IVec iv : indlist) {
			out.add(mapInds(iv, master));
		}		
		return out;		
	}
	
	static public IVec merge(List<IVec> indlist) {
		int treesize = 1;
		int height = 1;
		while (treesize < indlist.size()) {
			height += 1;
			treesize *= 2;
		}
		IVec [] tree = new IVec[height];
		for (IVec iv : indlist) {
			treeAdd(iv, tree);
		}
		IVec master = treeFlush(tree);
		return master;		
	}
	
	static public IVec mapInds(IVec a, IVec b)  { 
		IVec out = new IVec(a.length);
		int [] ad = a.data;
		int [] bd = b.data;
		int [] od = out.data;
		for (int i = 0; i < out.length; i++) {od[i] = -1;}
		int i = 0;
		int ito = 0;
		while (i < a.length && ito < b.length) { 
		  int xx = ad[i] - bd[ito];
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
	
  static public void treeAdd(IVec x, IVec [] tree) {
    if (x != null) {
    	IVec dd = x;
    	int j = 0;
    	while (tree[j] != null) {
    		dd = merge2(tree[j], dd);
    		tree[j] = null;
    		j += 1;
    	}
    	tree[j] = dd;
    }
  }
  
  static public IVec treeFlush(IVec [] tree) {
    int j = 0;
    IVec dd = null;
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
	
	static public IVec merge2(IVec a, IVec b) {
    int i = 0;
    int j = 0;
    int [] ad = a.data;
    int [] bd = b.data;
    int nout = 0;
    while (i < a.length && j < b.length) {
    	int comp = ad[i] - bd[j];
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
    IVec out = new IVec(nout);
    int [] od = out.data;
    i = 0;
    j = 0;
    nout = 0;
    while (i < a.length && j < b.length) {
    	int comp = ad[i] - bd[j];
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