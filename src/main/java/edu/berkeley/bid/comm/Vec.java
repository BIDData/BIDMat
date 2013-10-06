package edu.berkeley.bid.comm;

import java.util.List;
import java.util.LinkedList;
import java.lang.StringBuffer;


public class Vec {
	public float[] data;
	public int length;
	
	public Vec(int n) {
		data = new float[n];
		length = n;
	}
	
	public Vec(float [] d) {
		data = d;
		length = d.length;
	} 
	
	static public Vec row(float... d) {
		Vec out = new Vec(d.length);
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
	
	public void clear() {
		for (int i = 0; i < length; i++) {
			data[i] = 0;
		}
	}
	
	public int nnz() {
		int nn = 0;
		for (int i = 0; i < length; i++) {
			if (data[i] != 0) nn++;
		}
		return nn;
	}
	
	public void addTo(Vec b, IVec map) {
		if (length != map.length) {
			throw new RuntimeException("addTo dimensions mismatch "+length+" "+map.length);
		}
		float [] ad = data;
		float [] bd = b.data;
		int [] md = map.data;
		for (int i = 0; i < length; i++) {
			bd[md[i]] += ad[i];
		}
	}
	
	public Vec mapFrom(IVec map) {
		Vec out = new Vec(map.size());
		for (int i = 0; i < out.size(); i++) {
			out.data[i] = data[map.data[i]];
		}
		return out;
	}
	
}