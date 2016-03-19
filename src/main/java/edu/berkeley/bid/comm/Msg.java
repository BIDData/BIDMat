package edu.berkeley.bid.comm;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

public class Msg {
	byte [] buf;
	public int size;
	public int sender;
	public int receiver;
	public int tag;

	public Msg(int size0, int sender0, int receiver0, int tag0) {
		buf = new byte[4*size0];
		size = size0;
		sender = sender0;
		receiver = receiver0;
		tag = tag0;
	}

	public Msg(byte [] inbuf, int size0, int sender0, int receiver0, int tag0) {
		buf = new byte[4*size0];
		System.arraycopy(inbuf, 0, buf, 0, 4*size0);
		size = size0;
		sender = sender0;
		receiver = receiver0;
		tag = tag0;
	}
	
	public static String printStack(Exception e) { 
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		PrintStream ps = new PrintStream(baos);
		e.printStackTrace(ps);
		String str = baos.toString();
		ps.close();
		return str;
	}
}
