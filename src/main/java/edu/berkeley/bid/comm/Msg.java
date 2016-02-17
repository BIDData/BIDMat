package edu.berkeley.bid.comm;

public class Msg {
	byte [] buf;
	int size;
	int sender;
	int receiver;
	int tag;

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
}