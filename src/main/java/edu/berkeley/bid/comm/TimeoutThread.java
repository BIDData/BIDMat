package edu.berkeley.bid.comm;

import java.util.concurrent.Future;

public class TimeoutThread implements Runnable {
		Future <?> [] futures;
		int mtime;
		String canceller;
		Machine machine;

		public TimeoutThread(Machine machine0, int mtime0, Future <?> [] futures0, String canceller0) {		
			futures = futures0;
			mtime = mtime0;
			canceller = canceller0;
			machine = machine0;
		}

		public void run() {
			try {
				Thread.sleep(mtime);

				for (int i = 0; i < futures.length; i++) {
					if (futures[i] != null) {
						if (machine.trace > 1 && machine.network != null) {
							synchronized (machine.network) {
								System.out.format("Machine %d cancelling %s future %d\n", machine.imachine, canceller, i);
							}
						}
						futures[i].cancel(true);
					}
				}
			} catch (InterruptedException e) {
			}
		}
	}