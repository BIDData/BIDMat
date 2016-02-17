
package edu.berkeley.bid.comm;

// Third version of sparse allreduce
// Includes support for matrices, feature range limits, long feature indices
//
// 

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.Arrays;
import java.lang.Math;
import java.io.*;
import java.nio.*;
import java.net.*;
//import mpi.*;


public class AllReduceY {
	
	public Machine [] simNetwork = null;

	public AllReduceY(int M) {
		simNetwork = new Machine[M];
	}

	public void stop() {
		if (simNetwork != null) {
			for (int i = 0; i < simNetwork.length; i++) simNetwork[i].stop();
		}
	}

}
