
package edu.berkeley.bid.comm;

// Third version of sparse allreduce
// Includes support for matrices, feature range limits, long feature indices
//
// 


//import mpi.*;


public class AllReduceY {
	
	public Machine [] machines = null;

	public AllReduceY(int M) {
		machines = new Machine[M];
	}

	public void stop() {
		if (machines != null) {
			for (int i = 0; i < machines.length; i++) machines[i].stop();
		}
	}

}
