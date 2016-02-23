package BIDMat

import edu.berkeley.bid.comm._
import scala.collection.parallel._
import BIDMat.MatFunctions._
import BIDMat.SciFunctions._
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;


class AllReducer(val M:Int, val F:Int, val nnz:Int) {
  /*
   * Arguments: 
   *   M: number of machines
   *   F: number of features
   *   nnz: number of non-zero features (expected) per node
   *   stride: number of rows in the matrices to reduce
   *   trace: 0,1,2 increasing amount of trace printed
   *   replicate: Integer replication factor. Normally either 1 (no replication) or 2.
   *   deads: comma-separated list of dead node ids
   */
  
	var configDir = "/code/BIDMat/data/bestAllreduce/";
	
	var gmods:IMat = null;
	var gridmachines:IMat = null;        
	var groups:Groups = null;
	var network:Network = null;

	var nodeData:Array[FMat] = null;
	var nodeInds:Array[IMat] = null;
	var reducedData:Array[FMat] = null;
	var total:FMat = null;
	var netExecutor:ExecutorService = null;
	var doSim:Boolean = true;
	var configTimeout = 1000;
	var reduceTimeout = 1000;
	var doConfigReduce = false;
	var trace = 0;
	var replicate:Int = 1;
	var deadnodes:IMat = IMat(0,1);
  
  def runSim(stride:Int) = {
    if (gmods.asInstanceOf[AnyRef] == null) {
    	val clengths = loadIMat(configDir + "dims.imat.lz4");
    	val allgmods = loadIMat(configDir + "gmods.imat.lz4");
    	val allmachinecodes = loadIMat(configDir + "machines.imat.lz4");
    	gmods = allgmods(0->clengths(M-1), M-1);
    	gridmachines = allmachinecodes(0->M, M-1);
    }         
    groups = new Groups(M, gmods.data, gridmachines.data, 1000);
    network = new Network(M);
    netExecutor = Executors.newFixedThreadPool(M+2);

    nodeData = new Array[FMat](M);
    nodeInds = new Array[IMat](M);
    reducedData = new Array[FMat](M);
    
    total = zeros(stride, F);
    var totvals = 0L;

    for (i <- 0 until M) {                                         // Make some random, power-law data
    	val rv = rand(1, F);
    	val thresh = nnz/math.log(nnz)/row(1 to F);
    	nodeInds(i) = find(rv < thresh);
    	nodeData(i) = rand(stride, nodeInds(i).length);
    	val bufsize = (1.5 * nodeData(0).length).toInt;
    	network.machines(i) = new Machine(network, groups, i, M, bufsize, doSim, trace, replicate, null);
    	network.machines(i).configTimeout = configTimeout;
    	network.machines(i).reduceTimeout = reduceTimeout ;
    	totvals += stride * nodeInds(i).length;
    }

    val nreps =1;
    tic;
    val futures = (0 until M).toArray.map((i:Int) =>
      netExecutor.submit(new Runnable {def run:Unit = {
        	if (deadnodes.length == 0 || sum(deadnodes == i).v == 0) {    // Simulate dead nodes
        	  val result = if (doConfigReduce) {
        	  	network.machines(i).configReduce(nodeInds(i).data, nodeInds(i).data, nodeData(i).data, stride);  
        	  } else {
        	  	network.machines(i).config(nodeInds(i).data, nodeInds(i).data);
        	  	network.machines(i).reduce(nodeData(i).data, stride);
        	  }
        		reducedData(i) = new FMat(stride, nodeInds(i).length, result);
        	}
        }
      }));
    	//    	  if (irep == 17) network.simNetwork(i).trace=2; else network.simNetwork(i).trace=0;     
    for (i <- 0 until M) futures(i).get();                              // Block until all threads are done

    network.stop;
    netExecutor.shutdownNow();
    val tt= toc;
    println("Allreduce done, t=%4.3f msecs, BW=%4.3fGB/s" format (1000*tt/nreps, totvals*replicate*8*2/tt*nreps/1e9));

    // Compute ground truth
    for (i <- 0 until M) {
      total(?, nodeInds(i)) = total(?, nodeInds(i)) + nodeData(i);
    }
    var nerrors = 0;
    for (i <- 0 until M) {
      val diff = abs(reducedData(i) - total(?, nodeInds(i)));
      nerrors += sum(sum(diff > 1e-3)).v.toInt;
    }
    println("Checked %d nodes, %d errors" format (M, nerrors));
  }
  
}