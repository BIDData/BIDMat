package BIDMat

import edu.berkeley.bid.comm.{IVec,LVec,Vec}
import edu.berkeley.bid.comm._
import scala.collection.parallel._
import BIDMat.MatFunctions._
import BIDMat.SciFunctions._
import java.util.concurrent.CountDownLatch

class AllReducer {
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
  
	val configDir = "/code/BIDMat/data/bestAllreduce/";
	
	var gmods:IMat = null;
	var gridmachines:IMat = null;        
	var groups:Groups = null;
	var network:edu.berkeley.bid.comm.AllReduceY = null;

	var nodeData:Array[FMat] = null;
	var nodeInds:Array[IMat] = null;
	var reducedData:Array[FMat] = null;
	var total:FMat = null;
  
  def main(args:Array[String]) = {
    val M = args(0).toInt                               
    val F = args(1).toInt
    val nnz = args(2).toInt
    val stride = args(3).toInt
    val trace = if (args.length > 4) args(4).toInt else 0
    val replicate = if (args.length > 5) args(5).toInt else 1
    val deads = if (args.length > 6) args(6).split(",") else new Array[String](0)
    val deadnodes = IMat(deads.length, 1)
    for (i <- 0 until deads.length) deadnodes(i) = deads(i).toInt

    makeSim(M, F, nnz, stride, trace, replicate, deadnodes);
  }
  
  def makeSim(M:Int, F:Int, nnz:Int, stride:Int, trace:Int = 0, replicate:Int = 1, deadnodes:IMat = IMat(0,1)) = {
    val clengths = loadIMat(configDir + "dims.imat.lz4");
    val allgmods = loadIMat(configDir + "gmods.imat.lz4");
    val allmachines = loadIMat(configDir + "machines.imat.lz4");
    gmods = allgmods(0->clengths(M-1), M-1);
    gridmachines = allmachines(0->M, M-1);        
    groups = new Groups(M, gmods.data, gridmachines.data, 1000);
    network = new AllReduceY(M);

    nodeData = new Array[FMat](M);
    nodeInds = new Array[IMat](M);
    reducedData = new Array[FMat](M);
    
    total = zeros(stride, F);
    var totvals = 0L;

    for (i <- 0 until M) {
    	val rv = rand(1, F);
    	val thresh = nnz/math.log(nnz)/row(1 to F);
    	nodeInds(i) = find(rv > thresh);
    	nodeData(i) = rand(stride, nodeInds(i).length);
    	val bufsize = (1.5 * nodeData(0).length).toInt;
    	network.machines(i) = new Machine(network, groups, i, M, bufsize, false, trace, replicate, null);
    	totvals += stride * nodeInds(i).length;
    }

    System.setProperty("actors.corePoolSize", "%d" format M*replicate);
    System.setProperty("actors.maxPoolSize", "%d" format M*replicate);
    val nreps =1;
    tic;
    val par = (0 until M).par.map((i:Int) =>
    	//    	  if (irep == 17) network.simNetwork(i).trace=2; else network.simNetwork(i).trace=0;     
      if (deadnodes.length == 0 || sum(deadnodes == i).v == 0) {    // Simulate dead nodes
    	  network.machines(i).config(nodeInds(i).data, nodeInds(i).data);
    	  val result = network.machines(i).reduce(nodeData(i).data, stride);
      	reducedData(i) = new FMat(stride, nodeInds(i).length, result);
      });

    network.stop;
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