package BIDMat

import edu.berkeley.bid.comm.{IVec,LVec,Vec}
import edu.berkeley.bid.comm.AllReduceX
import scala.concurrent.future
import scala.concurrent.ExecutionContext.Implicits.global
import BIDMat.MatFunctions._
import BIDMat.SciFunctions._
import java.util.concurrent.CountDownLatch


object AllReduce {
  /*
   * Arguments: 
   *   n: array size
   *   dens: array density
   *   ks: comma-separated list of k values
   *   trace: 0,1,2 increasing amount of trace printed
   *   replicate: Integer replication factor. Normally either 1 (no replication) or 2.
   *   deads: comma-separated list of dead node ids
   */
  
  def main(args:Array[String]) = {
    val n = args(0).toInt
    val dens = args(1).toFloat
    val ks = args(2).split(",")
    val trace = if (args.length > 3) args(3).toInt else 0
    val replicate = if (args.length > 4) args(4).toInt else 1
    val deads = if (args.length > 5) args(5).split(",") else new Array[String](0)
    val deadnodes = IMat(deads.length, 1)
    for (i <- 0 until deads.length) deadnodes(i) = deads(i).toInt
    val allks = IMat(ks.length, 1)
    for (i <- 0 until ks.length) allks(i) = ks(i).toInt
    test1(n, dens, allks, trace, replicate, deadnodes)
  }
  
  def test1(n:Int, dens:Float, allks:IMat, trace:Int, replicate:Int, deadnodes:IMat) {
    val a = sprand(n, n, dens)
    val dd = spdiag(n)
    val ad = a + dd
    val b = rand(n, 1)
    makeSim(ad, b, allks, trace, replicate, deadnodes)
  }
  
  def makeSim(a:SMat, b:FMat, allks:IMat, trace:Int, replicate:Int, deadnodes:IMat) = {
    val M = allks.data.reduceLeft(_*_)
    val bufsize = 5 * (a.nnz / M) 
    val ioff = Mat.ioneBased
    val network = new AllReduceX(M*replicate)
    val rowvecs = new Array[IMat](M)
    val colvecs = new Array[IMat](M)
    val vvecs = new Array[Vec](M)
    val irows = new Array[IVec](M)
    val icols = new Array[IVec](M)
    val ivals = new Array[Vec](M)
    val retval = new Array[Vec](M*replicate)
    val smats = new Array[SMat](M)
    val (ii, jj, vv) = find3(a)
    val rr = min(M-1, IMat(rand(ii.nrows, 1)*M))
    val counts = accum(rr, 1, M)
    for (j <- 0 until replicate) {
    	for (i <- 0 until M) {
    		rowvecs(i) = IMat(counts(i), 1)
    		colvecs(i) = IMat(counts(i), 1)
    		vvecs(i) = new Vec(counts(i))
    		network.simNetwork(i+j*M) = new network.Machine(a.nrows, allks.data, i + j*M, M, bufsize, false, trace, replicate, null);
    	}
    }
    var i = 0
    counts.clear
    while (i < rr.length) {
      val iix = rr(i)
      val ic = counts(iix)
      rowvecs(iix)(ic) = ii(i)
      colvecs(iix)(ic) = jj(i)
      vvecs(iix).data(ic) = vv(i)
      counts(iix) = ic+1
      i += 1
    }
    var totvals = 0L
    for (i <- 0 until M) {
      val s = new SMat(a.nrows, a.ncols, rowvecs(i).length, SparseMat.incInds(rowvecs(i).data),
          SparseMat.compressInds(colvecs(i).data, a.ncols, new Array[Int](a.ncols+1), rowvecs(i).length), vvecs(i).data)
      val (ii1, jj1, vv1) = find3(sum(s,2))      
      irows(i) = new IVec(ii1.data)      
      ivals(i) = new Vec(vv1.data)
      icols(i) = new IVec(find(sum(s,1)).data)
      totvals += icols(i).size()
    }
    System.setProperty("actors.corePoolSize", "%d" format M*replicate)
    System.setProperty("actors.maxPoolSize", "%d" format M*replicate)
    val latch = new CountDownLatch(M*replicate)
    tic
    val nreps = 1
    	for (i <- 0 until M*replicate) {
//    	  if (irep == 17) network.simNetwork(i).trace=2; else network.simNetwork(i).trace=0;
    		future {        
    			if (deadnodes.length == 0 || sum(deadnodes == i).v == 0) {    // Simulate dead nodes
    				val i0 = i % M
    				network.simNetwork(i).config(irows(i0), icols(i0))
    				retval(i) = network.simNetwork(i).reduce(ivals(i0))
    				retval(i) = network.simNetwork(i).reduce(ivals(i0))
    			}
    			latch.countDown()
    		}
    	}
    	latch.await()
    network.stop
    val tt=toc
    println("Allreduce done, t=%4.3f msecs, BW=%4.3fGB/s" format (1000*tt/nreps, totvals*replicate*8*2/tt*nreps/1e9))
    val msum = new Vec(a.ncols)
    msum.clear
    for (i <- 0 until M) {
      ivals(i).addTo(msum, irows(i))
    }
    var nerrors = 0
    var nnodes = 0
    for (i <- 0 until M*replicate) {
    	if (retval(i) != null && retval(i).length > 0) {
    	  nnodes += 1
    		var j = 0 
    		val i0 = i % M
    		while (j < icols(i0).size()) {

    			val v1 = retval(i).data(j)
    			val v2 = msum.data(icols(i0).data(j))
    			if (math.abs(v1-v2)/math.max(1e-9, v2) > 1e-6) {
    				println("Bad value machine %d, pos %d, index %d, vals %f %f" format (i, j, icols(i0).data(j), v1, v2))
    				nerrors += 1
    			}
    			j += 1
    		}
      }
    }
    println("Checked %d nodes, %d errors" format (nnodes, nerrors))
  }
  
  
}
