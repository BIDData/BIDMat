package BIDMat

import edu.berkeley.bid.comm.{IVec,LVec,Vec}
import edu.berkeley.bid.comm.AllReduce
import scala.actors.Actor._
import BIDMat.MatFunctions._
import BIDMat.SciFunctions._
import java.util.concurrent.CountDownLatch


object AllReducer {
  
  def main(args:Array[String]) = {
    val n = args(0).toInt
    val dens = args(1).toFloat
    val ks = args(2).split(",")
    val trace = if (args.length > 3) args(3).toInt else 0
    val allks = IMat(ks.length, 1)
    for (i <- 0 until ks.length) allks(i) = ks(i).toInt
    test1(n, dens, allks, trace)
  }
  
  def test1(n:Int, dens:Float, allks:IMat, trace:Int) {
    val a = sprand(n, n, dens)
    val dd = spdiag(n)
    val ad = a + dd
    val b = rand(n, 1)
    makeSim(ad, b, allks, trace)
  }
  
  def makeSim(a:SMat, b:FMat, allks:IMat, trace:Int) = {
    val M = allks.data.reduceLeft(_*_)
    val bufsize = 5 * (a.nnz / M) 
    val ioff = Mat.ioneBased
    val network = new AllReduce(M)
    val rowvecs = new Array[IMat](M)
    val colvecs = new Array[IMat](M)
    val vvecs = new Array[Vec](M)
    val irows = new Array[IVec](M)
    val icols = new Array[IVec](M)
    val ivals = new Array[Vec](M)
    val retval = new Array[Vec](M)
    val smats = new Array[SMat](M)
    val (ii, jj, vv) = find3(a)
    val rr = IMat(rand(ii.nrows, 1)*M)
    val counts = accum(rr, 1, M)
    for (i <- 0 until M) {
      rowvecs(i) = IMat(counts(i), 1)
      colvecs(i) = IMat(counts(i), 1)
      vvecs(i) = new Vec(counts(i))
      network.simNetwork(i) = new network.Machine(a.nrows, allks.data, i, M, bufsize, true, trace)
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
    for (i <- 0 until M) {
      val s = new SMat(a.nrows, a.ncols, rowvecs(i).length, SparseMat.incInds(rowvecs(i).data),
          SparseMat.compressInds(colvecs(i).data, a.ncols, new Array[Int](a.ncols+1), rowvecs(i).length), vvecs(i).data)
      val (ii1, jj1, vv1) = find3(sum(s,2))      
      irows(i) = new IVec(ii1.data)      
      ivals(i) = new Vec(vv1.data)
      icols(i) = new IVec(find(sum(s,1)).data)
    }
    System.setProperty("actors.corePoolSize", "%d" format M)
    System.setProperty("actors.maxPoolSize", "%d" format M)
    val latch = new CountDownLatch(M)
    for (i <- 0 until M) {
      actor {
        network.simNetwork(i).config(irows(i), icols(i))
        retval(i) = network.simNetwork(i).reduce(ivals(i))
        latch.countDown()
      }
    }
    latch.await();
    println("Allreduce done")
    
    val msum = new Vec(a.ncols)
    msum.clear
    for (i <- 0 until M) {
      ivals(i).addTo(msum, irows(i))
    }
    var nerrors = 0
    for (i <- 0 until M) {
      var j = 0 
      while (j < icols(i).size()) {
        val v1 = retval(i).data(j)
        val v2 = msum.data(icols(i).data(j))
        if (Math.abs(v1-v2)/Math.max(1e-9, v2) > 1e-6) {
          println("Bad value machine %d, pos %d, index %d, vals %f %f" format (i, j, icols(i).data(j), v1, v2))
          nerrors += 1
        }
        j += 1
      }
    }
    println("Checking done, %d errors" format nerrors)
  }
  
  
}