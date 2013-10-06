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
    val allks = IMat(ks.length, 1)
    for (i <- 0 until ks.length) allks(i) = ks(i).toInt
    test1(n, dens, allks)
  }
  
  def test1(n:Int, dens:Float, allks:IMat) {
    val a = sprand(n, n, dens)
    val dd = spdiag(n)
    val ad = a + dd
    val b = rand(n, 1)
    makeSim(ad, b, allks)
  }
  
  def makeSim(a:SMat, b:FMat, allks:IMat) = {
    val M = allks.data.reduceLeft(_*_)
    val bufsize = 5 * (a.nnz / M) 
    val ioff = Mat.ioneBased
    val network = new AllReduce(M)
    val rowvecs = new Array[IMat](M)
    val colvecs = new Array[IMat](M)
    val vvecs = new Array[FMat](M)
    val irows = new Array[IVec](M)
    val icols = new Array[IVec](M)
    val smats = new Array[SMat](M)
    val (ii, jj, vv) = find3(a)
    val rr = IMat(rand(ii.nrows, 1)*M)
    val counts = accum(rr, 1, M)
    for (i <- 0 until M) {
      rowvecs(i) = IMat(counts(i), 1)
      colvecs(i) = IMat(counts(i), 1)
      vvecs(i) = FMat(counts(i), 1)
      network.simNetwork(i) = new network.Machine(a.nrows, allks.data, i, M, bufsize, true)
    }
    var i = 0
    counts.clear
    while (i < rr.length) {
      val iix = rr(i)
      val ic = counts(iix)
      rowvecs(iix)(ic) = ii(i)
      colvecs(iix)(ic) = jj(i)
      vvecs(iix)(ic) = vv(i)
      counts(iix) = ic+1
      i += 1
    }
    for (i <- 0 until M) {
      val s = new SMat(a.nrows, a.ncols, rowvecs(i).length, SparseMat.incInds(rowvecs(i).data),
          SparseMat.compressInds(colvecs(i).data, a.ncols, new Array[Int](a.ncols+1), rowvecs(i).length), vvecs(i).data)
      irows(i) = new IVec(find(sum(s,2)).data)
      icols(i) = new IVec(find(sum(s,1)).data)      
    }
    val latch = new CountDownLatch(M)
    for (i <- 0 until M) {
      actor {
        network.simNetwork(i).config(irows(i), icols(i))
        network.simNetwork(i).reduce(new Vec(vvecs(i).data))
        latch.countDown()
      }
    }
    latch.await();
    println("All done")
  }
  
  
}