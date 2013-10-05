package BIDMat

import edu.berkeley.bid.comm.{IVec,LVec,Vec}
import edu.berkeley.bid.comm.AllReduce
import scala.actors.Actor._
import BIDMat.MatFunctions._
import BIDMat.SciFunctions._



object AllReducer {
  
  def makeSim(a:SMat, b:FMat, allks:IMat) = {
    val M = allks.data.reduceLeft(_*_)
    val bufsize = 1024*1024;
    val network = new AllReduce(M)
    val rowvecs = new Array[IMat](M)
    val colvecs = new Array[IMat](M)
    val vvecs = new Array[FMat](M)
    val smats = new Array[SMat](M)
    val (ii, jj, vv) = find3(a)
    val rr = IMat(rand(ii.nrows, 1)*M)
    val counts = accum(rr, 1, M)
    for (i <- 0 until M) {
      rowvecs(i) = IMat(counts(i), 1)
      colvecs(i) = IMat(counts(i), 1)
      vvecs(i) = FMat(counts(i), 1)
      network.simNetwork(i) = new network.Machine(a.nrows, allks.data, i, M, bufsize, false)
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
    	smats(i) = new SMat(a.nrows, a.ncols, rowvecs(i).length, SparseMat.incInds(rowvecs(i).data), 
    			SparseMat.compressInds(new Array[Int](a.ncols+1), a.ncols, colvecs(i).data, rowvecs(i).length), vvecs(i).data)
    }
    for (i <- 0 until M) {
      actor {
        network.simNetwork(i).config(new IVec(rowvecs(i).data), new IVec(colvecs(i).data))
        network.simNetwork(i).reduce(new Vec(vvecs(i).data))
      }
    }
  }
  
  
}