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
    val (ii, jj, vv) = find3(a)
    val rr = IMat(rand(ii.nrows, 1)*M)
    for (i <- 0 until M) {
      val inds = find(rr == i)
      rowvecs(i) = ii(inds)
      colvecs(i) = jj(inds)
      vvecs(i) = vv(inds)
      network.simNetwork(i) = new network.Machine(a.nrows, allks.data, i, M, bufsize, false)
    }
    
    for (i <- 0 until M) {
      actor {
        network.simNetwork(i).config(new IVec(rowvecs(i).data), new IVec(colvecs(i).data))
        network.simNetwork(i).reduce(new Vec(vvecs(i).data))
      }
    }
  }
  
  
}