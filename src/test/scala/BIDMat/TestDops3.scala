package BIDMat

import Mat._
import DMat._
import FMat._
import scala.compat.Platform._


object TestDops3 {
  def main(args: Array[String]): Unit = {
    val n = 50000
    val k = 10
    val l = 1
    val a = FMat(n,k)
//    val a2 = FMat(k,n)
    val b = FMat(l,n)
    val d = FMat(k,l)
    val m = 30000
    val t0 = currentTime
    println("Starting up")
    for (i <- 0 until m) {
//      val c = b * a
//      val e = a * d
      val c = a t
    }
    val t1 = currentTime - t0
	println("time="+t1+"msec, gflops="+(2.0*k*n*l*m/t1/1e6))
  }
}

