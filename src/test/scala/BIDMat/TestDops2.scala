package BIDMat

import Mat._
import DMat._
import FMat._
import scala.compat.Platform._


object TestDops2 {
  def main(args: Array[String]): Unit = {
    val n = 50000
    val k = 20
    val l = 1
    val a = FMat(l,n)
    val b = FMat(n,k)
    val d = FMat(k,n)
    val e = FMat(n,l)
    val m = 30000
    val t0 = currentTime
    println("Starting up")
    for (i <- 0 until m) {
      val c = a * b
//      val c = d * e
    }
    val t1 = currentTime - t0
	println("time="+t1+"msec, gflops="+(2.0*k*n*l*m/t1/1e6))
  }
}

