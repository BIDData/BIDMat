package BIDMat


import DMat._
import IMat._
import FMat._
import scala.compat.Platform._


object TestDops {
  def main(args: Array[String]): Unit = {
	 val n = 2000
	 val a = IMat(n,n)
	 val b = IMat(n,n)
	 val t0 = currentTime
	 val m = 1000
	 println("starting up")
	 for (i <- 0 until m) {
		 val c = a + b
	 }
	 val t1 = currentTime - t0
	 println("time="+t1+" msec, gflops="+(n.doubleValue*n*m/t1/1e6))
  }
}

