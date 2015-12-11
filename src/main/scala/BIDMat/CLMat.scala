/**
 * OpenCL, single precision, dense matrix
 */

package BIDMat

import org.jocl.CL._
import org.jocl._

case class CLMat(nr: Int, nc: Int, var data: cl_mem, realsize: Long) extends Mat(nr, nc) {

  override def mytype = "CLMat"

}

object CLMat {

  

}
