//-*-coding:utf-8-*-
package BIDMat
import MatFunctions._
import edu.berkeley.bid.CBLAS._
import java.util.Arrays
import java.util.concurrent.atomic._
import scala.concurrent.future
import scala.concurrent.ExecutionContext.Implicits.global


case class F3D(val d1:Int, val d2:Int, val d3:Int, val data:Array[Float]) extends A3D(d1, d2, d3) { 

  def apply(indx:Int):Float = { 
    if (indx < length) { 
      data(indx)
    } else { 
      throw new RuntimeException("F3D index out of range")
    }
  }

  def apply(x:Int, y:Int, z:Int):Float = { 
    val indx = x + dim1 * (y + dim2 * z )
    if (indx < length) { 
      data(indx)
    } else { 
      throw new RuntimeException("F3D indices out of range")
    }
  }
  
  def apply(ii:IMat, jj:IMat, kk:IMat):F3D = {
    val nx = ii match {case aa:MatrixWildcard => dim1; case _ => ii.length}
    val ny = jj match {case aa:MatrixWildcard => dim2; case _ => jj.length}
    val nz = kk match {case aa:MatrixWildcard => dim3; case _ => kk.length}
    val out = F3D(nx, ny, nz)
    var k = 0
    while (k < nz) {
      val kval = kk match {
      case aa:MatrixWildcard => k
      case _ => kk.data(k)
      }
      var j = 0
      while (j < ny) {
        val jval = jj match {
          case aa:MatrixWildcard => j
          case _ => jj.data(j)
        }
        val outbase = nx * (j + ny * k)
        val inbase = dim1 * (jval + dim2 * kval)
        var i = 0
        ii match {
          case aa:MatrixWildcard => {
        	while (i < nz) {
        	  out.data(i + outbase) = data(i + inbase) 
        	  i += 1
        	}
          }
          case _ => {
        	while (i < nz) {
        	  out.data(i + outbase) = data(ii.data(i) + inbase)  
        	  i += 1
        	}
          }
        }
        j += 1
      }
      k += 1
    }
    out
  }
  
  def update(indx:Int, v:Float):F3D = { 
    if (indx < 0 || indx < length) { 
      data(indx) = v
    } else { 
      throw new RuntimeException("F3D index out of range")
    }
    this
  }

  def update(x:Int, y:Int, z:Int, v:Float):F3D = { 
    val indx = x + dim1 * (y + dim2 * z)
    if (indx < 0 || indx < length) { 
      data(indx) = v
    } else { 
      throw new RuntimeException("F3D indices out of range")
    }
    this
  }
  
  override def update(ii:IMat, jj:IMat, kk:IMat, vv:A3D):F3D = {
    vv match {
      case uu:F3D => update(ii, jj, kk, uu)
      case _ => throw new RuntimeException("F3D update illegal RHS")
    }
  }
  
  def update(ii:IMat, jj:IMat, kk:IMat, vv:F3D):F3D = {
    val nx = ii match {case aa:MatrixWildcard => dim1; case _ => ii.length}
    val ny = jj match {case aa:MatrixWildcard => dim2; case _ => jj.length}
    val nz = kk match {case aa:MatrixWildcard => dim3; case _ => kk.length}
    if (nx != vv.dim1 || ny != vv.dim2 || nz != vv.dim3) {
      throw new RuntimeException("F3D update dimensions mismatch")
    }
    var k = 0
    while (k < nz) {
      val kval = kk match {
      case aa:MatrixWildcard => k
      case _ => kk.data(k)
      }
      var j = 0
      while (j < ny) {
        val jval = jj match {
          case aa:MatrixWildcard => j
          case _ => jj.data(j)
        }
        val inbase = nx * (j + ny * k)
        val outbase = dim1 * (jval + dim2 * kval)
        var i = 0
        ii match {
          case aa:MatrixWildcard => {
        	while (i < nz) {
        	  data(i + outbase) = vv.data(i + inbase) 
        	  i += 1
        	}
          }
          case _ => {
        	while (i < nz) {
        	  data(ii.data(i) + outbase) = vv.data(i + inbase) 
        	  i += 1
        	}
          }
        }
        j += 1
      }
      k += 1
    }
    this
  }
  
  def update(ii:IMat, jj:IMat, kk:IMat, v:Float):F3D = {
    val nx = ii match {case aa:MatrixWildcard => dim1; case _ => ii.length}
    val ny = jj match {case aa:MatrixWildcard => dim2; case _ => jj.length}
    val nz = kk match {case aa:MatrixWildcard => dim3; case _ => kk.length}
    var k = 0
    while (k < nz) {
      val kval = kk match {
      case aa:MatrixWildcard => k
      case _ => kk.data(k)
      }
      var j = 0
      while (j < ny) {
        val jval = jj match {
          case aa:MatrixWildcard => j
          case _ => jj.data(j)
        }
        val outbase = dim1 * (jval + dim2 * kval)
        var i = 0
        ii match {
          case aa:MatrixWildcard => {
        	while (i < nz) {
        	  data(i + outbase) = v 
        	  i += 1
        	}
          }
          case _ => {
        	while (i < nz) {
        	  data(ii.data(i) + outbase) = v
        	  i += 1
        	}
          }
        }
        j += 1
      }
      k += 1
    }
    this
  }


  def permute(p1:Int, p2:Int, p3:Int):F3D = permute(irow(p1,p2,p3))

  def permute(perm0:IMat):F3D = { 
    if (perm0.length != 3) { 
      throw new RuntimeException("F3D permute bad permutation ")
    }
    val perm = perm0.copy
    val dims = irow(dim1, dim2, dim3)
    val iperm = invperm(perm)
    var out = F3D(dims(iperm(0)), dims(iperm(1)), dims(iperm(2)))
    var out2 = F3D(dims(iperm(0)), dims(iperm(1)), dims(iperm(2)))
    System.arraycopy(data, 0, out.data, 0, length)
    for (i <- 2 until 0 by -1) { 
      if (perm(i) != i) { 
        val (d1, d2, d3) = A3D.getDims(i, perm, dims)
        if (d1 > 1 && d2 > 1) { 
          spermute(d1, d2, d3, out.data, out2.data)
          val tmp = out2
          out2 = out
          out = tmp
        }
        A3D.rotate(i, perm, dims)
      } 
    }
    out
  }

}

object F3D {
  
  def apply(dim1:Int, dim2:Int, dim3:Int):F3D = new F3D(dim1, dim2, dim3, new Array[Float](dim1*dim2*dim3))

  def apply(dim1:Int, dim2:Int, dim3:Int, f:FMat):F3D = { 
    if (dim1*dim2*dim3 == f.length) { 
      new F3D(dim1, dim2, dim3, f.data)
    } else { 
      throw new RuntimeException("F3D input matrix size mismatch")
    }
  } 
  
}






