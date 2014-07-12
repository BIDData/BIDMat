//-*-coding:utf-8-*-
package BIDMat
import MatFunctions._
import edu.berkeley.bid.CBLAS._
import java.util.Arrays
import java.util.concurrent.atomic._
import scala.concurrent.future
import scala.concurrent.ExecutionContext.Implicits.global


case class A3D(val dim1:Int, val dim2:Int, val dim3:Int, val data:Array[Float]) { 

  def length = dim1*dim2*dim3

  def size() = length;

  private var _GUID = Mat.myrand.nextLong
  
  def setGUID(v:Long):Unit = {_GUID = v}
  
  def GUID:Long = _GUID

  def apply(indx:Int):Float = { 
    if (indx < length) { 
      data(indx)
    } else { 
      throw new RuntimeException("A3D index out of range")
    }
  }

  def apply(x:Int, y:Int, z:Int):Float = { 
    val indx = x + dim1 * (y + dim2 * z )
    if (indx < length) { 
      data(indx)
    } else { 
      throw new RuntimeException("A3D indices out of range")
    }
  }
  
  def apply(ii:IMat, jj:IMat, kk:IMat):A3D = {
    val nx = ii match {case aa:MatrixWildcard => dim1; case _ => ii.length}
    val ny = jj match {case aa:MatrixWildcard => dim2; case _ => jj.length}
    val nz = kk match {case aa:MatrixWildcard => dim3; case _ => kk.length}
    val out = A3D(nx, ny, nz)
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
  
  def update(ii:IMat, jj:IMat, kk:IMat, vv:A3D):A3D = {
    val nx = ii match {case aa:MatrixWildcard => dim1; case _ => ii.length}
    val ny = jj match {case aa:MatrixWildcard => dim2; case _ => jj.length}
    val nz = kk match {case aa:MatrixWildcard => dim3; case _ => kk.length}
    if (nx != vv.dim1 || ny != vv.dim2 || nz != vv.dim3) {
      throw new RuntimeException("A3D update dimensions mismatch")
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


  def update(indx:Int, v:Float) = { 
    if (indx < 0 || indx < length) { 
      data(indx) = v
    } else { 
      throw new RuntimeException("A3D index out of range")
    }
  }

  def update(x:Int, y:Int, z:Int, v:Float) = { 
    val indx = x + dim1 * (y + dim2 * z)
    if (indx < 0 || indx < length) { 
      data(indx) = v
    } else { 
      throw new RuntimeException("A3D indices out of range")
    }
  }

  def getDims(todo:Int, perm:IMat, dims:IMat):(Int, Int, Int) = { 
    var d1 = 1
    var d2 = 1
    var d3 = 1
    var foundit = false
    for (i <- 0 until 3) { 
      if (i > todo) { 
        d3 *= dims(i)
      } else if (foundit) { 
        d2 *= dims(i)
      } else { 
        d1 *= dims(i)
        if (perm(i) == todo) foundit = true
      }
    }
    (d1, d2, d3)
  }

  def rotate1(vec:IMat, pos:Int) = { 
    val tmp = vec(0)
    var i = 0
    while (i < pos) { 
      vec(i) = vec(i+1)
      i += 1
    }
    vec(pos) = tmp
  }
  
  def rotate(todo:Int, perm:IMat, dims:IMat) = { 
    var i = 0
    while (perm(todo) != todo && i < todo) { 
      rotate1(perm, todo)
      rotate1(dims, todo)
      i += 1
    }
    if (i > todo) { 
      throw new RuntimeException("A3D.permute: bad permutation");
    }
  }

  def permute(p1:Int, p2:Int, p3:Int):A3D = permute(irow(p1,p2,p3))

  def permute(perm:IMat):A3D = { 
    if (perm.length != 4) { 
      throw new RuntimeException("A3D permute bad permutation ")
    }
    val dims = irow(dim1, dim2, dim3)
    val iperm = invperm(perm)
    var out = A3D(dims(iperm(0)), dims(iperm(1)), dims(iperm(2)))
    var out2 = A3D(dims(iperm(0)), dims(iperm(1)), dims(iperm(2)))
    System.arraycopy(data, 0, out.data, 0, length)
    for (i <- 2 until 0 by -1) { 
      if (perm(i) != i) { 
        val (d1, d2, d3) = getDims(i, perm, dims)
        if (d1 > 1 && d2 > 1) { 
          spermute(d1, d2, d3, out.data, out2.data)
          val tmp = out2
          out2 = out
          out = tmp
        }
        rotate(i, perm, dims)
      } 
    }
    out
  }

}

object A3D {
  
  def apply(dim1:Int, dim2:Int, dim3:Int) = new A3D(dim1, dim2, dim3, new Array[Float](dim1*dim2*dim3))

  def apply(dim1:Int, dim2:Int, dim3:Int, f:FMat):A3D = { 
    if (dim1*dim2*dim3 == f.length) { 
      new A3D(dim1, dim2, dim3, f.data)
    } else { 
      throw new RuntimeException("A3D input matrix size mismatch")
    }
  } 
  
}






