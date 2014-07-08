//-*-coding:utf-8-*-
package BIDMat
import MatFunctions._
import edu.berkeley.bid.CBLAS._
import java.util.Arrays
import java.util.concurrent.atomic._
import scala.concurrent.future
import scala.concurrent.ExecutionContext.Implicits.global


case class Tens4(val dim1:Int, val dim2:Int, val dim3:Int, val dim4:Int, val data:Array[Float]) { 

  def length = dim1*dim2*dim3*dim4

  def size() = length;

  private var _GUID = Mat.myrand.nextLong
  
  def setGUID(v:Long):Unit = {_GUID = v}
  
  def GUID:Long = _GUID

  def apply(indx:Int):Float = { 
    if (indx < length) { 
      data(indx)
    } else { 
      throw new RuntimeException("Tens4 index out of range")
    }
  }

  def apply(x:Int, y:Int, z:Int, w:Int):Float = { 
    val indx = x + dim1 * (y + dim2 * (z + dim3 * w))
    if (indx < length) { 
      data(indx)
    } else { 
      throw new RuntimeException("Tens4 indices out of range")
    }
  }


  def update(indx:Int, v:Float) = { 
    if (indx < 0 || indx < length) { 
      data(indx) = v
    } else { 
      throw new RuntimeException("Tens4 index out of range")
    }
  }

  def update(x:Int, y:Int, z:Int, w:Int, v:Float) = { 
    val indx = x + dim1 * (y + dim2 * (z + dim3 * w))
    if (indx < 0 || indx < length) { 
      data(indx) = v
    } else { 
      throw new RuntimeException("Tens4 indices out of range")
    }
  }

  def getDims(todo:Int, perm:IMat, dims:IMat):(Int, Int, Int) = { 
    var d1 = 1
    var d2 = 1
    var d3 = 1
    var foundit = false
    for (i <- 0 until 4) { 
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
      throw new RuntimeException("Tens4.permute: bad permutation");
    }
  }

  def permute(p1:Int, p2:Int, p3:Int, p4:Int):Tens4 = permute(irow(p1,p2,p3,p4))

  def permute(perm:IMat):Tens4 = { 
    if (perm.length != 4) { 
      throw new RuntimeException("Tens4 permute bad permutation ")
    }
    val dims = irow(dim1, dim2, dim3, dim4)
    val iperm = invperm(perm)
    var out = Tens4(dims(iperm(0)), dims(iperm(1)), dims(iperm(2)), dims(iperm(3)))
    var out2 = Tens4(dims(iperm(0)), dims(iperm(1)), dims(iperm(2)), dims(iperm(3)))
    System.arraycopy(data, 0, out.data, 0, length)
    for (i <- 3 until 0 by -1) { 
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

object Tens4 {
  
  def apply(dim1:Int, dim2:Int, dim3:Int, dim4:Int) = new Tens4(dim1, dim2, dim3, dim4, new Array[Float](dim1*dim2*dim3*dim4))

  def apply(dim1:Int, dim2:Int, dim3:Int, dim4:Int, f:FMat):Tens4 = { 
    if (dim1*dim2*dim3*dim4 == f.length) { 
      new Tens4(dim1, dim2, dim3, dim4, f.data)
    } else { 
      throw new RuntimeException("Tens4 input matrix size mismatch")
    }
  } 
  
}






