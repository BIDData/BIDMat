//-*-coding:utf-8-*-
package BIDMat
import MatFunctions._
import edu.berkeley.bid.CBLAS._
import java.util.Arrays
import java.util.concurrent.atomic._
import scala.concurrent.future
import scala.concurrent.ExecutionContext.Implicits.global


abstract class A3D(val dim1:Int, val dim2:Int, val dim3:Int) { 

  def length = dim1*dim2*dim3

  def size() = length;

  private var _GUID = Mat.myrand.nextLong
  
  def setGUID(v:Long):Unit = {_GUID = v}
  
  def GUID:Long = _GUID

  def apply(indx:Int):Float 

  def apply(x:Int, y:Int, z:Int):Float
  
  def apply(ii:IMat, jj:IMat, kk:IMat):A3D

  def update(indx:Int, v:Float):A3D

  def update(x:Int, y:Int, z:Int, v:Float):A3D
  
  def update(ii:IMat, jj:IMat, kk:IMat, v:Float):A3D
  
  def update(ii:IMat, jj:IMat, kk:IMat, vv:A3D):A3D 


}

object A3D {
  
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

  
}






