//-*-coding:utf-8-*-
package BIDMat
import MatFunctions._
import edu.berkeley.bid.CBLAS._
import java.util.Arrays
import java.util.concurrent.atomic._
import scala.concurrent.future
import scala.concurrent.ExecutionContext.Implicits.global


abstract class A4D(val dim1:Int, val dim2:Int, val dim3:Int, val dim4:Int) { 

  def length = dim1*dim2*dim3*dim4

  def size() = length;

  private var _GUID = Mat.myrand.nextLong
  
  def setGUID(v:Long):Unit = {_GUID = v}
  
  def GUID:Long = _GUID

  def apply(indx:Int):Float 

  def apply(x:Int, y:Int, z:Int, w:Int):Float 
  
  def apply(x:IMat, y:IMat, z:IMat, w:IMat):A4D

  def update(indx:Int, v:Float):A4D

  def update(x:Int, y:Int, z:Int, w:Int, v:Float):A4D
  
  def update(x:IMat, y:IMat, z:IMat, w:IMat, v:A4D):A4D

}

object A4D {

  
}






