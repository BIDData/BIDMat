package BIDMat
import scala.reflect._

object MatHDF5 {
  def addOne(ii:Array[Int]) = {
  	if (Mat.ioneBased == 1) {
  		var i = 0;
  		while (i < ii.length) {
  			ii(i) += 1;
  			i += 1;
  		}
  	}
  }

  def subOne(ii:Array[Int]) = {
  	if (Mat.ioneBased == 1) {
  		var i = 0;
  		while (i < ii.length) {
  			ii(i) = ii(i) - 1;
  			i += 1;
  		}
  	}
  }
}
