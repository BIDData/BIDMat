package BIDMat

import edu.berkeley.bid.CBLAS._
import edu.berkeley.bid.LAPACK._
import edu.berkeley.bid.SPBLAS._
import edu.berkeley.bid.UTILS._


import java.util.Arrays
import java.util.Comparator
import scala.math._
import scala.reflect._
import scala.math.Numeric._
import scala.collection.mutable.ArrayBuffer


/*
 * TMat's are composed of tiles, given as an Array of Mats and their companion subindices.
 * This allows us to use a novel multiplication strategy.
 *
 * x(i) is the Mat.oneBased col index of the top left corner of tile i.
 * y(i) is the Mat.oneBased row index of the top left corner of tile i.
 *
 * y should be passed in sorted, and contiguous x's with the same
 * y value should also be sorted
 *
 * We will check for this in the future.
 * Also, there is currently no overlap checking.
 * 
 * Everything is assumed single precision, so maybe TMat should go back to the original name:
 * TFMat
 * 
 * tiles(i) is the i-th tile as a generic Mat, but again, we assume single precision throughout
 * 
 */


class TMat 
      ( nr: Int, 
        nc: Int, 
        x : Array[Int], 
        y : Array[Int], 
        val tiles : Array[Mat] ) extends Mat(nr, nc) {

  require(x.length == y.length, "x.length must equal y.length")
  require(x.length == tiles.length, "x.length must equal tiles.length")
  

  override def mytype = "TMat"

  /*
   * Apply a (Mat, scalar) => Mat to a TMat tilewise.
   *
   */  

  def ggMatOpScalarF(b : Float, f : (Mat, Float) => Mat, oldmat:TMat) : TMat = {
    var i = 0
   
    var out = TMat.newOrCheckTMat(nrows,ncols,x,y,tiles,oldmat)

    while (i < tiles.length) {
      Mat.nflops += tiles(i).length
      out.tiles(i) = f(tiles(i), b)
      i += 1
    }

    out
  }

  /*
   * Apply a (scalar,scalar) => scalar elementwise within tiles to combine
   * two TMats.
   *
   */
  
/*
  def ggMatOpF(aa : TMat, f : (Float, Float) => Float, oldmat:TMat) : TMat = {
    var i = 0
   
    var out = TMat.newOrCheckTMat(nrows,ncols,x,y,tiles,oldmat)

    while (i < tiles.length) {
      Mat.nflops += tiles(i).length

      i += 1
    }

    out
  }
*/

  /*
   * An implementation of slicing
   * Surprisingly tricky
   *  
   */

  def gapply(rowinds:IMat, colinds:IMat): TMat = {
    var out:TMat = null

    /*
     * these could all probably be immutable Arrays
     * and the slicing implemented with a map, but
     * that seems a bit unnatural given the explicit copying
     *
     */

    var xInds = ArrayBuffer.empty[Int]
    var yInds = ArrayBuffer.empty[Int]
    val mats = ArrayBuffer.empty[Mat]

    val off = Mat.oneBased 
    var i = 0

    rowinds match {
    	case dummy:MatrixWildcard => {
    	   colinds match {
  	      case dummy2:MatrixWildcard => {
		   while (i < tiles.length) {
		     mats += tiles(i).copy
                     i += 1
                   }
                   
                    new TMat(nrows,ncols,x,y,mats.toArray)
		  }
	   }

	}
/*        case _ => {



        }
*/
    }
  }

  def full() : FMat = full(null)

  def full(mat:Mat) : FMat  = {
    val out = FMat.newOrCheckFMat(nrows, ncols, mat, GUID, "full".hashCode)
    out.clear

    var i = 0

    while (i < tiles.length) {
     // use contents to make this generic
     tiles(i) match {
       case fMat : FMat =>          
                val rowInds:IMat = IMat(fMat.nrows,1,(y(i) to y(i)+fMat.nrows).toArray)
                val colInds:IMat = IMat(fMat.ncols,1,(x(i) to x(i)+fMat.ncols).toArray)

                out._update(rowInds,colInds,fMat)
                i += 1
       case _ => { throw new RuntimeException("no match in TMat.full for tiles " + i); }
        }

    }
    out
  }

/*
 * tMult takes advantage of the tiled structure of TMat
 * 
 * We call tileMult repeatedly
 * tileMult(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:FMat, broff:Int, bcoff:Int, c:FMat, croff:Int, ccoff:Int)
 * 
 * Result is an FMat
 *
 */ 
 
def tMult(a:Mat, outmat:Mat) : FMat =  {
      if (ncols == 1 && nrows == 1){ // left scalar multiply
  		val out = FMat.newOrCheckFMat(a.nrows, a.ncols, outmat, GUID, a.GUID, "tMult".##)
  		Mat.nflops += a.length
 
		var i = 0

/*                tiles(0) match { 
                   case fMat : FMat => 
                       val dvar = fMat.data(0)
  		
                       while (i < a.length) {
                                out.data(i) = numeric.times(dvar,a.data(i))
         	         	i += 1
   	         	}	
                }

*/
  		out			  
      } else if (a.ncols == 1 && a.nrows == 1){ 
                val out = FMat.newOrCheckFMat(nrows, ncols, outmat, GUID, a.GUID, "tMult".##)
  		
  		var i = 0
/*  		val dvar = a.data(0)

  		while (i < tiles.length) {
                        tiles(i) match {
                           case fMat : FMat =>
		               val rowInds : IMat = IMat(fMat.nrows,1,(y(i) to y(i)+fMat.nrows).toArray)
                               val colInds : IMat = IMat(fMat.ncols,1,(x(i) to x(i)+fMat.ncols).toArray)
  			       out.update(rowInds,colInds,fMat.fDMult(a,out(rowInds,colInds)))

                	   }
        	      	i += 1
  		}
 */
  		out
      } else if (ncols == a.nrows) {
         	var out = FMat.newOrCheckFMat(nrows, a.ncols, outmat, GUID, a.GUID, "tMult".##)
                var tmp = FMat.newOrCheckFMat(nrows, a.ncols, null)

                var i = 0

                while (i < tiles.length) {
                  var m = tiles(i)
                  tmp.clear

		  Mat.nflops += 2L * m.length * a.ncols
                  if (!Mat.useMKL) {
                    // not sure
                    out 
		   } else {
                          m.tileMult(m.nrows, a.ncols, m.ncols, 0, 0, a, x(i), 0, tmp, y(i), 0); 
                          // println("m:   " + m);
                          // println("a:   " + a);
                          // println("x(i): " + x(i));
                          // println("y(i): " + y(i));
                          // println("tmp: " + tmp);
                          out += tmp
  	            }
                  i+= 1			 
	         }
                out
      }	else throw new RuntimeException("dimension mismatch")    

  }

}
 
object TMat {

 def apply( nr:Int, 
            nc:Int, 
            xInds:Array[Int], 
            yInds:Array[Int], 
            data:Array[Mat] ) = 
    new TMat(nr, nc, xInds, yInds, data)


 /*
  * This function is very unsafe at the moment.
  * In the future it should check that outmat is an instance of TMat, safely
  * then should check that the tile structure matches appropriately.
  *
  * Also we should define the recycle function for TMat's to grow them
  * to the right size if need be.
  *
  * The clone function is used because we expect to pass another TMat's tiles
  * in to this function. This is to deal with the fact that the memory layout 
  * of TMat is not known apriori ( e.g. without knowing the types and dimensions
  * of each tile ).
  *
  */

 def newOrCheckTMat( nr:Int, 
                     nc:Int, 
                     xInds: Array[Int], 
                     yInds: Array[Int],
                     data: Array[Mat],
                     outmat:Mat ) : TMat = {
    if (outmat.asInstanceOf[AnyRef] == null || (outmat.nrows == 0 && outmat.ncols == 0)) {
      TMat(nr, nc, xInds, yInds, data.clone()) 
    } else {
        outmat.asInstanceOf[TMat]
      }
  }
  
  def newOrCheckTMat( nr:Int, 
                      nc:Int, 
                      xInds: Array[Int], 
                      yInds: Array[Int],
                      data: Array[Mat],
                      outmat:Mat, matGuid:Long, opHash:Int ) : TMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckTMat(nr, nc, xInds, yInds, data, outmat)
    } else {
      val key = (matGuid, opHash)
      val res = Mat.cache2(key)
      if (res != null) {
      	newOrCheckTMat(nr, nc, xInds, yInds,data,res)
      } else {
        val omat = newOrCheckTMat(nr, nc, xInds, yInds, data, null)
        Mat.cache2put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckTMat( nr:Int, 
                      nc:Int, 
                      xInds: Array[Int],
                      yInds: Array[Int],
                      data: Array[Mat],
                      outmat:Mat, guid1:Long, guid2:Long, opHash:Int) : TMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckTMat(nr, nc, xInds, yInds, data, outmat)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      if (res != null) {
      	newOrCheckTMat(nr, nc, xInds, yInds, data, Mat.cache3(key))
      } else {
        val omat = newOrCheckTMat(nr, nc, xInds, yInds, data, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
  }
   
  def newOrCheckTMat( nr:Int, 
                      nc:Int, 
                      xInds: Array[Int],
                      yInds: Array[Int],
                      data: Array[Mat],
                      outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):TMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckTMat(nr, nc, xInds, yInds, data, outmat)
    } else {
      val key = (guid1, guid2, guid3, opHash)
      val res = Mat.cache4(key)
      if (res != null) {
      	newOrCheckTMat(nr, nc, xInds, yInds, data, Mat.cache4(key))
      } else {
        val omat = newOrCheckTMat(nr, nc, xInds, yInds, data, null)
        Mat.cache4put(key, omat)
        omat
      }
    }
  }
}