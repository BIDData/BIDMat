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

/*
 * TMat's are composed of vertical tiles, given as an Array of (Dense)Mats.  
 * This allows us to use a different multiplication strategy than for DenseMats.
 *
 * x(i) is the Mat.oneBased col index of the top left corner of tile i.
 * y(i) is the Mat.oneBased row index of the top left corner of tile i.
 *
 * y should be passed in sorted, and contiguous x's with the same
 * y value should also be sorted
 *
 * Does not yet check for overlaps!
 * 
 * Everything is assumed single precision, so maybe TMat should go back to the original name:
 * TFMat
 * 
 * tiles(i) is the i-th tile as a generic Mat, but again, we assume single precision throughout
 * 
 */



class TMat 
      (nr: Int, nc: Int, var x : Array[Int], var y : Array[Int], tiles : Array[Mat])  extends Mat(nr, nc) {
  require(x.length == y.length, "x.length must equal y.length")
  require(x.length == tiles.length, "x.length must equal tiles.length")
  

  override def mytype = "TMat"

  def apply(r0: Int, c0: Int): Float = {
    val off = Mat.oneBased
    val r = r0 - off
    val c = c0 - off

    if (r < 0 || r >= nrows || c < 0 || c >= ncols) {
	throw new IndexOutOfBoundsException("("+(r+off)+","+(c+off)+") vs ("+nrows+","+ncols+")");
    } else {
	
       var rindex = math.max(Mat.ibinsearch(r,y,0,y.length),0)

       var tmp = rindex
       while (tmp > 0 && y(tmp) == y(rindex)) 
          tmp -= 1

       rindex = tmp    

       var cindex = rindex
       var ret = 0.0f
       var found = false
        
      /*
       * similar to saddleback search
       * review this later
       * 
       */

       while (rindex < y.length && y(rindex) <= r) {
         while (cindex < x.length && x(cindex) <= c) { 
          tiles(cindex) match {
	   case fMat : FMat =>           
                        if ( (x(cindex) + fMat.ncols > c) && (c-x(cindex) < fMat.ncols) && (r-y(cindex) < fMat.nrows) && (y(cindex) <= r) ) {
	 
                          found = true
                          return ret
	                }  else {  
                           cindex += 1 
                        }
           }                        
 	 }
        cindex = rindex	
        rindex += 1
       }

       if (found)
          ret
       else  
          throw new IndexOutOfBoundsException("not found");
    }
  }

  def full() : FMat = full(null)

  def full(mat:Mat) : FMat  = {
    val out = FMat.newOrCheckFMat(nrows, ncols, mat, GUID, "full".hashCode)
    out.clear

    var i = 0

    while (i < y.length) {
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
                          m.tileMult(m.nrows, a.ncols, m.ncols, 0, 0, a, x(i), 0, tmp, y(i), 0) 
                          out += tmp
  	            }
                  i+= 1			 
	         }
                out
      }	else throw new RuntimeException("dimension mismatch")    

  }
  
}
