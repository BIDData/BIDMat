package BIDMat
import edu.berkeley.bid.CBLAS._
import edu.berkeley.bid.LAPACK._
import edu.berkeley.bid.SPBLAS._
import edu.berkeley.bid.UTILS._


import java.util.Arrays
import java.util.Comparator
import scala.math._
import scala.reflect._

/*
 * TMat's are composed of vertical tiles, given as an Array of (Dense)Mats.  
 * This allows us to use a different multiplication strategy than for DenseMats.
 *
 * x(i) is the Mat.oneBased col index of the top left corner of tile i.
 * y(i) is the Mat.oneBased row index of the top left corent of tile i.
 *
 * y should be passed in sorted, and contiguous sections of x with the same
 * y value should also be sorted
 *
 * Does not yet check for overlaps!
 * 
 * if implemented..
 *     width(i) is the width of the i-th tile.
 *     height(i) is the height of the i-th tile.
 * doesn't seem strictly necessary as data has this information
 * 
 * tiles(i) is the i-th tile (as an FMat for now - need to change to Mat)
 * 
 */

class TMat (nr: Int, nc: Int, var x : Array[Int], var y : Array[Int], tiles : Array[FMat]) extends Mat(nr, nc) {
  require(x.length == y.length, "x.length must equal y.length")
  require(x.length == tiles.length, "x.length must equal tiles.length")
  

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
	  if ( (x(cindex) + tiles(cindex).ncols > c) && (c-x(cindex) < tiles(cindex).ncols) && (r-y(cindex) < tiles(cindex).nrows) && (y(cindex) <= r) ) {
           ret = tiles(cindex).apply(r-y(cindex),c-x(cindex)) 
           found = true
           return ret
	  }
           cindex += 1
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

  def full(mat:Mat): FMat  = {
    val out = FMat.newOrCheckFMat(nrows, ncols, mat)

    var i = 0
    var m = tiles(0)

    while (i < y.length) {
     m = tiles(i)
     val r:IMat = IMat(m.nrows,1,(y(i) to y(i)+m.nrows).toArray)
     val c:IMat = IMat(m.ncols,1,(x(i) to x(i)+m.ncols).toArray)

     out._update(r,c,m)
     i += 1
    }
    out
  }

/*
 * tMult takes advantage of the tiled structure of TFMat
 * 
 * We call tileMult repeatedly
 * tileMult(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:FMat, broff:Int, bcoff:Int, c:FMat, croff:Int, ccoff:Int)
 * 
 */ 
 
def tMult(a:FMat, outmat:Mat) : FMat =  {
      if (ncols == 1 && nrows == 1){ // left scalar multiply
  		val out = FMat.newOrCheckFMat(a.nrows, a.ncols, outmat, GUID, a.GUID, "dMult".##)
  		Mat.nflops += a.length
 
		var i = 0
  		val dvar = tiles(0).data(0)
  		while (i < a.length) {
  			out.data(i) = dvar * a.data(i)
  			i += 1
  		}	
  		out			  
      } else if (a.ncols == 1 && a.nrows == 1){ // right scalar multiply. should refactor to return a TFMat later
                val out = FMat.newOrCheckFMat(nrows, ncols, outmat, GUID, a.GUID, "dMult".##)
  		
  		var i = 0
  		val dvar = a.data(0)
                var m = tiles(0)

  		while (i < tiles.length) {
                        m = tiles(i)
                        val r:IMat = IMat(m.nrows,1,(y(i) to y(i)+m.nrows).toArray)
                        val c:IMat = IMat(m.ncols,1,(x(i) to x(i)+m.ncols).toArray)

  			out._update(r,c,m.fDMult(a,out(r,c)))
  			i += 1
  		}
 
  		out
      } else if (ncols == a.nrows) {
  		var out = FMat.newOrCheckFMat(nrows, a.ncols, outmat, GUID, a.GUID, "dMult".##)
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
                       // println("y(i): " + y(i))
                       // println("x(i): " + x(i))
                       out += m.tileMult(m.nrows, a.ncols, m.ncols, 0, 0, a, x(i), 0, tmp, y(i), 0)
                       // println("tmp:\n" + tmp)
                       // println("out:\n" + out)
  		  }
                 i+= 1  			 
		} 
               out
      }	else throw new RuntimeException("dimension mismatch")    

  }
  
}
