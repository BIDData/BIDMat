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
 * y value should also be sorted (thus the order is lexicographic with priority on y)
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
  
  def tMatOpF(aa : TMat, f : (Float, Float) => Float) : TMat = tMatOpF(aa,f,null)

  def tMatOpF(aa : TMat, f : (Float, Float) => Float, oldmat:TMat) : TMat = {
  
    var i = 0
    var out = TMat.newOrCheckTMat(nrows,ncols,x,y,tiles,oldmat)

    while (i < tiles.length) {
      Mat.nflops += tiles(i).length

      tiles(i) match {
        case fMat : FMat => {   
          // no cacheing for the moment
          out.tiles(i) = fMat.ffMatOp(aa.tiles(i), f, null) 
        }

        case sMat : SMat => {
          // there must be a better way than unsafe type coercions..
          // reflection is a nightmare

          out.tiles(i) = sMat.ssMatOp((aa.tiles(i)).asInstanceOf[SMat], f, oldmat.tiles(i))    
        }

        case gMat : GMat => { 
          throw new RuntimeException("GSMat not yet supported")
//          out.tiles(i) = gMat.gOp((aa.tiles(i)).asInstanceOf[GMat], oldmat.tiles(i), f)
        }

        case gSMat : GSMat => { 
          throw new RuntimeException("GSMat not yet supported")
        }
      }

      i += 1
    }

    out
  }


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

  def full() : Mat = full(null)

  def full(mat:Mat) : Mat  = {

    val out = tiles(0) match {
      case t0 : FMat => FMat.newOrCheckFMat(nrows, ncols, mat, GUID, "full".hashCode)
      case t0 : GMat => GMat.newOrCheckGMat(nrows, ncols, mat, GUID, "full".hashCode)
    }
    out.clear

    var i = 0

    while (i < tiles.length) {
     // use contents to make this generic
     tiles(i) match {
       case fMat : FMat =>          
                val rowInds:IMat = IMat(fMat.nrows,1,(y(i) to y(i)+fMat.nrows).toArray)
                val colInds:IMat = IMat(fMat.ncols,1,(x(i) to x(i)+fMat.ncols).toArray)

                out(rowInds,colInds) = fMat
                i += 1
       case gMat : GMat =>          
                val rowInds:IMat = IMat(gMat.nrows,1,(y(i) to y(i)+gMat.nrows).toArray)
                val colInds:IMat = IMat(gMat.ncols,1,(x(i) to x(i)+gMat.ncols).toArray)

                out(rowInds,colInds) = gMat
                i += 1
       case _ => { throw new RuntimeException("no match in TMat.full for tiles " + i); }
        }

    }
    out
  }


 /* colslice: 
  *
  * tiles completely to the left of the slice range are unaffected 
  *
  * tiles to the right of the slice range get their x coordinate reset only
  *
  * tiles in the slice range are sliced and their x coordinate is reset
  *
  */ 

/*
  override def colslice(left:Int, right:Int, omat:Mat) : TMat = {
    val ioff = Mat.ioneBased

    /* 
     * calculate new yInds, xInds
     */ 

    var i = 0
    var j = 0

    // should be logarithmically many indices, so no cacheing necessary

    var newXinds = new Array(tiles.length)
    var newYinds = new Array(tiles.length) 
    var used = new Array(tiles.length)

    while (i < tiles.length) { 
      if (xInds[i]+data[i].ncols < left) {
     // not present at all
    
      } else if (xInds[i] < left) {
     // do slice

       newXinds[j] = xInds[i]+data[i].ncols-left
       newYinds[j] = Yinds[i]

       used[i] = true
       j++
      } else if (xInds[i] >  right) {

      // also do nothing. case analysis just for clarity
      }

      i++
    }   

    /*
     * try to use cacheing, copy and return
     */ 

     var newData = omata.data.zipWithIndex.collect {
       case (x,i) if i % 3 == 0 => x
     }

     var out = TMat.newOrCheckTMat(right-left+1, ncols, newXinds, newYinds, omat.data)

      
  } 
*/


/*
 * tMult takes advantage of the tiled structure of TMat
 * 
 * We call tileMult repeatedly
 * tileMult(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:FMat, broff:Int, bcoff:Int, c:FMat, croff:Int, ccoff:Int)
 * 
 * Result is a Mat
 *
 * I couldn't figure out how to do tMult in place, so we need to use a tmp matrix.
 * This makes the caching a bit more complicated than in the other mults
 *
 */ 
 
def tMult(a:Mat, outmat:Mat, tmpmat: Mat) : Mat =  {
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
               a match { 
                  case aa : FMat => { 
                     
                     var (out,tmp) = 
                       	(FMat.newOrCheckFMat(nrows, a.ncols, outmat, GUID, a.GUID, "tMult".##), FMat.newOrCheckFMat(nrows, a.ncols, tmpmat, GUID, a.GUID, "tMult".##));
                                       
                     var i = 0
               
                     while (i < tiles.length) {
                          var m = tiles(i)
                          tmp.clear

            	 	  Mat.nflops += 2L * m.length * a.ncols
                          if (!Mat.useMKL) {
                            out  // not sure
              		  } else {
                            m.tileMult(m.nrows, a.ncols, m.ncols, 0, 0, a, x(i), 0, tmp, y(i), 0); 
                            out += tmp
        	          }
                        i+= 1			 
  	             }
                     out
                   }
                  case aa : SMat => { 
                     
                     var (out,tmp) = 
                       	(FMat.newOrCheckFMat(nrows, a.ncols, outmat, GUID, a.GUID, "tMult".##), FMat.newOrCheckFMat(nrows, a.ncols, tmpmat, GUID, a.GUID, "tMult".##));
                                       
                     var i = 0
               
                     while (i < tiles.length) {
                          var m = tiles(i)
                          tmp.clear

            	 	  Mat.nflops += 2L * m.length * a.ncols
                          if (!Mat.useMKL) {
                            out  // not sure
              		  } else {
                            m.tileMult(m.nrows, a.ncols, m.ncols, 0, 0, a, x(i), 0, tmp, y(i), 0); 
                            out += tmp
        	          }
                        i+= 1			 
  	             }
                     out
                   }

                  case aa : GMat => { 

                     var (out,tmp) = 
                       	(GMat.newOrCheckGMat(nrows, a.ncols, outmat, GUID, a.GUID, "tMult".##), GMat.newOrCheckGMat(nrows, a.ncols, tmpmat, GUID, a.GUID, "tMult".##));
                                       
                     var i = 0
                     out.clear

                     while (i < tiles.length) {
                          var m = tiles(i)
                          tmp.clear

            	 	  Mat.nflops += 2L * m.length * a.ncols
                          if (!Mat.useMKL) {
                            out  // not sure
              		  } else {
                            m.tileMult(m.nrows, a.ncols, m.ncols, 0, 0, a, x(i), 0, tmp, y(i), 0); 
                            out += tmp
        	          }
                        i+= 1			 
  	             }
                     out
                   }

                  case aa : GSMat => { 

                     var (out,tmp) = 
                       	(GMat.newOrCheckGMat(nrows, a.ncols, outmat, GUID, a.GUID, "tMult".##), GMat.newOrCheckGMat(nrows, a.ncols, tmpmat, GUID, a.GUID, "tMult".##));
                                       
                     var i = 0
                     out.clear

                     while (i < tiles.length) {
                          var m = tiles(i)
                          tmp.clear

            	 	  Mat.nflops += 2L * m.length * a.ncols
                          if (!Mat.useMKL) {
                            out  // not sure
              		  } else {
                            m.tileMult(m.nrows, a.ncols, m.ncols, 0, 0, a, x(i), 0, tmp, y(i), 0); 
                            out += tmp
        	          }
                        i+= 1			 
  	             }
                     out
                   }
                 }
      }	else throw new RuntimeException("dimension mismatch")    
  
  }

  def ~ (b: GMat) = new TGPair(this,b);
  def ~ (b: FMat) = new TFPair(this,b);
  def ~ (b: GSMat) = new TGSPair(this,b);
  def ~ (b: SMat) = new TSPair(this,b);

  override def ~ (b : Mat) : Pair = b match {
    case bb:GMat => new TGPair(this,bb);
    case bb:FMat => new TFPair(this,bb);
    case bb:GSMat => new TGSPair(this,bb);
    case bb:SMat => new TSPair(this,bb);
  }

  def * (a : FMat) = tMult(a,null,null);
  def * (a : GMat) = tMult(a,null,null);
  def * (a : SMat) = tMult(a,null,null);
  def * (a : GSMat) = tMult(a,null,null);


  def *@ (b : TMat) = tMatOpF(b, (x,y) => x*y, null)
  def / (b : TMat) = tMatOpF(b, (x,y) => x/y, null)
  
}

class TGPair(val left:TMat, val right:GMat) extends Pair {

}

class TFPair(val left:TMat, val right:FMat) extends Pair {

}

class TGSPair(val left:TMat, val right:GSMat) extends Pair {

}

class TSPair(val left:TMat, val right:SMat) extends Pair {

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