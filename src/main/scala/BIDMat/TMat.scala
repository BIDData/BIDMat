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

import jcuda._
import jcuda.runtime._
import jcuda.runtime.JCuda._
import jcuda.runtime.cudaMemcpyKind._
import jcuda.jcublas._
import jcuda.jcublas.JCublas._
import jcuda.jcusparse._


/*
 * TMat's are composed of tiles, given as an Array of Mats and their companion subindices.
 * This allows us to use a novel multiplication strategy.
 *
 * x(i) is the Mat.oneBased col index of the top left corner of tile i.
 * y(i) is the Mat.oneBased row index of the top left corner of tile i.
 *
 * y should be passed in sorted, and contiguous x's with the same
 * y values should also be sorted (thus the order is lexicographic with priority on y)
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
        val x : Array[Int], 
        val y : Array[Int], 
        val tiles : Array[Mat] ) extends Mat(nr, nc) {

  require(x.length == y.length, "x.length must equal y.length")
  require(x.length == tiles.length, "x.length must equal tiles.length")
  

  override def mytype = "TMat"

  /*
   * Apply a (Mat, scalar) => Mat to a TMat tilewise.
   *
   */  

  def tMatOpScalarF(b : Float, f : (Mat, Float) => Mat, oldmat:TMat) : TMat = {
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
   *  tOp 
   *  
   *  Map over the tiles, applying a type specialized, ideally fast, Op
   *  This is more for the GPU types which can't take arbitrary lambdas
   *  but have ops denoted by number.
   *
   */ 

  def tOp(a : TMat, omat : Mat, op : Int) = {

  }

  /*
   * Apply a generate elementwise op on the tile level
   * probably slow.
   *
   * In most cases this can be applied only to CPU mats. 
   */
 
  def tOp(a : TMat, omat : Mat, op : (Mat,Mat) => Mat) : TMat = {
    var (newTiles,newXInds,newYInds) = if (omat.asInstanceOf[AnyRef] == null)
                                     {
                                       (Array.ofDim[Mat](tiles.length), 
                                        Array.ofDim[Int](tiles.length),
                                        Array.ofDim[Int](tiles.length))
                                     } else {
                                        (omat.asInstanceOf[TMat].tiles,
                                         omat.asInstanceOf[TMat].x,
                                         omat.asInstanceOf[TMat].y)
                                     }
   
    for (i <- 0 to (tiles.length-1)) {
       newTiles(i) = op(tiles(i), a.tiles(i))
       newXInds(i) = x(i)                                    // TODO add shape checking
       newYInds(i) = y(i)      
    }

    TMat(nr,nc,newXInds,newYInds,newTiles)
  }

  def tOp(a : Mat, omat : Mat, op : (Mat,Mat) => Mat) : Mat = {
    a match {
      case aa : TMat => tOp(aa,omat,op)
    }
  }

  /*
   * tMatOpF 
   * Apply a (scalar,scalar) => scalar elementwise within tiles to combine
   * two TMats.
   * 
   * Need to implement cacheing
   */
  
  def tMatOpF(aa : TMat, f : (Float, Float) => Float) : TMat = tMatOpF(aa,f,null)

  def tMatOpF(aa : TMat, f : (Float, Float) => Float, oldmat:TMat) : TMat = {
  
    var i = 0
    var out = TMat.newOrCheckTMat(nrows,ncols,x,y,tiles,oldmat)

    while (i < tiles.length) {
      Mat.nflops += tiles(i).length

      tiles(i) match {
        case fMat : FMat => {   
          out.tiles(i) = fMat.ffMatOp(aa.tiles(i), f, null) 
        }

        case sMat : SMat => {
          out.tiles(i) = sMat.ssMatOp((aa.tiles(i)).asInstanceOf[SMat], f, oldmat.tiles(i))
                                                                           // ^^ can cause null pointer excp    
        }

        case gMat : GMat => { 
          throw new RuntimeException("GMat not yet supported")
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

  override def set(f:Float) = {
    for (i <- 0 to (tiles.length-1)) {
      tiles(i).set(f)
    }
    this
  }
  
  var cacheFull:Mat = null

  def full() : Mat = {
      cacheFull = full(cacheFull)
      cacheFull
  }

  def full(mat:Mat) : Mat  = {

    val out = tiles(0) match {
      case t0 : FMat => FMat.newOrCheckFMat(nrows, ncols, mat, GUID, "full".hashCode)
      case t0 : GMat => GMat.newOrCheckGMat(nrows, ncols, mat, GUID, "full".hashCode)
    }
    out.clear

    var i = 0
      val a=SciFunctions.GPUmem._1

    while (i < tiles.length) {
     tiles(i) match {
       case fMat : FMat =>          
                val rowInds:IMat = IMat(fMat.nrows,1,(y(i) to y(i)+fMat.nrows).toArray)
                val colInds:IMat = IMat(fMat.ncols,1,(x(i) to x(i)+fMat.ncols).toArray)

                out(rowInds,colInds) = fMat
                i += 1
       case gMat : GMat =>          
                //val rowInds:IMat = IMat(gMat.nrows,1,(y(i) to y(i)+gMat.nrows).toArray)
                //val colInds:IMat = IMat(gMat.ncols,1,(x(i) to x(i)+gMat.ncols).toArray)
                //The previous approach will cause memory leak when casting into GIMat
                val rowInds:GIMat = GIMat(MatFunctions.irow(y(i) to y(i)+gMat.nrows-1))
                val colInds:GIMat = GIMat(MatFunctions.irow(x(i) to x(i)+gMat.ncols-1))
                out(rowInds,colInds) = gMat
                i += 1
       case _ => { throw new RuntimeException("no match in TMat.full for tiles " + i); }
        }
    }
      val aa=SciFunctions.GPUmem._1
      if (a!=aa)
        println(a,aa,out.GUID)

    out
  }

  def theta(x : Int) = {
    if (x < 0) 0 else 1 
  }

/* 
 * Colslice implementation. Needs cacheing review
 *
 */

override def colslice(left:Int, right:Int, omat:Mat) : TMat = {
    val ioff = Mat.ioneBased

    /* 
     * calculate new yInds, xInds
     */ 

    var i = 0
    var j = 0

    /*
     * Haven't worried about cacheing here because there shouldn't be too many indices
     */ 

    var newXinds = Array.ofDim[Int](tiles.length)
    var newYinds = Array.ofDim[Int](tiles.length)
    var newTiles = Array.ofDim[Mat](tiles.length)

    while (i < tiles.length) { 
      if (x(i)+tiles(i).ncols < left) {

     // these are tiles not present in result TMat
     // maybe feed these discarded matrices into cache later?
           
      } else if (x(i) < right) {
     // do slice

       newXinds(j) = (x(i)-left)*theta(x(i)-left)
       newYinds(j) = y(i)

       val localLeftSlice =  (left - x(i))*theta(left-x(i))
       val localRightSlice = (right - x(i))+(x(i)+tiles(i).ncols-right)*theta(right -x(i)- tiles(i).ncols)

       newTiles(j) = tiles(i).colslice(localLeftSlice,localRightSlice,null)  // could cache here too ?

       j += 1
      } else if (x(i) >= right) {

      // these tiles are also not present. case analysis just for clarity
      }

      i += 1
    }   

    TMat.newOrCheckTMat(nrows, right-left, newXinds, newYinds, newTiles,null) // <-- note cache omitted
  } 

/*
 * tMult is the method that the TMat container exists to implement. It's a fast matrix multiplication
 * filling in a subset of hypothetical dense matrix resulting in a Mat. 
 *
 * We call tileMult repeatedly
 * tileMult(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:FMat, broff:Int, bcoff:Int, c:FMat, croff:Int, ccoff:Int)
 * 
 * This is the TMat * FMat, TMat * SMat version
 * TMat = Mat * Mat is a static method below
 * 
 */ 
 
def tMult(a:Mat, outmat:Mat) : Mat =  {
   if (ncols == a.nrows) {
               a match { 
                  case aa : FMat => { 
                     var out = FMat.newOrCheckFMat(nrows, a.ncols, outmat, GUID, a.GUID, "tMult".##)
                     var i = 0
               
                     out.clear 

                     while (i < tiles.length) {
                          var m = tiles(i)
            	 	  Mat.nflops += 2L * m.length * a.ncols
                          if (!Mat.useMKL) {
                            out  // not sure
              		  } else {
                            m.tileMult(m.nrows, a.ncols, m.ncols, 0, 0, a, x(i), 0, out, y(i), 0)
        	          }
                        i+= 1			 
  	             }
                     out
                   }
                  case aa : SMat => { 
                     var out = FMat.newOrCheckFMat(nrows, a.ncols, outmat, GUID, a.GUID, "tMult".##)         
                     var i = 0

                     out.clear
               
                     while (i < tiles.length) {
                          var m = tiles(i)

            	 	  Mat.nflops += 2L * m.length * a.ncols
                          if (!Mat.useMKL) {
                            out  // not sure
              		  } else {
                            m.tileMult(m.nrows, a.ncols, m.ncols, 0, 0, a, x(i), 0, out, y(i), 0); 
        	          }
                        i+= 1			 
  	             }
                     out
                   }

                  case aa : GMat => { 
                     var out = GMat.newOrCheckGMat(nrows, a.ncols, outmat, GUID, a.GUID, "tMult".##)         
                     var i = 0

                     out.clear

                     while (i < tiles.length) {
                          var m = tiles(i)

            	 	  Mat.nflops += 2L * m.length * a.ncols
                          if (!Mat.useMKL) {
                            out  // not sure
              		  } else {
                            m.tileMult(m.nrows, a.ncols, m.ncols, 0, 0, a, x(i), 0, out, y(i), 0); 
        	          }
                        i+= 1			 
  	             }
                     out
                   }

                  case aa : GSMat => { 
                     var out = GMat.newOrCheckGMat(nrows, a.ncols, outmat, GUID, a.GUID, "tMult".##)    
                     var i = 0

                     out.clear

                     while (i < tiles.length) {
                          var m = tiles(i)

            	 	  Mat.nflops += 2L * m.length * a.ncols
                          if (!Mat.useMKL) {
                            out  // not sure
              		  } else {
                            m.tileMult(m.nrows, a.ncols, m.ncols, 0, 0, a, x(i), 0, out, y(i), 0); 
        	          }
                        i+= 1			 
  	             }
                     out
                   }
                 }
      }	else throw new RuntimeException("dimension mismatch")      
  }

def tMultT(a:Mat, outmat:Mat) : Mat =  {
   if (ncols == a.ncols) {
               a match { 
                  case aa : FMat => { 
                     var out = FMat.newOrCheckFMat(nrows, a.nrows, outmat, GUID, a.GUID, "tMultT".##)
                     var i = 0
               
                     out.clear 

                     while (i < tiles.length) {
                          var m = tiles(i)
            	 	  Mat.nflops += 2L * m.length * a.ncols
                          if (!Mat.useMKL) {
                            out  // not sure
              		  } else {
                            m.tileMultT(m.nrows, a.nrows, m.ncols, 0, 0, a, 0, x(i), out, y(i), 0)
        	          }
                        i+= 1			 
  	             }
                     out
                   }
                  case aa : SMat => { 
                     var out = FMat.newOrCheckFMat(nrows, a.nrows, outmat, GUID, a.GUID, "tMultT".##)         
                     var i = 0

                     out.clear
               
                     while (i < tiles.length) {
                          var m = tiles(i)

            	 	  Mat.nflops += 2L * m.length * a.ncols
                          if (!Mat.useMKL) {
                            out  // not sure
              		  } else {
                            m.tileMultT(m.nrows, a.nrows, m.ncols, 0, 0, a, 0, x(i), out, y(i), 0); 
        	          }
                        i+= 1			 
  	             }
                     out
                   }

                  case aa : GMat => { 
                     var out = GMat.newOrCheckGMat(nrows, a.nrows, outmat, GUID, a.GUID, "tMultT".##)         
                     var i = 0

                     out.clear

                     while (i < tiles.length) {
                          var m = tiles(i)

            	 	  Mat.nflops += 2L * m.length * a.ncols
                          if (!Mat.useMKL) {
                            out  // not sure
              		  } else {
                            m.tileMultT(m.nrows, a.ncols, m.ncols, 0, 0, a, 0, x(i), out, y(i), 0); 
        	          }
                        i+= 1			 
  	             }
                     out
                   }

                  case aa : GSMat => { 
                     var out = GMat.newOrCheckGMat(nrows, a.ncols, outmat, GUID, a.GUID, "tMultT".##)    
                     var i = 0

                     out.clear

                     while (i < tiles.length) {
                          var m = tiles(i)

            	 	  Mat.nflops += 2L * m.length * a.ncols
                          if (!Mat.useMKL) {
                            out  // not sure
              		  } else {
                            m.tileMultT(m.nrows, a.ncols, m.ncols, 0, 0, a, 0, x(i), out, y(i), 0); 
        	          }
                        i+= 1			 
  	             }
                     out
                   }
                 }
      }	else throw new RuntimeException("dimension mismatch")      
  }

  /*
   * The sum method for TMats
   *
   */

  def sum(n: Int, oldmat: Mat) : Mat = {
  // check if it's GPU or CPU, then iterate over tiles
  // calling sum on each, then aggregating
    val ioff = Mat.ioneBased

    val nn = if (n > 0) n else if (nrows == 1) 2 else 1
    
    var (out,tmp) = 
      tiles(0) match {
        case b:GMat => {
          ( GMat.newOrCheckGMat(if (nn==1) 1 else nrows, if (nn==1) ncols else 1, oldmat, GUID, n, "sum".##),
            GMat.newOrCheckGMat(if (nn==1) 1 else nrows, if (nn==1) ncols else 1, null) // <-- FIXME
          )
        }
        case b:GSMat => { 
          ( GMat.newOrCheckGMat(if (nn==1) 1 else nrows, if (nn==1) ncols else 1, oldmat, GUID, n, "sum".##),
            GMat.newOrCheckGMat(if (nn==1) 1 else nrows, if (nn==1) ncols else 1, null)
          )
        }
        case b:FMat => { 
         ( FMat.newOrCheckFMat(if (nn==1) 1 else nrows, if (nn==1) ncols else 1, oldmat, GUID, n, "sum".##),
           FMat.newOrCheckFMat(if (nn==1) 1 else nrows, if (nn==1) ncols else 1, null)
         )
        }
        case b:SMat => { 
         ( FMat.newOrCheckFMat(if (nn==1) 1 else nrows, if (nn==1) ncols else 1, oldmat, GUID, n, "sum".##),
           FMat.newOrCheckFMat(if (nn==1) 1 else nrows, if (nn==1) ncols else 1, null)
         )
        }
    }

    var i = 0
    while (i < tiles.length) {
       tmp <-- SciFunctions.sum(tiles(i),n)    // should pass tmp as "oldmat" here, I think, for caching
       // println("tiles(i): " + tiles(i))
       // println("tmp: " + tmp)

       val indexTuple = nn match {
         case 1 => (0,x(i))
         case _ => (y(i),0)
       }

       val offsetTuple = nn match {
         case 1 => (0,tiles(i).ncols-1)
         case _ => (tiles(i).nrows-1,0)
       }

       // println("index tuple: " + indexTuple)
       // println("offset tuple: " + offsetTuple)

       out(MatFunctions.irow(indexTuple._1 to (indexTuple._1 + offsetTuple._1)), 
           MatFunctions.irow(indexTuple._2 to (indexTuple._2 + offsetTuple._2))) += 
            tmp(MatFunctions.irow(0 to offsetTuple._1), MatFunctions.irow(0 to offsetTuple._2))
       i += 1
    }
    out
  }
  
  def op_tmat(t:TMat,o:TMat,op:Mop,op_num:Int) = {
      var out = TMat.newOrCheckTMat(nrows,ncols,x,y,tiles,o,GUID,t.GUID,op_num)
        var i=0
        while (i<tiles.length){
            op.op(tiles(i),t.tiles(i),out.tiles(i))
            i+=1
        }
        out      
  }
  
  def op_mat(t:Mat,o:TMat,op:Mop,op_num:Int) = {
      var out = TMat.newOrCheckTMat(nrows,ncols,x,y,tiles,o,GUID,t.GUID,op_num)
      var i=0
      while (i<tiles.length){
            op.op(tiles(i),t,out.tiles(i))
            i+=1
        }  
      out
    }
    
  def op(m:Mat,o:TMat,op:Mop,op_num:Int) = 
    m match {
        case t:TMat=>op_tmat(t,o,op,op_num)
        case _=>op_mat(m,o,op,op_num)
    }
    
  override def ~ (b: Mat) = b match { 
    case bb:TMat => new TPair(this,bb);
    case bb:Mat => new TTPair(this,bb);
  }
  

  override def zeros(nr: Int, nc: Int) = {
      if (nr == 1 && nc == 1) tiles(0).zeros(1,1)
      else
        TMat.zeros( nr,
                 nc,
                 x,
                 y,
                 TMat.cloneTiles(tiles))
  }

  override def ones(nr: Int, nc: Int) = {
     if (nr == 1 && nc == 1) tiles(0).ones(1,1)
      else
     TMat.ones ( nr,
                 nc,
                 x,
                 y,
                 TMat.cloneTiles(tiles))
  }
  
  override def toString:String = {
      "TMat(%d,%d,%d)" format (nr,nc,x.length)
  }
  
  override def clear = {
      var i = 0
      while (i < tiles.length) {
          tiles(i).clear
          i+=1
      }
      this
  }
  
  override def madd(b:Mat,c:Mat,at:Boolean,bt:Boolean) = {
      (b,c) match {
          case (gb:GMat,gc:GMat)=>
            madd(gb,gc,at,bt)
      }
  }
  
  def madd(b:GMat,c:GMat,at:Boolean,bt:Boolean) = {
      var i=0
      while (i<tiles.length) {
        bt match {
            case false=>{
                at match{
                    case true=>{
                    //cublasSgemm('t', 'n', ncols, a.ncols, nrows, 1.0f, data, nrows, a.data, a.nrows, 0f, out.data, out.nrows)
                        val a=tiles(i).asInstanceOf[GMat]
                        Mat.nflops += 2L * a.length * b.ncols
                        cublasSgemm('t', 'n', a.ncols, b.ncols, a.nrows, 1.0f, a.data, a.nrows, b.data.withByteOffset(Sizeof.FLOAT.toLong*(y(i))), b.nrows, 1.0f, c.data.withByteOffset(Sizeof.FLOAT.toLong*(x(i))), c.nrows)
                    }
                    case false=>{
                        val a=tiles(i).asInstanceOf[GMat]
                        Mat.nflops += 2L * a.length * b.ncols
                        cublasSgemm('n', 'n', a.nrows, b.ncols, a.ncols, 1.0f, a.data, a.nrows, b.data.withByteOffset(Sizeof.FLOAT.toLong*(x(i))), b.nrows, 1.0f, c.data.withByteOffset(Sizeof.FLOAT.toLong*(y(i))), c.nrows)
                    }
                }
            }
            case true=>{
                at match {
                    case true=>{
                        val a=tiles(i).asInstanceOf[GMat]
                        Mat.nflops += 2L * a.length * b.nrows
                        cublasSgemm('t', 't', a.ncols, b.nrows, a.nrows, 1.0f, a.data, a.nrows, b.data.withByteOffset(Sizeof.FLOAT.toLong*(y(i))*b.nrows), b.nrows, 1.0f, c.data.withByteOffset(Sizeof.FLOAT.toLong*(x(i))), c.nrows)
                    }
                    case false=>{
                        val a=tiles(i).asInstanceOf[GMat]
                        Mat.nflops += 2L * a.length * b.nrows
                        cublasSgemm('n', 't', a.nrows, b.nrows, a.ncols, 1.0f, a.data, a.nrows, b.data.withByteOffset(Sizeof.FLOAT.toLong*(x(i))*b.nrows), b.nrows, 1.0f, c.data.withByteOffset(Sizeof.FLOAT.toLong*(y(i))), c.nrows)
                    }           
                }
            }
        }
        i+=1  
      }
      c
  }
  
  
  def * (a : FMat) = this.tMult(a,null);
  def * (a : GMat) = this.tMult(a,null);
  def * (a : SMat) = this.tMult(a,null);
  def * (a : GSMat) = this.tMult(a,null);

  def *^ (a : FMat) = this.tMultT(a,null);
  def *^ (a : GMat) = this.tMultT(a,null);
  def *^ (a : SMat) = this.tMultT(a,null);
  def *^ (a : GSMat) = this.tMultT(a,null);

  override def * (a : Mat) = a match {
    case aa:FMat => this.tMult(a,null);
    case aa:GMat => this.tMult(a,null);
    case aa:SMat => this.tMult(a,null);
    case aa:GSMat => this.tMult(a,null); 
    case _ => throw new RuntimeException("no match in tMult");
  } 

  override def *^ (a : Mat) = a match {
    case aa:FMat => this.tMultT(a,null);
    case aa:GMat => this.tMultT(a,null);
    case aa:SMat => this.tMultT(a,null);
    case aa:GSMat => this.tMultT(a,null); 
    case _ => throw new RuntimeException("no match in tMultT");
  } 

  //def *@ (b : TMat) = tMatOpF(b, (x,y) => x*y, null)     // will need to specialize all of these for speed
  //def / (b : TMat) = tMatOpF(b, (x,y) => x/y, null)

  import TMat.BinOp._
  override def ^ (b : Float) = op_mat(GMat(b),null,Mop_Pow,op_pow)//tMatOpScalarF(b, (x,y) => x^y, null) // and these too
  override def *@ (b : Float) = op_mat(GMat(b),null,Mop_ETimes,op_mul)//tMatOpScalarF(b, (x,y) => x*y, null)
  
  override def *@ (m:Mat):Mat = op(m,null,Mop_ETimes,op_mul)/*{ 
        m match {
            case t:TMat=>
            case g:GMat=>op_mat(g,null,Mop_ETimes,op_mul)
        }
  }*/


  override def ^ (a : Mat) = op(a,null,Mop_Pow,op_pow)//tOp(a, null, (x:Mat,y:Mat) => x^y) 
  override def + (a : Mat) = op(a,null,Mop_Plus,op_add)//tOp(a, null, (x:Mat,y:Mat) => x^y) 
}

class TPair(val omat:Mat, val mat:TMat) extends Pair {
  import TMat.BinOp._
  override def * (a : Mat) = mat.tMult(a,null) // fix caching 
  override def ^ (a : Mat) = mat.tOp(a, omat, (x:Mat,y:Mat) => x^y) 
  override def *@ (a: Mat) = (omat) match {
      case (o:TMat)=>mat.op(a,o,Mop_ETimes,op_mul)
      //case (g:GMat,o:TMat)=>mat.op_mat(g,o)
  }
  override def + (a:Mat) = (omat) match {
    case (o:TMat)=>mat.op(a,o,Mop_Plus,op_add)    
  }
  
  override def / (a:Mat) = (omat) match {
    case (o:TMat)=>mat.op(a,o,Mop_EDiv,op_div)    
  }
  
  override def + (f:Float) = omat match {
      case o:TMat => mat.op_mat(GMat(f),o,Mop_Plus,op_add)
  }
}

class TTPair(val omat:Mat, val mat:Mat) extends Pair {
  override def * (a : Mat) = TMat.tMult(mat,a,omat.asInstanceOf[TMat]) 
  override def *^ (a : Mat) = TMat.tMultT(mat,a,omat.asInstanceOf[TMat]) 
}
 
object TMat {
    
  object BinOp {
  	val op_add=0
  	val op_sub=1
	  val op_mul=2
	  val op_div=3
	  val op_gt=4
	  val op_lt=5
	  val op_eq=6
	  val op_ge=7
	  val op_le=8
	  val op_ne=9
	  val op_max=10
	  val op_min=11
	  val op_atan2=12
	  val op_pow=13
  }  
    
  def cloneTiles(tiles:Array[Mat]) = {
      val out = tiles.clone()
      var i = 0
      while(i<tiles.length){
        out(i) = tiles(i).zeros(tiles(i).nrows,tiles(i).ncols)
        i+=1
      }
      out
  }

 def apply( nr:Int, 
            nc:Int, 
            xInds:Array[Int], 
            yInds:Array[Int], 
            data:Array[Mat] ) = {
                println("new",nr,nc)
    new TMat(nr, nc, xInds, yInds, data)
            }

  def TMatGPU ( tmat : TMat, omat : Mat ) : TMat = { 
   var i = 0
   val out = newOrCheckTMat(tmat.nrows, tmat.ncols, tmat.x, tmat.y, tmat.tiles, omat);
     
   while (i < tmat.tiles.length) {
    var tmp = tmat.tiles(i) match {
       case aa:FMat => GMat(aa)
       case aa:SMat => GSMat(aa)
       case aa:GMat => aa
       case aa:GSMat => aa
     }
    out.tiles(i) = tmp
    i +=1 
   }
   out
  }

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
      TMat(nr, nc, xInds, yInds, cloneTiles(data)) 
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

  def zeros ( nr: Int,
              nc: Int,
              xInds: Array[Int],
              yInds: Array[Int],
              data: Array[Mat] ) = {

    var i = 0
    while (i < data.length) {
      data(i).clear
      i += 1
    }

    new TMat( nr, nc, xInds, yInds, data)
  }

  def ones (  nr: Int,
              nc: Int,
              xInds: Array[Int],
              yInds: Array[Int],
              data: Array[Mat] ) = {

    var i = 0
    while (i < data.length) {
      data(i) = data(i).ones(data(i).nrows,data(i).ncols)
      i += 1
    }

    new TMat( nr, nc, xInds, yInds, data)
  }


  def tMult ( left: Mat, right: Mat, omat : TMat) : TMat = {
    var i = 0
    while (i < omat.tiles.length) {
      omat.tiles(i).clear
      left.tileMult( omat.tiles(i).nrows,   
                     omat.tiles(i).ncols,   
                     left.ncols,            
                     omat.y(i),             
                     0, 
                     right, 
                     0,
                     omat.x(i),             
                     omat.tiles(i),         
                     0,
                     0 )
      i += 1
    }
    omat
  }

  def tMultT ( left: Mat, right: Mat, omat : TMat) : TMat = {
    var i = 0
    while (i < omat.tiles.length) {
      omat.tiles(i).clear
      left.tileMultT( omat.tiles(i).nrows,  
                      omat.tiles(i).ncols,  
                      left.ncols,  
                      omat.y(i),   
                      0, 
                      right, 
                      omat.x(i),                         0,
                      omat.tiles(i),
                      0,
                      0 )
      i += 1
     }
     omat
   }
}