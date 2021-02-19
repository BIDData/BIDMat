package BIDMat

import edu.berkeley.bid.VML._
import edu.berkeley.bid.VSL
import edu.berkeley.bid.VSL._
import edu.berkeley.bid.CBLAS._
import edu.berkeley.bid.RAND;
import edu.berkeley.bid.RAND._;
import edu.berkeley.bid.SLATEC;
import java.util.Random._;
import SciState._
import org.apache.commons.math3.special._
import org.apache.commons.math3.distribution._
import org.apache.commons.math3.util.FastMath
import org.apache.commons.math3.random.RandomDataGenerator;


object LFunctions {
  import GMat.BinOp._
	
  /** Accumulate (row, col, value) tuples from inds \\ vals. nr and nc are row and column bounds */
  def accum(inds:IMat, vals:LMat, nr:Int, nc:Int):LMat = {
    (inds, vals) match {
//      case (ginds:GIMat, fvals:LMat) => GLFunctions.accum(ginds, GLMat(fvals), null, nr, nc);
//      case (finds:IMat, gvals:GLMat) => GLFunctions.accum(GIMat(finds), gvals, null, nr, nc);
      case _ => LMat(DenseMat.accum(inds, vals, nr, nc))
    }
  }
  
  /** Accumulate (row, col, value) tuples from inds \\ vals. nr is row and bounds, ncols = 1 */
  
  def accum(inds:IMat, vals:LMat, nr:Int):LMat = accum(inds, vals, nr, 1);
  
  /** Accumulate (row, value) tuples from inds \\ vals. Inds can be a vector or two-column matrix */
  def accum(inds:IMat, vals:LMat) = LMat(DenseMat.accum(inds, vals, 0, 0))
  
  /** Accumulate (row, col, value) tuples from inds \\ vals. nr and nc are row and column bounds */
  def accum(inds:IMat, v:Long, nr:Int, nc:Int) = {
    inds match {
//      case ginds:GIMat => GLFunctions.accum(ginds, v, null, nr, nc);
      case _ => LMat(DenseMat.accum(inds, LMat.lelem(v), nr, nc));
    }
  }   
  
  /** Accumulate (row, col, value) tuples from inds \\ vals. nr is row and bounds, ncols = 1 */
  def accum(inds:IMat, v:Long, nr:Int):LMat = accum(inds, v, nr, 1);
  
  /** Accumulate (row, value) tuples from inds \\ vals. Inds can be a vector or two-column matrix */
  def accum(inds:IMat, v:Long) = LMat(DenseMat.accum(inds, LMat.lelem(v), 0, 0));
  
	def min(a:LMat, b:LMat, out:Mat) = {
	  (a, b) match {
//	    case (aa:GLMat, bb:LMat) => GLFunctions.min(aa, GLMat(b), out);
//	    case (aa:LMat, bb:GLMat) => GLFunctions.min(GLMat(a), bb, out);
	    case _ => a.iiMatOpv(b, LMat.vecMinFun, op_min, out);
	  }
	}
	
	def max(a:LMat, b:LMat, out:Mat) = {
	  (a, b) match {
//	    case (aa:GLMat, bb:LMat) => GLFunctions.max(aa, GLMat(b), out);
//	    case (aa:LMat, bb:GLMat) => GLFunctions.max(GLMat(a), bb, out);
	    case _ => a.iiMatOpv(b, LMat.vecMaxFun, op_max, out);
	  }
	}
	
  def min(a:LMat, b:Int, out:Mat) = {
	  a match {
//	    case aa:GLMat=> GLFunctions.min(aa, GLMat.elem(b), out);
	    case _ => a.iiMatOpScalarv(b, LMat.vecMinFun, out);
	  }
	}
	
	def max(a:LMat, b:Int, out:Mat) = {
	  a match {
//	    case aa:GLMat=> GLFunctions.max(aa, GLMat.elem(b), out);
	    case _ => a.iiMatOpScalarv(b, LMat.vecMaxFun, out);
	  }
	}

	def maxi(a:LMat, n:Int, out:Mat) = {
	  a match {
//	    case aa:GLMat => GLFunctions.maxi(aa, n, out);
	    case _ => a.iiReduceOpv(n, LMat.idFun, LMat.vecMaxFun, out);
	  }
	}	
	
  def mini(a:LMat, n:Int, out:Mat) = {
	  a match {
//	    case aa:GLMat => GLFunctions.mini(aa, n, out);
	    case _ => a.iiReduceOpv(n, LMat.idFun, LMat.vecMinFun, out);
	  }
	}	
  
  def sum(a:LMat, n:Int, out:Mat) = {
	  a match {
//	    case aa:GLMat => GLFunctions.sum(aa, n, out);
	    case _ => a.iiReduceOpv(n, LMat.idFun, LMat.vecAddFun, out);
	  }
	}	
	
  def prod(a:LMat, n:Int, out:Mat) = {
	  a match {
//	    case aa:GLMat => GLFunctions.prod(aa, n, out);
	    case _ => a.iiReduceOpv(n, LMat.idFun, LMat.vecMulFun, out);
	  }
	}	

  
  def cumsum(a:LMat, n:Int, out:Mat) = {
		  a match {
//	    case aa:GLMat => GLFunctions.cumsum(aa, n, out);
	    case _ => a.iiReduceAll(n, LMat.idFun, LMat.sumFun, out);
		  }
  }
  
  def maxi2(a:LMat,d:Int):(LMat,IMat) = {
    a match {
//      case aa:GLMat => GLFunctions.maxi2(aa, null, null, d);
      case _ => {
    	  val (m,ii)=a.ggOpt2(d,LMat.gtPred); 
    	  (LMat(m), ii)
      }
    }
  }
  
  def mini2(a:LMat,d:Int):(LMat,IMat) = {
    a match {
//      case aa:GLMat => GLFunctions.mini2(aa, null, null, d);
      case _ => {
    	  val (m,ii)=a.ggOpt2(d,LMat.ltPred); 
    	  (LMat(m), ii)
      }
    }
  }

}