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


object IFunctions {
  import GMat.BinOp._
	
  /** Accumulate (row, col, value) tuples from inds \\ vals. nr and nc are row and column bounds */
  def accum(inds:IMat, vals:IMat, nr:Int, nc:Int):IMat = {
    (inds, vals) match {
//      case (ginds:GIMat, fvals:IMat) => GIFunctions.accum(ginds, GIMat(fvals), null, nr, nc);
//      case (finds:IMat, gvals:GIMat) => GIFunctions.accum(GIMat(finds), gvals, null, nr, nc);
      case _ => IMat(DenseMat.accum(inds, vals, nr, nc))
    }
  }
  
  /** Accumulate (row, col, value) tuples from inds \\ vals. nr is row and bounds, ncols = 1 */
  
  def accum(inds:IMat, vals:IMat, nr:Int):IMat = accum(inds, vals, nr, 1);
  
  /** Accumulate (row, value) tuples from inds \\ vals. Inds can be a vector or two-column matrix */
  def accum(inds:IMat, vals:IMat) = IMat(DenseMat.accum(inds, vals, 0, 0))
  
  /** Accumulate (row, col, value) tuples from inds \\ vals. nr and nc are row and column bounds */
  def accum(inds:IMat, v:Int, nr:Int, nc:Int) = {
    inds match {
//      case ginds:GIMat => GIFunctions.accum(ginds, v, null, nr, nc);
      case _ => IMat(DenseMat.accum(inds, IMat.ielem(v), nr, nc));
    }
  }   
  
  /** Accumulate (row, col, value) tuples from inds \\ vals. nr is row and bounds, ncols = 1 */
  def accum(inds:IMat, v:Int, nr:Int):IMat = accum(inds, v, nr, 1);
  
  /** Accumulate (row, value) tuples from inds \\ vals. Inds can be a vector or two-column matrix */
  def accum(inds:IMat, v:Int) = IMat(DenseMat.accum(inds, IMat.ielem(v), 0, 0));
  
	def min(a:IMat, b:IMat, out:Mat) = {
	  (a, b) match {
//	    case (aa:GIMat, bb:IMat) => GIFunctions.min(aa, GIMat(b), out);
//	    case (aa:IMat, bb:GIMat) => GIFunctions.min(GIMat(a), bb, out);
	    case _ => a.iiMatOpv(b, IMat.vecMinFun, op_min, out);
	  }
	}
	
	def max(a:IMat, b:IMat, out:Mat) = {
	  (a, b) match {
//	    case (aa:GIMat, bb:IMat) => GIFunctions.max(aa, GIMat(b), out);
//	    case (aa:IMat, bb:GIMat) => GIFunctions.max(GIMat(a), bb, out);
	    case _ => a.iiMatOpv(b, IMat.vecMaxFun, op_max, out);
	  }
	}
	
  def min(a:IMat, b:Int, out:Mat) = {
	  a match {
//	    case aa:GIMat=> GIFunctions.min(aa, GIMat.elem(b), out);
	    case _ => a.iiMatOpScalarv(b, IMat.vecMinFun, out);
	  }
	}
	
	def max(a:IMat, b:Int, out:Mat) = {
	  a match {
//	    case aa:GIMat=> GIFunctions.max(aa, GIMat.elem(b), out);
	    case _ => a.iiMatOpScalarv(b, IMat.vecMaxFun, out);
	  }
	}

	def maxi(a:IMat, n:Int, out:Mat) = {
	  a match {
//	    case aa:GIMat => GIFunctions.maxi(aa, n, out);
	    case _ => a.iiReduceOpv(n, IMat.idFun, IMat.vecMaxFun, out);
	  }
	}	
	
  def mini(a:IMat, n:Int, out:Mat) = {
	  a match {
//	    case aa:GIMat => GIFunctions.mini(aa, n, out);
	    case _ => a.iiReduceOpv(n, IMat.idFun, IMat.vecMinFun, out);
	  }
	}	
  
  def sum(a:IMat, n:Int, out:Mat) = {
	  a match {
//	    case aa:GIMat => GIFunctions.sum(aa, n, out);
	    case _ => a.iiReduceOpv(n, IMat.idFun, IMat.vecAddFun, out);
	  }
	}	
	
  def prod(a:IMat, n:Int, out:Mat) = {
	  a match {
//	    case aa:GIMat => GIFunctions.prod(aa, n, out);
	    case _ => a.iiReduceOpv(n, IMat.idFun, IMat.vecMulFun, out);
	  }
	}	

  
  def cumsum(a:IMat, n:Int, out:Mat) = {
		  a match {
//	    case aa:GIMat => GIFunctions.cumsum(aa, n, out);
	    case _ => a.iiReduceAll(n, IMat.idFun, IMat.sumFun, out);
		  }
  }
  
  def maxi2(a:IMat,d:Int):(IMat,IMat) = {
    a match {
//      case aa:GIMat => GIFunctions.maxi2(aa, null, null, d);
      case _ => {
    	  val (m,ii)=a.ggOpt2(d,IMat.gtPred); 
    	  (IMat(m), ii)
      }
    }
  }
  
  def mini2(a:IMat,d:Int):(IMat,IMat) = {
    a match {
//      case aa:GIMat => GIFunctions.mini2(aa, null, null, d);
      case _ => {
    	  val (m,ii)=a.ggOpt2(d,IMat.ltPred); 
    	  (IMat(m), ii)
      }
    }
  }

}