
package BIDMat

import jcuda._
import jcuda.runtime._
import jcuda.runtime.JCuda._
import jcuda.runtime.cudaMemcpyKind._
import jcuda.jcublas._
import jcuda.jcublas.JCublas2._
import jcuda.jcusparse._
import scala.concurrent.Future
import scala.concurrent.ExecutionContext.Implicits.global
import scala.util.hashing.MurmurHash3
import edu.berkeley.bid.MurmurHash3.MurmurHash3_x64_64
import edu.berkeley.bid.CUMAT
import java.io.ObjectOutputStream
import java.io.ObjectInputStream

@SerialVersionUID(100L)
class GMat(dims0:Array[Int], @transient var pdata:Pointer, val realsize:Long) extends FMat(dims0, null) {
  import GMat.BinOp._
  
  /** 2D Constructor */
  def this(nr:Int, nc:Int, pdata:Pointer, realsize:Long) = this(Array(nr, nc), pdata, realsize);
    
  override def mytype = "GMat";
  
  override def nnz = length;
  
  override def t = {
    val out = GMat.newOrCheckGMat(ncols, nrows, null, GUID, "t".##)
    CUMAT.transpose(this.pdata, nrows, out.pdata, ncols, nrows, ncols)
    out
  }
  
  /** hold indices in GPU mem */
  val ginds = new Array[GIMat](ndims);

  /** Access values and contents */
  override def dv:Double =
    if (nrows > 1 || ncols > 1) {
      throw new RuntimeException("Matrix should be 1x1 to extract value")
    } else {
      toFMat(null).data(0)
    }
  
  override def fv:Float =
    if (nrows > 1 || ncols > 1) {
      throw new RuntimeException("Matrix should be 1x1 to extract value")
    } else {
      toFMat(null).data(0)
    }
  
  override def contents() = {
    val out = new GMat(length, 1, pdata, realsize);
    out.setGUID(MurmurHash3.mix(MurmurHash3.mix(length, 1), (GUID*7897889).toInt));
    out
  }
  
  override def view(nr:Int, nc:Int):GMat = {
    if (1L * nr * nc > realsize) {
      throw new RuntimeException("view dimensions too large")
    }
    if (nr == nrows && nc == ncols) {
      this
    } else {
    	val out = new GMat(nr, nc, pdata, realsize);
    	out.setGUID(MurmurHash3.mix(MurmurHash3.mix(nr, nc), (GUID*3145341).toInt));
    	out
    }
  } 
    
  var saveMe:FMat = null;
  
  var printVal:Array[Float] = null;
  
  override def printOne(i:Int):String = {
    if (printVal.asInstanceOf[AnyRef] == null) printVal = new Array[Float](1);
    GMat.GPUtoCPUarraycopy(pdata, i,  printVal, 0, 1, "printOne");
    val v = printVal(0);
    if (v % 1 == 0 && math.abs(v) < 1e10) {
      "%d" format v.intValue
    } else {
      "%.5g" format v
    }
  }
  
  private def writeObject(out:ObjectOutputStream):Unit = {
    saveMe = FMat(this);
  	out.defaultWriteObject();
  }
  
  private def readObject(in:ObjectInputStream):Unit = {
    in.defaultReadObject();
    val gpu = SciFunctions.getGPU;
    GFunctions.setGPU(myGPU);
    pdata = GMat(saveMe).pdata;
    GFunctions.setGPU(gpu);
    saveMe = null;
  }
  
  def getIndxWrapper(i:Int):GIMat = {
    if (ginds(i).asInstanceOf[AnyRef] == null) {
      ginds(i) = GIMat(1,1)
    }
    ginds(i);
  }
  
  def getIndexValue(i:Int, v:Int):GIMat = {
    getIndxWrapper(i);
    ginds(i).set(v);
    ginds(i);
  }
  
  def getIndexMat(i:Int, v:IMat):GIMat = {
    if (v.length == 0 || v.length > 1) {
      GIMat(v);
    } else {
    	getIndxWrapper(i).set(v.v);
      ginds(i);
    }
  }
   
  def safePointer(ind:GIMat):Pointer = {
    if (ind.asInstanceOf[AnyRef] == null) {
      GMat.nullPointer;
    } else {
      ind.pdata;
    }
  }
  
  /** 1D access */
  
  override def apply(ind:Int):Float = {
  	val tmp = new Array[Float](1);
    GMat.GPUtoCPUarraycopy(pdata, ind, tmp, 0, 1, "GMat apply");
    tmp(0)
  }
  
  /** 2D access */
  
  override def apply(i:Int, j:Int):Float = {
    val tmp = new Array[Float](1);
    GMat.GPUtoCPUarraycopy(pdata, i + nrows * j, tmp, 0, 1, "GMat apply");
    tmp(0)
  }
  
  /** ND access */
  
  override def applyv(inds:Array[Int]):Float = {
    val indx = ND.linearize(inds, dims.data);
    val tmp = new Array[Float](1);
    GMat.GPUtoCPUarraycopy(pdata, indx, tmp, 0, 1, "GMat apply");
    tmp(0);
  }
  
    /** apply to an index IMat, and mirror its structure in the result */
  /* should be implemented in subclasses */
  
  override def apply(inds:IMat):GMat = {
    	inds match {
    	case aa:MatrixWildcard => {
    		val out = GMat.newOrCheckGMat(length, 1, null, GUID, inds.GUID, "apply(?)".##);
    		GMat.GPUtoGPUarraycopy(pdata, 0,  out.pdata, 0, length, "GMat IMat apply" );
    		out
    	}
    	case _ => {
    		val newinds = getIndexMat(0, inds);
    		val out = GMat.newOrCheckGMat(inds.dims, null, GUID, inds.GUID, "apply IMat".##);
    		val err = CUMAT.copyFromInds(pdata, out.pdata, safePointer(newinds), inds.length);
        if (err != 0) throw new RuntimeException("GMat apply(I) error" + cudaGetErrorString(err));
        out;
    	}
    	}
    }
 
//  override def apply(i1:IMat):GMat = applyi(Array(i1), null);
  override def apply(i1:IMat, i2:IMat):GMat = applyi(Array(i1, i2), null);
  override def apply(i1:IMat, i2:Int):GMat = applyi(Array(i1, IMat.ielem(i2)), null);
  override def apply(i1:Int, i2:IMat):GMat = applyi(Array(IMat.ielem(i1), i2), null);
  
  override def applyi(inds:Array[IMat]):GMat = applyi(inds, null);
 
  def applyi(inds:Array[IMat], omat:Mat):GMat = { 
    if (inds.length > 2 && inds.length != _dims.length) throw new RuntimeException("GMat applyi dims must match");
    val mydims = if (inds.length == 1) Array(length) else  if (inds.length == 2) Array(nrows, ncols) else _dims;
    val newdims = new Array[Int](inds.length)
    val newinds = new Array[GIMat](inds.length)
    for (i <- 0 until inds.length) {
      inds(i) match {
        case aa:MatrixWildcard => {
          newdims(i) = mydims(i); 
        }
        case _ => {
          newdims(i) = inds(i).length;
          newinds(i) = getIndexMat(i, inds(i));
        }
      }
    }
    val out = GMat.newOrCheckGMat(newdims, omat, GUID, ND.hashGUIDs(inds), "apply".##);
    inds.length match {
    case 1 => {
        val err = CUMAT.copyFromInds(pdata, out.pdata, safePointer(newinds(0)), newdims(0));
        if (err != 0) throw new RuntimeException("GMat apply(I) error: " + cudaGetErrorString(err));
      }
      case 2 => {
        val err = CUMAT.copyFromInds2D(pdata, dims(0), out.pdata, newdims(0), safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1));
        if (err != 0) throw new RuntimeException("GMat apply(I, J) error: " + cudaGetErrorString(err));
      }
      case 3 => {
        val err = CUMAT.copyFromInds3D(pdata, dims(0), dims(1), out.pdata, newdims(0), newdims(1), 
            safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1), safePointer(newinds(2)), newdims(2));
        if (err != 0) throw new RuntimeException("GMat apply(I, J, K) error: " + cudaGetErrorString(err));
      }
      case 4 => {
        val err = CUMAT.copyFromInds4D(pdata, dims(0), dims(1), dims(2), out.pdata, newdims(0), newdims(1), newdims(2),
            safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1), safePointer(newinds(2)), newdims(2), safePointer(newinds(3)), newdims(3));
        if (err != 0) throw new RuntimeException("GMat apply(I, J, K, L) error: " + cudaGetErrorString(err));
      }   
      case _ => throw new RuntimeException("GMat slice access with more than 4 indices not supported");
    }
    out;
  }
  
  /** 1D update */
  
  override def update(i:Int, v:Float):GMat = {
    val tmp = new Array[Float](1);
    tmp(0) = v;
    GMat.CPUtoGPUarraycopy(tmp, 0, pdata, i, 1, "GMat update");
    this
  }
  
  /** 2D update */
  
  override def update(i:Int, j:Int, v:Float):GMat = {
    val tmp = new Array[Float](1);
    tmp(0) = v;
    GMat.CPUtoGPUarraycopy(tmp, 0, pdata, i + j * nrows, 1, "GMat update");
    this
  }
  
  /** ND update */
  
  override def updatev(inds:Array[Int], v:Float):GMat = {
    val indx = ND.linearize(inds, _dims); 
    val tmp = Array[Float](v);
    GMat.CPUtoGPUarraycopy(tmp, 0, pdata, indx, 1, "GMat update");
    this
  } 
 
  /** ND sliced updates */
  
  override def update(inds:IMat, vv:FMat):GMat = updatei(inds, GMat(vv));
  
  override def update(iv:IMat, jv:IMat, b:FMat):FMat = updatei(Array(iv, jv), GMat(b));
  override def update(iv:IMat, j:Int, b:FMat):FMat = updatei(Array(iv, IMat.ielem(j)), GMat(b));
  override def update(i:Int, jv:IMat, b:FMat):FMat = updatei(Array(IMat.ielem(i), jv), GMat(b));
  
  override def update(inds:IMat, vv:Mat):FMat = updatei(inds, GMat(vv));
  
  override def update(iv:IMat, jv:IMat, b:Mat):FMat = updatei(Array(iv, jv), GMat(b));
  override def update(iv:IMat, j:Int, b:Mat):FMat = updatei(Array(iv, IMat.ielem(j)), GMat(b));
  override def update(i:Int, jv:IMat, b:Mat):FMat = updatei(Array(IMat.ielem(i), jv), GMat(b));
 
//  override def update(i1:IMat, vv:Float):FMat = updatei(Array(i1), vv);
  override def update(i1:IMat, i2:IMat, vv:Float):FMat = updatei(Array(i1, i2), vv);
 
  override def updatei(inds:Array[IMat], vv:FMat):GMat = updatei(inds, GMat(vv));
  
  def updatei(inds:IMat, vv:GMat):GMat = {
		val newinds = getIndexMat(0, inds);
		val err = inds match {
		case aa:MatrixWildcard => {
			if (vv.length != length) throw new RuntimeException("GMat column update length mismatch")
			CUMAT.copyToInds(vv.pdata, pdata, safePointer(newinds), length);
		}
		case _ => {
			if (inds.length != vv.length) throw new RuntimeException("GMat column update length mismatch")
			CUMAT.copyToInds(vv.pdata, pdata, safePointer(newinds), inds.length);
		}
		}
		if (err != 0) throw new RuntimeException("GMat update (I)=v error " + cudaGetErrorString(err));
		this;
  }

  def updatei(inds:Array[IMat], vv:GMat):GMat = {
		 if (inds.length > 2 && inds.length != _dims.length) throw new RuntimeException("GMat update dims must match")
		 val mydims = if (inds.length == 2) Array(nrows, ncols) else _dims;
		 val newdims = new Array[Int](inds.length)
		 val newinds = new Array[GIMat](inds.length)
		 var j = 0
		 for (i <- 0 until inds.length) {
			 inds(i) match {
			 case aa:MatrixWildcard => {
				 newdims(i) = mydims(i); 
			 }
			 case _ => {
				 newdims(i) = inds(i).length;
				 newinds(i) = getIndexMat(i, inds(i));
			 }
			 }
		 }
		 ND.checkDims("GMat update:", ND.trimDims(newdims), ND.trimDims(vv._dims));
		 inds.length match {
		 case 2 => {
			 val err = CUMAT.copyToInds2D(vv.pdata, vv.dims(0), pdata, nrows, 
					 safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1));
			 if (err != 0) throw new RuntimeException("GMat update (I, J)=V error " + cudaGetErrorString(err));
		 }
		 case 3 => {
			 val err = CUMAT.copyToInds3D(vv.pdata, vv.dims(0), vv.dims(1), pdata, dims(0), dims(1), 
					 safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1), safePointer(newinds(2)), newdims(2));
			 if (err != 0) throw new RuntimeException("GMat update (I, J, K)=V error " + cudaGetErrorString(err));
		 }
		 case 4 => {
			 val err = CUMAT.copyToInds4D(vv.pdata, vv.dims(0), vv.dims(1), vv.dims(2), pdata, dims(0), dims(1), dims(2),
					 safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1), safePointer(newinds(2)), newdims(2), safePointer(newinds(3)), newdims(3));
			 if (err != 0) throw new RuntimeException("GMat udpate (I, J, K, L)=V error " + cudaGetErrorString(err));
		 }
		 case _ => throw new RuntimeException("GMat slice access with more than 4 indices not supported");
		 }
		 this
 }
  
  override def update(inds:IMat, vv:Float):GMat = {
    val newinds = getIndexMat(0, inds);
    val err = inds match {
        case aa:MatrixWildcard => {
          CUMAT.fillToInds(vv, pdata, safePointer(newinds), length);
        }
        case _ => {
        	CUMAT.fillToInds(vv, pdata, safePointer(newinds), inds.length);
        }
    }
    if (err != 0) throw new RuntimeException("GMat update (I)=v error " + cudaGetErrorString(err));
    this;
 }
 
  override def updatei(inds:Array[IMat], vv:Float):GMat = {
		if (inds.length > 2 && inds.length != _dims.length) throw new RuntimeException("GMat update dims must match");
		val mydims = if (inds.length == 2) Array(nrows, ncols) else _dims;
    val newdims = new Array[Int](inds.length);
    val newinds = new Array[GIMat](inds.length);
    for (i <- 0 until inds.length) {
      inds(i) match {
        case aa:MatrixWildcard => {
          newdims(i) = mydims(i); 
        }
        case _ => {
          newdims(i) = inds(i).length;
          newinds(i) = getIndexMat(i, inds(i));
        }
      }
    }
    inds.length match {
    case 2 => {
      val err = CUMAT.fillToInds2D(vv,  pdata, nrows, 
          safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1));
      if (err != 0) throw new RuntimeException("GMat update (I, J)=v error " + cudaGetErrorString(err));
    }
    case 3 => {
      val err = CUMAT.fillToInds3D(vv, pdata, dims(0), dims(1), 
          safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1), safePointer(newinds(2)), newdims(2));
      if (err != 0) throw new RuntimeException("GMat update (I, J, K)=v error " + cudaGetErrorString(err));
    }
    case 4 => {
      val err = CUMAT.fillToInds4D(vv, pdata, dims(0), dims(1), dims(2),
          safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1), safePointer(newinds(2)), newdims(2), safePointer(newinds(3)), newdims(3));
      if (err != 0) throw new RuntimeException("GMat udpate (I, J, K, L)=v error " + cudaGetErrorString(err));
    }
    case _ => throw new RuntimeException("GMat slice update with more than 4 indices not supported");
    }
    this
  }
  
    
  /** explicit ND access */
  /* see if superclass methods work
  override def apply(i1:Int, i2:Int, i3:Int):Float = apply(Array(i1, i2, i3));
  override def apply(i1:Int, i2:Int, i3:Int, i4:Int):Float = apply(Array(i1, i2, i3, i4));
  override def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int):Float = apply(Array(i1, i2, i3, i4, i5));
  override def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int):Float = apply(Array(i1, i2, i3, i4, i5 ,i6));
*/

  /** ND slicing */
  /* see if superclass methods work

  override def apply(i1:IMat, i2:IMat, i3:IMat):GMat = apply(Array(i1, i2, i3), null);
  override def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat):GMat = apply(Array(i1, i2, i3, i4), null);
  override def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat):GMat = apply(Array(i1, i2, i3, i4, i5), null);
  override def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat):GMat = apply(Array(i1, i2, i3, i4, i5), null);
  
  */
  
  
  /** Type helpers for 2D updating with Ints */
  /* Should come from superclass
  override def update(i:Int, b:Double):FMat = update(i, b.toFloat);
  override def update(i:Int, b:Int):FMat = update(i, b.toFloat);   
    
  override def update(i:Int, j:Int, b:Double):FMat = update(i, j, b.toFloat);
  override def update(i:Int, j:Int, b:Int):FMat = update(i, j, b.toFloat);
  */
  
  /** ND single element updates */
  /* should come from superclass
  override def update(i1:Int, i2:Int, i3:Int, vv:Float):FMat = update(Array(i1, i2, i3), vv)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, vv:Float):FMat = update(Array(i1, i2, i3, i4), vv)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, vv:Float):FMat = update(Array(i1, i2, i3, i4, i5), vv)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, vv:Float):FMat = update(Array(i1, i2, i3, i4, i5, i6), vv)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, vv:Float):FMat = update(Array(i1, i2, i3, i4, i5, i6, i7), vv)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, i8:Int, vv:Float):FMat = update(Array(i1, i2, i3, i4, i5, i6, i7, i8), vv)
 
  override def update(i1:Int, i2:Int, i3:Int, vv:Double):FMat = update(Array(i1, i2, i3), vv.toFloat)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, vv:Double):FMat = update(Array(i1, i2, i3, i4), vv.toFloat)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, vv:Double):FMat = update(Array(i1, i2, i3, i4, i5), vv.toFloat)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, vv:Double):FMat = update(Array(i1, i2, i3, i4, i5, i6), vv.toFloat)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, vv:Double):FMat = update(Array(i1, i2, i3, i4, i5, i6, i7), vv.toFloat)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, i8:Int, vv:Double):FMat = update(Array(i1, i2, i3, i4, i5, i6, i7, i8), vv.toFloat)
 
  override def update(i1:Int, i2:Int, i3:Int, vv:Int):FMat = update(Array(i1, i2, i3), vv.toFloat)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, vv:Int):FMat = update(Array(i1, i2, i3, i4), vv.toFloat)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, vv:Int):FMat = update(Array(i1, i2, i3, i4, i5), vv.toFloat)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, vv:Int):FMat = update(Array(i1, i2, i3, i4, i5, i6), vv.toFloat)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, vv:Int):FMat = update(Array(i1, i2, i3, i4, i5, i6, i7), vv.toFloat)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, i8:Int, vv:Int):FMat = update(Array(i1, i2, i3, i4, i5, i6, i7, i8), vv.toFloat)
  */
  
  /* should be inherited from superclass
  
  override def update(i1:IMat, i2:IMat, i3:IMat, vv:FMat):GMat = update(Array(i1, i2, i3), GMat(vv));
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, vv:FMat):GMat = update(Array(i1, i2, i3, i4), GMat(vv));
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, vv:FMat):GMat = update(Array(i1, i2, i3, i4, i5), GMat(vv));
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, vv:FMat):GMat = update(Array(i1, i2, i3, i4, i5, i6), GMat(vv));

  override def update(iv:IMat, jv:IMat, b:Mat):FMat = update(Array(iv, jv), GMat(b));
  override def update(iv:IMat, j:Int, b:Mat):FMat = update(Array(iv, IMat.ielem(j)), GMat(b));
  override def update(i:Int, jv:IMat, b:Mat):FMat = update(Array(IMat.ielem(i), jv), GMat(b));
  
  override def update(i1:IMat, i2:IMat, i3:IMat, vv:Mat):GMat = update(Array(i1, i2, i3), GMat(vv));
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, vv:Mat):GMat = update(Array(i1, i2, i3, i4), GMat(vv));
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, vv:Mat):GMat = update(Array(i1, i2, i3, i4, i5), GMat(vv));
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, vv:Mat):GMat = update(Array(i1, i2, i3, i4, i5, i6), GMat(vv));
  */

  override def colslice(a:Int, b:Int):GMat =  colslice(a, b, null, 0);
    
  override def colslice(a:Int, b:Int, omat:Mat):GMat = colslice(a, b, omat, 0);
  
  override def colslice(a:Int, b:Int, omat:Mat, c:Int):GMat = {
	val newdims = _dims.clone;
    newdims(dims.length-1) = b-a+c;
    val out = if (omat.asInstanceOf[AnyRef] != null && omat.isInstanceOf[GMat] && omat.ncols >= b-a+c && omat.nrows == nrows) { 
      omat.asInstanceOf[GMat]
    } else { 
      GMat.newOrCheckGMat(newdims, omat, GUID, a, b, c, "colslice".##);
    }
    cudaMemcpy(out.pdata.withByteOffset(1L*c*nrows*Sizeof.FLOAT), pdata.withByteOffset(1L*a*nrows*Sizeof.FLOAT), 1L*(b-a)*nrows*Sizeof.FLOAT, cudaMemcpyDeviceToDevice);
    cudaStreamSynchronize(Mat.SyncMethod);
    val err = cudaGetLastError;
    if (err != 0) throw new RuntimeException("GMat colslice() error " + cudaGetErrorString(err));
    out
  }
  
  override def colslice(a:Int, b:Int, omat:Mat, c:Int, pb:Boolean):GMat = colslice(a, b, omat, c);
  
  val myGPU = SciFunctions.getGPU;
  
    /** reshaping */

  override def reshape(newdims:Int*):GMat = reshape(newdims.toArray)
  
  override def reshape(newdims:Array[Int]):GMat = {
    if (newdims.reduce(_*_) == length) {
      val out = GMat.newOrCheckGMat(newdims, null, GUID, ND.hashInts(newdims), "reshape".##);
      cudaMemcpy(out.pdata, pdata, 1L*llength*Sizeof.FLOAT, cudaMemcpyDeviceToDevice);
      cudaStreamSynchronize(Mat.SyncMethod);
      val err = cudaGetLastError;
      if (err != 0) throw new RuntimeException("GMat reshape() error " + cudaGetErrorString(err));
      out
    } else {
      throw new RuntimeException("GMat reshape total length doesnt match")
    }
  }
  
  override def reshapeView(newdims:Int*):GMat = reshapeView(newdims.toArray)
  
  override def reshapeView(newdims:Array[Int]):GMat = {
    if (newdims.reduce(_*_) == length) {
      val out = new GMat(newdims, pdata, llength);
      out.setGUID(MurmurHash3_x64_64(newdims.map(_.toLong) :+ GUID, "reshapeView".##));
      out
    } else {
      throw new RuntimeException("GMat reshapeView total length doesnt match")
    }
  }
  
  override def reshapeView(adims:IMat):GMat = reshapeView(adims.data);

  override def reshapeTrim(newdims:Int*):GMat = reshapeTrim(newdims.toArray)
  
  override def reshapeTrim(newdims:Array[Int]):GMat = {
    if (newdims.reduce(_*_) <= realsize) {
      val out = new GMat(newdims, pdata, realsize);
      out.setGUID(MurmurHash3_x64_64(newdims.map(_.toLong) :+ GUID, "reshapeTrim".##));
      out
    } else {
      throw new RuntimeException("GMat reshapeTrim total length too large")
    }
  }
  
  override def reshapeTrim(adims:IMat):GMat = reshapeTrim(adims.data);

  /** transpose */
  override def transpose(dims:Array[Int]):GMat = _transpose(MatFunctions.irow(dims));

  override def transpose(perm:IMat):GMat = _transpose(perm);
  
  override def _transpose(perm:IMat):GMat = { 
    val nd = _dims.length
    if (perm.length != nd) { 
      throw new RuntimeException("GMat transpose bad permutation ")
    }
    if (ND.isIdentity(perm)) {
    	this
    } else {
    	val xdims = MatFunctions.irow(_dims);
    	val iperm = MatFunctions.invperm(perm);
    	val pdims = xdims(perm).data;
    	var out = GMat.newOrCheckGMat(pdims, null, GUID, ND.hashInts(pdims), "transpose".##);
    	var out2 = GMat.newOrCheckGMat(pdims, null, GUID, ND.hashInts(pdims), "transpose1".##);
    	cudaMemcpy(out.pdata, pdata, 1L*Sizeof.FLOAT*length, cudaMemcpyDeviceToDevice);
    	cudaStreamSynchronize(Mat.SyncMethod);
    	for (i <- (nd - 1) until 0 by -1) { 
    		if (iperm(i) != i) { 
    			val (d1, d2, d3) = ND.getDims(i, iperm, xdims);
    			if (d1 > 1 && d2 > 1) { 
    				//         println("spermute %d %d %d" format (d1,d2,d3))
    				CUMAT.spermute(d1, d2, d3, out.pdata, out2.pdata);
    				val tmp = out2;
    				out2 = out;
    				out = tmp;
    			}
    			ND.rotate(i, iperm, xdims);
    		} 
    	}
    	out;
    }
  }
  
  override def fromNHWCtoNCHW:GMat = {
    if (dims.length != 4) throw new RuntimeException("fromNHWCtoNCHW ndims must be 4");
    transpose(MatFunctions.irow(1,2,0,3)).reshapeView(dims);
  }
  
  override def fromNCHWtoNHWC:GMat = {
    if (dims.length != 4) throw new RuntimeException("fromNCHWtoNHWC ndims must be 4");
    reshapeView(MatFunctions.irow(dims(1), dims(2), dims(0), dims(3))).transpose(MatFunctions.irow(2,0,1,3));
  }
  
  override def clear = {
  	cudaMemset(pdata, 0, Sizeof.FLOAT*length)
  	cudaStreamSynchronize(Mat.SyncMethod)
  	this    
  }
  
  override def set(v:Float):GMat = {
    CUMAT.setval(pdata, v, length)
    cudaStreamSynchronize(Mat.SyncMethod);
    this
  }
  
 
  override def zeros(nr:Int, nc:Int) = GMat.zeros(nr, nc);
  
  override def zeros(nr:Int, nc:Int, nnz:Int) = GMat.zeros(nr, nc);
  
  override def zeros(dims:IMat):GMat = {
    GMat.zeros(dims)
  }
  
  override def ones(nr:Int, nc:Int) = GMat.ones(nr, nc);
  
  override def ones(dims:IMat) = GMat.ones(dims);
  
  override def zero = GMat.zeros(1, 1);
  
  override def one = GMat.ones(1, 1);
  
  override def izeros(m:Int, n:Int) = {
    GIMat.izeros(m,n)
  }
  
  override def izeros(dims:IMat) = {
    GIMat.izeros(dims);
  }
  
  override def iones(m:Int, n:Int) = {
    GIMat.iones(m,n)
  }
  
  override def iones(dims:IMat) = {
    GIMat.iones(dims)
  }
  
  def horzcat(aa:FMat, omat:Mat) = {
	  val a = GMat(aa);
	  if (nrows != a.nrows)
		  throw new RuntimeException("GMat \\ row dims not equal");
	  val out = GMat.newOrCheckGMat(nrows, ncols+a.ncols, omat, GUID, a.GUID, "horzcat".##);
	  cudaMemcpy(out.pdata, pdata, 1L*length*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyDeviceToDevice);
	  cudaStreamSynchronize(Mat.SyncMethod);
	  var err = cudaGetLastError;
	  if (err != 0) throw new RuntimeException("GMat horzcat() error " + cudaGetErrorString(err));
	  cudaMemcpy(out.pdata.withByteOffset(1L*length*Sizeof.FLOAT), a.pdata, 1L*a.length*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
	  cudaStreamSynchronize(Mat.SyncMethod);
	  err = cudaGetLastError;
	  if (err != 0) throw new RuntimeException("GMat horzcat() error " + cudaGetErrorString(err));
	  out
  }
  
  def vertcat(aa:FMat, omat:Mat) = {
	  val a = GMat(aa);
	  if (ncols != a.ncols)
		  throw new RuntimeException("GMat on row dims not equal");
	  val out = GMat.newOrCheckGMat(nrows+a.nrows, ncols, omat, GUID, a.GUID, "vertcat".##);
	  cudaMemcpy2D(out.pdata, 1L*out.nrows*Sizeof.FLOAT, pdata, 1L*nrows*Sizeof.FLOAT, 1L*nrows*Sizeof.FLOAT, 1L*ncols, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
	  cudaStreamSynchronize(Mat.SyncMethod);
	  var err = cudaGetLastError;
	  if (err != 0) throw new RuntimeException("GMat vertcat() error " + cudaGetErrorString(err));
	  cudaMemcpy2D(out.pdata.withByteOffset(1L*nrows*Sizeof.FLOAT), 1L*out.nrows*Sizeof.FLOAT, a.pdata, 1L*a.nrows*Sizeof.FLOAT, 1L*a.nrows*Sizeof.FLOAT,  1L*a.ncols, cudaMemcpyKind.cudaMemcpyDeviceToDevice);
	  cudaStreamSynchronize(Mat.SyncMethod);
	  err = cudaGetLastError;
	  if (err != 0) throw new RuntimeException("GMat vertcat() error " + cudaGetErrorString(err));
	  out
  }
  
  def getHandle = {
    val igpu = Array(0);
    jcuda.runtime.JCuda.cudaGetDevice(igpu)
    GFunctions.cublasHandles(igpu(0)).asInstanceOf[cublasHandle];
  }

  def GMult(aa:FMat, oldmat:Mat):GMat = {
    val a = GMat(aa);
    if (ncols == 1 && nrows == 1) {
      val out = GMat.newOrCheckGMat(a.nrows, a.ncols, oldmat, GUID, a.GUID, "GMult1".##);
      Mat.nflops += 1L * a.length;
      val err = CUMAT.applyop(pdata, nrows, ncols, a.pdata, a.nrows, a.ncols, out.pdata, GMat.BinOp.op_mul);
      if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.applyop " + cudaGetErrorString(err))};
      out
    } else if (a.ncols == 1 && a.nrows == 1) {
      val out = GMat.newOrCheckGMat(nrows, ncols, oldmat, GUID, a.GUID, "GMult2".##);
      Mat.nflops += 1L * length;
      val err = CUMAT.applyop(pdata, nrows, ncols, a.pdata, a.nrows, a.ncols, out.pdata, GMat.BinOp.op_mul);
      if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.applyop " + cudaGetErrorString(err))}
      out;
    } else if (ncols == a.nrows) {
    	val out = GMat.newOrCheckGMat(nrows, a.ncols, oldmat, GUID, a.GUID, "GMult".##);
    	Mat.nflops += 2L * length * a.ncols;
    	if (nrows == 1) {
    		//        cublasSgemv('t', a.nrows, a.ncols, 1.0f, a.pdata, nrows, pdata, 1, 0f, out.pdata, 1)
    		out.clear;
    		CUMAT.dmv(a.pdata, a.nrows, a.ncols, pdata, out.pdata, 1);
 //   		Thread.sleep(0);
    		cudaStreamSynchronize(Mat.SyncMethod);
    		val err = cudaGetLastError;
    		if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.dmv " + cudaGetErrorString(err))}
    	} else if (a.ncols == 1) {
    		//        cublasSgemv('n', nrows, ncols, 1.0f, pdata, nrows, a.pdata, 1, 0f, out.pdata, 1)
    		out.clear;
    		CUMAT.dmv(pdata, nrows, ncols, a.pdata, out.pdata, 0);
//    		Thread.sleep(0)
    		cudaStreamSynchronize(Mat.SyncMethod);
    		val err = cudaGetLastError;
    		if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.dmv " + cudaGetErrorString(err))}
    	} else {
    		cublasSgemm(getHandle, cublasOperation.CUBLAS_OP_N, cublasOperation.CUBLAS_OP_N, nrows, a.ncols, ncols, 
    				GMat.pONE, pdata, nrows, a.pdata, a.nrows, GMat.pZERO, out.pdata, nrows);
    		if (length * a.ncols > GMat.multYieldSize) Thread.`yield`;
    		cudaStreamSynchronize(Mat.SyncMethod);
    		val err = cudaGetLastError;
    		if (err != 0) {
    			println("device is %d" format SciFunctions.getGPU);
    			throw new RuntimeException("Cublas error in * "+err);
    		}
    	}

    	out 
    } else throw new RuntimeException("dimensions mismatch (%d %d), (%d %d)" format (nrows, ncols, a.nrows, a.ncols));
  }
  
  override def mult(bb:SMat, c:FMat, at:Boolean, bt:Boolean):FMat = {
    val b = GSMat(bb);
    (at, bt) match {
      case (false, false) => GSMult(b, c);
      case (false, true) => GSMultT(b, c);
      case _ => throw new RuntimeException("mult unsupported options SMat, FMat %b %b" format (at, bt));
    }
  }
  
  override def mult(bb:FMat, c:FMat, at:Boolean, bt:Boolean):FMat = {
    val b = GMat(bb);
    (at, bt) match {
      case (false, false) => GMult(b, c);
      case (false, true) => GMultT(b, c);
      case (true, false) => GTMult(b, c);
      case _ => throw new RuntimeException("mult unsupported options FMat, FMat %b %b" format (at, bt));
    }
  }
  
  override def mult(b:Mat, c:Mat, at:Boolean, bt:Boolean):Mat = {
    (b, c) match {
      case (bb:FMat, cc:FMat) => mult(bb, cc, at, bt);
      case (bb:SMat, cc:FMat) => mult(bb, cc, at, bt);
      case _ => throw new RuntimeException("mult unsupported types %s %s" format (b.mytype, c.mytype));
    }
    c
  }
  
  override def mult(b:Mat, c:Mat):Mat = mult(b, c, false, false);
  
   
  
  override def madd(bb:FMat, cc:FMat, at:Boolean, bt:Boolean):GMat = {
	  val b = GMat(bb);
	  val c = GMat(cc);
	  val (arows, acols, atrans) = if (at) (ncols, nrows, cublasOperation.CUBLAS_OP_T) else (nrows, ncols, cublasOperation.CUBLAS_OP_N);
	  val (brows, bcols, btrans) = if (bt) (b.ncols, b.nrows, cublasOperation.CUBLAS_OP_T) else (b.nrows, b.ncols, cublasOperation.CUBLAS_OP_N);
	  if (acols != brows || arows != c.nrows || bcols != c.ncols) {
		  throw new RuntimeException("madd bad dimensions (%d %d) (%d %d) (%d %d)" format (arows, acols, brows, bcols, c.nrows, c.ncols));
	  }
	  Mat.nflops += 2L * arows * bcols * acols;
	  cublasSgemm(getHandle, atrans, btrans,	arows, bcols, acols, GMat.pONE, pdata, nrows, b.pdata, b.nrows, GMat.pONE, c.pdata, c.nrows);
	  if (1L * arows * bcols * acols > GMat.multYieldSize) Thread.`yield`;
	  cudaStreamSynchronize(Mat.SyncMethod); 
	  c
  }
  
  override def madd(b:FMat, c:FMat):GMat = madd(b, c, false, false);
  
  override def madd(b:SMat, c:FMat, bt:Boolean, ct:Boolean):GMat = {
    (bt, ct) match {
      case (false, false) => madd(b, c);
      case (false, true) => maddT(b, c);
      case _ => throw new RuntimeException("madd unsupported options GSMat, GMat %b %b" format (bt, ct));
    }
  }
  
  override def madd(b:FMat,c:TMat):TMat = madd(b,c,false,false)

  import BIDMat.IMatWildcard
 
  override def madd(bb:FMat,c:TMat,at:Boolean,bt:Boolean):TMat = {
    val b = GMat(bb);
    for (i <- 0 until c.tiles.length) {
      val m = c.tiles(i);
    	if (!at) {
    		if (!bt) {
    			tileMult(m.nrows,m.ncols,ncols,c.y(i),0,b,0,c.x(i),m,0,0);
    		}	else {
    			tileMultNT(m.nrows,m.ncols,ncols,c.y(i),0,b,c.x(i),0,m,0,0);
    		}
    	} else {
    		if (!bt) {
    			tileMultTN(m.nrows,m.ncols,nrows,0,c.y(i),b,0,c.x(i),m,0,0);
    		}	else {
    			tileMultTT(m.nrows,m.ncols,nrows,0,c.y(i),b,c.x(i),0,m,0,0);
    		}
    	}
    }
    c;
  }
  
  override def madd(b:Mat, c:Mat, at:Boolean, bt:Boolean):Mat = {
  	(b, c) match {
  	case (bb:FMat, cc:FMat) => madd(bb, cc, at, bt);
  	case (bb:SMat, cc:FMat) => madd(bb, cc, at, bt);
  	case (bb:FMat,cc:TMat) => madd(bb,cc,at,bt);
  	case (bb:SMat,cc:TMat) => madd(bb,cc,at,bt);
  	case _ => throw new RuntimeException("madd unsupported types %s %s" format (b.mytype, c.mytype));
  	}
  	c
  }
  
  override def madd(b:Mat, c:Mat):Mat = madd(b, c, false, false);
  
  override def blockmult(bb:FMat, cc:FMat, nblocks:Int, at:Boolean, bt:Boolean, cfact:Float):FMat = {
    val b = GMat(bb);
    val c = GMat(cc);
 
    val (anrows, ancols) = if (dims.length >= 3) {
      (dims(0), dims(1))
    } else {
      (nrows/nblocks, ncols)
    }
    val (bnrows, bncols) = if (b.dims.length >= 3) {
      (b.dims(0), b.dims(1))
    } else {
      (b.nrows/nblocks, b.ncols)
    }
    val (cnrows,cncols) = if (c.dims.length >= 3) {
      (c.dims(0), c.dims(1))
    } else {
      (c.nrows/nblocks, c.ncols)
    }
    blockGemm(if (at) 1 else 0, if (bt) 1 else 0, cnrows, cncols, if (at) anrows else ancols, 1f, 0, anrows, anrows*ancols,
    		b, 0, bnrows, bnrows*bncols, cfact, c, 0, cnrows, cnrows*cncols, nblocks);
    c
  }
  
  override def blockmult(b:FMat, c:FMat, nblocks:Int, at:Boolean, bt:Boolean):FMat = blockmult(b, c, nblocks, at, bt, 0f);
 
  override def blockmult(b:Mat, c:Mat, nblocks:Int):Mat = blockmult(b, c, nblocks, false, false);
  
  override def blockmult(b:Mat, c:Mat, nblocks:Int, at:Boolean, bt:Boolean):Mat = {
    (b, c) match {
      case (bb:FMat, cc:FMat) => blockmult(bb, cc, nblocks, at, bt);
      case _ => throw new RuntimeException("blockmult unsupported types %s %s" format (b.mytype, c.mytype));
    }
    c
  }
  
  override def blockmult2(bb:FMat, cc:FMat, nblocks:Int, at:Boolean, bt:Boolean, cfact:Float):FMat = {
    val b = GMat(bb);
    val c = GMat(cc);
 
    val anrows = dims(0)
    val astep = dims(1)
    val ancols = dims(2)
    val bnrows = b.dims(0)
    val bstep = b.dims(1)
    val bncols = b.dims(2)
    val cnrows = c.dims(0)
    val cstep = c.dims(1)
    val cncols = c.dims(2)
    if (dims.length == 3) {
    	blockGemm(if (at) 1 else 0, if (bt) 1 else 0, cnrows, cncols, if (at) anrows else ancols, 1f, 0, anrows*astep, anrows,
    			b, 0, bnrows*bstep, bnrows, cfact, c, 0, cnrows*cstep, cnrows, nblocks);
    } else {
      val reps2 = dims.data.slice(3, dims.length).reduce(_*_);
      blockGemm4D(if (at) 1 else 0, if (bt) 1 else 0, cnrows, cncols, if (at) anrows else ancols, 1f, 0, anrows*astep, anrows, anrows*ancols*astep,
    			b, 0, bnrows*bstep, bnrows, bnrows*bncols*bstep, cfact, c, 0, cnrows*cstep, cnrows, cnrows*cncols*cstep, nblocks, reps2);
    }
    c
  }
  
  override def blockmult2(b:FMat, c:FMat, nblocks:Int, at:Boolean, bt:Boolean):FMat = blockmult2(b, c, nblocks, at, bt, 0f);
 
  override def blockmult2(b:Mat, c:Mat, nblocks:Int):Mat = blockmult2(b, c, nblocks, false, false);
  
  override def blockmult2(b:Mat, c:Mat, nblocks:Int, at:Boolean, bt:Boolean):Mat = {
    (b, c) match {
      case (bb:FMat, cc:FMat) => blockmult2(bb, cc, nblocks, at, bt);
      case _ => throw new RuntimeException("blockmult2 unsupported types %s %s" format (b.mytype, c.mytype));
    }
    c
  }
  
  override def blockmadd(b:FMat, c:FMat, nblocks:Int, at:Boolean, bt:Boolean):FMat = blockmult(b, c, nblocks, at, bt, 1f);
 
  override def blockmadd(b:Mat, c:Mat, nblocks:Int):Mat = blockmadd(b, c, nblocks, false, false);
  
  override def blockmadd(b:Mat, c:Mat, nblocks:Int, at:Boolean, bt:Boolean):Mat = {
    (b, c) match {
      case (bb:FMat, cc:FMat) => blockmadd(bb, cc, nblocks, at, bt);
      case _ => throw new RuntimeException("blockmadd unsupported types %s %s" format (b.mytype, c.mytype));
    }
    c
  }
  
  override def blockmadd2(b:FMat, c:FMat, nblocks:Int, at:Boolean, bt:Boolean):FMat = blockmult2(b, c, nblocks, at, bt, 1f);
 
  override def blockmadd2(b:Mat, c:Mat, nblocks:Int):Mat = blockmadd2(b, c, nblocks, false, false);
  
  override def blockmadd2(b:Mat, c:Mat, nblocks:Int, at:Boolean, bt:Boolean):Mat = {
    (b, c) match {
      case (bb:FMat, cc:FMat) => blockmadd2(bb, cc, nblocks, at, bt);
      case _ => throw new RuntimeException("blockmadd2 unsupported types %s %s" format (b.mytype, c.mytype));
    }
    c
  }

  
  def GMultT(aa:FMat, oldmat:Mat):GMat = {
    val a = GMat(aa);
    if (ncols == a.ncols) {
      val out = GMat.newOrCheckGMat(nrows, a.nrows, oldmat, GUID, a.GUID, "GMultT".##)
      Mat.nflops += 2L * length * a.nrows
      cublasSgemm(getHandle, cublasOperation.CUBLAS_OP_N, cublasOperation.CUBLAS_OP_T, nrows, a.nrows, ncols, 
          GMat.pONE, pdata, nrows, a.pdata, a.nrows, GMat.pZERO, out.pdata, nrows);
      if (1L * length * a.nrows > GMat.multYieldSize) Thread.`yield`;
      cudaStreamSynchronize(Mat.SyncMethod);
      val err = cudaGetLastError
      if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("Cublas error in xT " + cudaGetErrorString(err))
      }
      out
    } else throw new RuntimeException("dimensions mismatch")
  }
  
  /**
   * Tile multiplication of a(this) * b into c. Tile coordinates are (r,c,height,width):
   * a: (aroff, acoff, nr, kk)
   * b: (broff, bcoff, kk, nc)
   * c: (croff, ccoff, nr, nc)
   * 
   * Note: c is not cleared by the kernel, and the result is added to it. 
   */
  
  def tileMult(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:GMat, broff:Int, bcoff:Int, c:GMat, croff:Int, ccoff:Int):GMat = {
    if (aroff < 0 || acoff < 0 || broff < 0 || bcoff < 0 || croff < 0 || ccoff < 0 || nr < 0 || nc < 0 || kk < 0) {
    	throw new RuntimeException("tileMul: cant have negative offsets or dimensions");
    } else if (aroff + nr > nrows || acoff + kk > ncols || broff + kk > b.nrows || bcoff + nc > b.ncols || croff + nr > c.nrows || ccoff + nc > c.ncols) {
      throw new RuntimeException("tileMult: tile strays outside matrix dimensions");
    } else {
      Mat.nflops += 2L * nr * nc * kk;
    	cublasSgemm(getHandle, cublasOperation.CUBLAS_OP_N, cublasOperation.CUBLAS_OP_N, nr, nc, kk, GMat.pONE, 
    	    pdata.withByteOffset(Sizeof.FLOAT.toLong*(aroff+acoff*nrows)), nrows, 
    	    b.pdata.withByteOffset(Sizeof.FLOAT.toLong*(broff+bcoff*b.nrows)), b.nrows, GMat.pONE, 
      		c.pdata.withByteOffset(Sizeof.FLOAT.toLong*(croff+ccoff*c.nrows)), c.nrows);
    	cudaStreamSynchronize(Mat.SyncMethod);
      c;
    }
  }
  
  def tileMultNT(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:GMat, broff:Int, bcoff:Int, c:GMat, croff:Int, ccoff:Int):GMat = {
    if (aroff < 0 || acoff < 0 || broff < 0 || bcoff < 0 || croff < 0 || ccoff < 0 || nr < 0 || nc < 0 || kk < 0) {
    	throw new RuntimeException("tileMultNT: cant have negative offsets or dimensions");
    } else if (aroff + nr > nrows || acoff + kk > ncols || broff + nc > b.nrows || bcoff + kk > b.ncols || croff + nr > c.nrows || ccoff + nc > c.ncols) {
      throw new RuntimeException("tileMultNT: tile strays outside matrix dimensions");
    } else {
    	Mat.nflops += 2L * nr * nc * kk;
    	cublasSgemm(getHandle, cublasOperation.CUBLAS_OP_N, cublasOperation.CUBLAS_OP_T,	nr, nc, kk, GMat.pONE, 
    	    pdata.withByteOffset(Sizeof.FLOAT.toLong*(aroff+acoff*nrows)), nrows, 
    	    b.pdata.withByteOffset(Sizeof.FLOAT.toLong*(broff+bcoff*b.nrows)), b.nrows, GMat.pONE, 
      		c.pdata.withByteOffset(Sizeof.FLOAT.toLong*(croff+ccoff*c.nrows)), c.nrows);
    	cudaStreamSynchronize(Mat.SyncMethod);
      c;
    }
  }
  
  def tileMultTN(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:GMat, broff:Int, bcoff:Int, c:GMat, croff:Int, ccoff:Int):GMat = {
    if (aroff < 0 || acoff < 0 || broff < 0 || bcoff < 0 || croff < 0 || ccoff < 0 || nr < 0 || nc < 0 || kk < 0) {
    	throw new RuntimeException("tileMultTN: cant have negative offsets or dimensions");
    } else if (aroff + kk > nrows || acoff + nr > ncols || broff + kk > b.nrows || bcoff + nc > b.ncols || croff + nr > c.nrows || ccoff + nc > c.ncols) {
      throw new RuntimeException("tileMultTN: tile strays outside matrix dimensions");
    } else {
    	Mat.nflops += 2L * nr * nc * kk;
    	cublasSgemm(getHandle, cublasOperation.CUBLAS_OP_T, cublasOperation.CUBLAS_OP_N, nr, nc, kk, GMat.pONE, 
    	    pdata.withByteOffset(Sizeof.FLOAT.toLong*(aroff+acoff*nrows)), nrows, 
    	    b.pdata.withByteOffset(Sizeof.FLOAT.toLong*(broff+bcoff*b.nrows)), b.nrows, GMat.pONE, 
      		c.pdata.withByteOffset(Sizeof.FLOAT.toLong*(croff+ccoff*c.nrows)), c.nrows);
    	cudaStreamSynchronize(Mat.SyncMethod);
      c;
    }
  }
  
  def tileMult(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:GSMat, broff:Int, bcoff:Int, c:GMat, croff:Int, ccoff:Int):GMat = {
    if (aroff < 0 || acoff < 0 || broff < 0 || bcoff < 0 || croff < 0 || ccoff < 0 || nr < 0 || nc < 0 || kk < 0) {
    	throw new RuntimeException("tileMult: cant have negative offsets or dimensions");
    } else if (aroff + nr > nrows || acoff + kk > ncols || broff + kk > b.nrows || bcoff + nc > b.ncols || croff + nr > c.nrows || ccoff + nc > c.ncols) {
      throw new RuntimeException("tileMult: tile strays outside matrix dimensions");
    } else {
    	Mat.nflops += 2L * nr * b.nnz;
    	val err = CUMAT.dsmultTile(nr, nc, kk, b.nnz,  
    			pdata.withByteOffset(Sizeof.FLOAT.toLong*(aroff+acoff*nrows)), nrows, 
    	    b.pdata, b.pir, b.pic, broff, bcoff, 
      		c.pdata.withByteOffset(Sizeof.FLOAT.toLong*(croff+ccoff*c.nrows)), c.nrows, 0);
    	cudaStreamSynchronize(Mat.SyncMethod);
    	if (err != 0) {
    		throw new RuntimeException("CUMAT.tileMult error " + cudaGetErrorString(err))
    	}
      c;
    }
  }
  
  def tileMultNT(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:GSMat, broff:Int, bcoff:Int, c:GMat, croff:Int, ccoff:Int):GMat = {
    if (aroff < 0 || acoff < 0 || broff < 0 || bcoff < 0 || croff < 0 || ccoff < 0 || nr < 0 || nc < 0 || kk < 0) {
    	throw new RuntimeException("tileMultNT: cant have negative offsets or dimensions");
    } else if (aroff + nr > nrows || acoff + kk > ncols || broff + nc > b.nrows || bcoff + kk > b.ncols || croff + nr > c.nrows || ccoff + nc > c.ncols) {
      throw new RuntimeException("tileMultNT: tile strays outside matrix dimensions");
    } else {
    	Mat.nflops += 2L * nr * b.nnz * kk / b.ncols;
    	val err = CUMAT.dsmultTile(nr, nc, kk, b.nnz,  
    			pdata.withByteOffset(Sizeof.FLOAT.toLong*(aroff+acoff*nrows)), nrows, 
    	    b.pdata, b.pir, b.pic, broff, bcoff, 
      		c.pdata.withByteOffset(Sizeof.FLOAT.toLong*(croff+ccoff*c.nrows)), c.nrows, 1);
    	cudaStreamSynchronize(Mat.SyncMethod);
    	if (err != 0) {
    		throw new RuntimeException("CUMAT.tileMultT error " + cudaGetErrorString(err))
    	}
      c;
    }
  }
  
  override def tileMult(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:Mat, broff:Int, bcoff:Int, c:Mat, croff:Int, ccoff:Int):GMat = {
    (b, c) match {
      case (sb:GSMat, fc:GMat) => tileMult(nr, nc, kk, aroff, acoff, sb, broff, bcoff, fc, croff, ccoff);
      case (fb:GMat, fc:GMat) => tileMult(nr, nc, kk, aroff, acoff, fb, broff, bcoff, fc, croff, ccoff);
      case _ => throw new RuntimeException("tileMult couldnt match matrix types")
    }
  }
  
  override def tileMultNT(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:Mat, broff:Int, bcoff:Int, c:Mat, croff:Int, ccoff:Int):GMat = {
    (b, c) match {
      case (sb:GSMat, fc:GMat) => tileMultNT(nr, nc, kk, aroff, acoff, sb, broff, bcoff, fc, croff, ccoff);
      case (fb:GMat, fc:GMat) => tileMultNT(nr, nc, kk, aroff, acoff, fb, broff, bcoff, fc, croff, ccoff);
      case _ => throw new RuntimeException("tileMultT couldnt match matrix types")
    }
  }
  
  override def tileMultTN(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:Mat, broff:Int, bcoff:Int, c:Mat, croff:Int, ccoff:Int):GMat = {
    (b, c) match {
//      case (sb:GSMat, fc:GMat) => tileMultTN(nr, nc, kk, aroff, acoff, sb, broff, bcoff, fc, croff, ccoff);
      case (fb:GMat, fc:GMat) => tileMultTN(nr, nc, kk, aroff, acoff, fb, broff, bcoff, fc, croff, ccoff);
      case _ => throw new RuntimeException("tileMultT couldnt match matrix types")
    }
  }
  
  def GTMult(aa:FMat, oldmat:Mat):GMat = {
    val a= GMat(aa);
    if (nrows == a.nrows) {
      val out = GMat.newOrCheckGMat(ncols, a.ncols, oldmat, GUID, a.GUID, "GMultT".##)
      Mat.nflops += 2L * length * a.ncols
      cublasSgemm(getHandle, cublasOperation.CUBLAS_OP_T, cublasOperation.CUBLAS_OP_N, ncols, a.ncols, nrows, 
          GMat.pONE, pdata, nrows, a.pdata, a.nrows, GMat.pZERO, out.pdata, out.nrows);
      if (1L * length * a.ncols > GMat.multYieldSize) Thread.`yield`;
      cudaStreamSynchronize(Mat.SyncMethod);
      val err = cudaGetLastError
      if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("Cublas error in Tx " + cudaGetErrorString(err))
      }
      out
    } else throw new RuntimeException("dimensions mismatch")
  }
  
  def GSMult(aa:SMat, oldmat:Mat):GMat = {
    val a = GSMat(aa);
    if (ncols != a.nrows) {
      throw new RuntimeException("GSMult dimensions mismatch (%d %d) (%d %d)" format (nrows, ncols, a.nrows, a.ncols))
    }
    val out = GMat.newOrCheckGMat(nrows, a.ncols, oldmat, GUID, a.GUID, "GSMult".##);
    out.clear;
    madd(a, out);
  }
  
  override def madd(aa:SMat, oo:FMat):GMat = {
    val a = GSMat(aa);
    val out = GMat(oo);
    if (ncols != a.nrows || nrows != out.nrows || a.ncols != out.ncols) {
      throw new RuntimeException("GSMadd dimensions mismatch (%d %d) (%d %d) (%d %d)" format (nrows, ncols, a.nrows, a.ncols, out.nrows, out.ncols))
    }
    Mat.nflops += 2L * nrows * a.nnz;  
/*      if (nrows == 1) {                    // Alas, throws "too many resources requested for launch" with large a.nrows
      	val handle = GSMat.getHandle       // Also gives erroneous values
      	val descra = GSMat.getDescr
        var err = JCusparse.cusparseScsrmv(handle, cusparseOperation.CUSPARSE_OPERATION_NON_TRANSPOSE,
        		ncols, a.ncols, 1.0f, descra,	a.pdata, a.jc, a.ir, pdata, 0, out.pdata)
        cudaStreamSynchronize(Mat.SyncMethod)()
        if (err == 0) err = cudaGetLastError
        if (err != 0) {
        	println("device is %d" format SciFunctions.getGPU)
        	throw new RuntimeException("Cuda error in GSMult " + cudaGetErrorString(err))
        }
      } else { */
    val err = CUMAT.dsmult(nrows, a.ncols, a.nnz, pdata, a.pdata, a.pir, a.pic, out.pdata);
    if (err != 0) throw new RuntimeException("GMult: CUDA kernel error in CUMAT.dsmult " + cudaGetErrorString(err));
    //      }
    out;
  }
  
  def GSMultT(aa:SMat, oldmat:Mat):GMat = {
    val a = GSMat(aa);
  	if (ncols != a.ncols) { 
  		throw new RuntimeException("GSMult dimensions mismatch (%d %d) (%d %d)" format (nrows, ncols, a.ncols, a.nrows))
  	}
  	val out = GMat.newOrCheckGMat(nrows, a.nrows, oldmat, GUID, a.GUID, "GSMultT".##);
  	out.clear;
  	maddT(a, out);
  }
  
  override def maddT(aa:SMat, oo:FMat):GMat = {
    val a = GSMat(aa);
    val out = GMat(oo);
    if (ncols != a.ncols || nrows != out.nrows || a.nrows != out.ncols) {
      throw new RuntimeException("GSMadd dimensions mismatch (%d %d) (%d %d) (%d %d)" format (nrows, ncols, a.nrows, a.ncols, out.nrows, out.ncols))
    }
    Mat.nflops += 2L * nrows * a.nnz;
    val err = CUMAT.dsmultT(nrows, a.ncols, a.nnz, pdata, a.pdata, a.pir, a.pic, out.pdata);
    if (err != 0) throw new RuntimeException("GMult: CUDA kernel error in CUMAT.dsmultT " + cudaGetErrorString(err));
    out
  }
  
  def GMST(aa:FMat, oldmat:Mat):GMat = {
    val a = GMat(aa);
    if (ncols == a.ncols) {
      val out = GMat.newOrCheckGMat(nrows, a.nrows, oldmat, GUID, a.GUID, "GMST".##)
      Mat.nflops += 2L * nrows * a.nrows * ncols
      out.clear
      val err = CUMAT.maxsumx(pdata, nrows, a.pdata, a.nrows, out.pdata, nrows, ncols, nrows, a.nrows);
      if (err != 0) throw new RuntimeException("GMult: CUDA kernel error in CUMAT.maxsumx " + cudaGetErrorString(err))
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }
  
  override def kron(aa:FMat, oldmat:Mat):GMat = {
    val a = GMat(aa);
    val out = GMat.newOrCheckGMat(nrows * a.nrows, ncols * a.ncols, oldmat, GUID, a.GUID, "kron".##);
    Mat.nflops += 1L * out.nrows * out.ncols;
    val err = CUMAT.kron(pdata, a.pdata, out.pdata, nrows, ncols, a.nrows, a.ncols);
    if (err != 0) throw new RuntimeException("kron: CUDA kernel error in CUMAT.kron " + cudaGetErrorString(err));
    out;
  }
  
  def gOp(aa:FMat, oldmat:Mat, op:Int):GMat = {
    val a = GMat(aa);
    val (nr, nc, nra, nca) = ND.compatibleGDims(_dims, aa._dims, "DenseMat Op");
    val dims = ND.maxDims(_dims, aa._dims);
    val out = GMat.newOrCheckGMat(dims, oldmat, GUID, aa.GUID, op.hashCode);
    Mat.nflops += scala.math.max(length, a.length);
    val err = CUMAT.applyop(pdata, nr, nc, a.pdata, nra, nca, out.pdata, op);
    if (err != 0) {throw new RuntimeException("CUDA kernel error %d in CUMAT.applyop"  format err)}
    out
  }
  
  override def dot(aa:FMat, oldmat:Mat):GMat = {
  		val a = GMat(aa);
  		ND.checkDims("dot", dims, a.dims);
  		val odims = IMat.iones(1, dims.length); 
  		odims(dims.length-1) = a.ncols;
  		val out = GMat.newOrCheckGMat(odims, oldmat, GUID, a.GUID, "dot".##);
  		Mat.nflops += 2L * length;
  		val err = CUMAT.reducebin1op(nrows, ncols, pdata, a.pdata, out.pdata, op_mul, op_add);
  		if (err != 0) {throw new RuntimeException("GMat dot: CUDA kernel error in CUMAT.reducebin1op " + cudaGetErrorString(err))}
  		out;
  }
  
  override def dot (a:FMat):GMat = dot(a, null)
  
  override def dotr (aa:FMat, oldmat:Mat):GMat = {
    val a = GMat(aa);
    ND.checkDims("dotr", dims, a.dims);
    val odims = a.dims.copy;
    odims(odims.length-1) = 1;
    val out = GMat.newOrCheckGMat(odims, oldmat, GUID, a.GUID, "dotr".##);
    Mat.nflops += 2L * length;
    val err = CUMAT.reducebin2op(nrows, ncols, pdata, a.pdata, out.pdata, op_mul, op_add);
    if (err != 0) {throw new RuntimeException("GMat dotr: CUDA kernel error in CUMAT.reducebin2op " + cudaGetErrorString(err))}
    out;
  }
  
  override def dotr (a:FMat):GMat = dotr(a, null)
  
  override def ddot (a:Mat):Double = 
  	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("ddot dims not compatible")
  	} else {
  	  a match {
  	  case aa:GMat => {
  	    val result = Array(0f)
  	    cublasSdot(getHandle, length, pdata, 1, aa.pdata, 1, Pointer.to(result))
  	  	cudaStreamSynchronize(Mat.SyncMethod)
  	  	val err = cudaGetLastError
  	  	if (err != 0) {
  	  		println("device is %d" format SciFunctions.getGPU)
  	  		throw new RuntimeException("Cublas error in ddot " + cudaGetErrorString(err))
  	  	}
  	  result(0);
  	  }
  	  }
  	}
  
  override def ddot(a:FMat):Double = {
    val am:Mat = a;
    ddot(am);
  }
  
  def reduceOp(oldmat:Mat, dir:Int, initval:Float, op:Int):GMat = {
    if (dir == 1 || (dir == 0 && nrows > 1)) {
      val out = GMat.newOrCheckGMat(1, ncols, oldmat, GUID, 1, op) 
      out.clear
      val err = CUMAT.reduce1op(nrows, ncols, pdata, out.pdata, initval, op)
      if (err != 0) {throw new RuntimeException("CUDA kernel error in CUMAT.reduce1op " + cudaGetErrorString(err))}
      Mat.nflops += length
      out
    } else if (dir == 2 || dir == 0) {
      val out = GMat.newOrCheckGMat(nrows, 1, oldmat, GUID, 2, op)  
      out.clear
      val err = CUMAT.reduce2op(nrows, ncols, pdata, out.pdata, initval, op)
      if (err != 0) {throw new RuntimeException("CUDA kernel error in CUMAT.reduce2op " + cudaGetErrorString(err))}
      Mat.nflops += length
      out
    } else {
      throw new RuntimeException("dimension must be 1 or 2")
    }
  }

  def toFMat(a:Mat):FMat = {
    val out = FMat.newOrCheckFMat(nrows, ncols, a, GUID, "toFMat".##)
    cublasGetVector(nrows*ncols, Sizeof.FLOAT, pdata, 1, Pointer.to(out.data), 1)
    cudaStreamSynchronize(Mat.SyncMethod);
    val err = cudaGetLastError
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cublas error in toFMat " + cudaGetErrorString(err))
    }
    out
  }
  
  override def copyTo(a:FMat):FMat = {
    ND.checkDims("copyTo", dims, a.dims);
//  		val a = out.recycle(nrows, ncols, 0)
    cublasGetVector(length, Sizeof.FLOAT, pdata, 1, Pointer.to(a.data), 1)
    cudaStreamSynchronize(Mat.SyncMethod);
    val err = cudaGetLastError;
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU);
    	throw new RuntimeException("Cublas error in copyTo " + cudaGetErrorString(err));
    }
    a;
  }
  
  def copyTo(a:GIMat):GIMat = {
    ND.checkDims("copyTo", dims, a.dims)
    val err = CUMAT.floatToInt(pdata, a.pdata, length);
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU);
    	throw new RuntimeException("error in copyTo " + cudaGetErrorString(err));
    }
    a
  }
  
  def vecAdd(fromi:Int, b:GMat, toi:Int, n:Int):GMat = {
    val bb = b.pdata.withByteOffset(toi*Sizeof.FLOAT);
    CUMAT.applyop(pdata.withByteOffset(fromi*Sizeof.FLOAT), n, 1, bb, n, 1, bb, GMat.BinOp.op_add);
    b;
  }
  
  override def vecAdd(fromi:Int, b:Mat, toi:Int, n:Int):Mat = {
    b match {
      case bb:GMat => vecAdd(fromi, bb, toi, n);
    }
  }
  
  def tileCopy(fromrow:Int, fromcol:Int, to:GMat, torow:Int, tocol:Int, height:Int, width:Int):GMat = {
    val toindx = torow + tocol * to.nrows;
    val fromindx = fromrow + fromcol * nrows;
    cudaMemcpy2D(to.pdata.withByteOffset(toindx * Sizeof.FLOAT), to.nrows*Sizeof.FLOAT, pdata.withByteOffset(fromindx * Sizeof.FLOAT), nrows*Sizeof.FLOAT,
        height*Sizeof.FLOAT, width, cudaMemcpyKind.cudaMemcpyDeviceToDevice);  
    to
  }
  
  override def tileCopy(fromrow:Int, fromcol:Int, to:Mat, torow:Int, tocol:Int, height:Int, width:Int):FMat = {
    tileCopy(fromrow, fromcol, to.asInstanceOf[GMat], torow, tocol, height, width);
  }
  
  def copyFrom(aa:FMat):GMat = {
  	ND.checkDims("GMat copyFrom FMat", dims, aa.dims);
    aa match {
      case in:GMat => cudaMemcpy(pdata, in.pdata, 1L*nrows*ncols*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyDeviceToDevice); 
      case in:FMat => cudaMemcpy(pdata, Pointer.to(in.data), 1L*nrows*ncols*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyHostToDevice);
    }
    cudaStreamSynchronize(Mat.SyncMethod);
    val err = cudaGetLastError;
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU);
    	throw new RuntimeException("Cublas error in copyFrom " + cudaGetErrorString(err));
    }
    this
  }
  
  def copyTo(a:GMat):GMat = {
	  ND.checkDims("GMat copyTo GMat", dims, a.dims); 
    cudaMemcpy(a.pdata, pdata, 1L*length*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyDeviceToDevice);
    cudaStreamSynchronize(Mat.SyncMethod);
    val err = cudaGetLastError;
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU);
    	throw new RuntimeException("Cublas error in copyTo " + cudaGetErrorString(err));
    }
    a
  }
  
  override def copyTo(out:Mat):Mat = {
    out match {
      case a:GMat => copyTo(a)
      case a:GIMat => copyTo(a)
      case a:FMat => copyTo(a)
      case a:TMat => copyTo(a)
    }
  }
  
  override def copy() = {
    val out = GMat.newOrCheckGMat(dims, null, GUID, "GMat.copy".##)
    copyTo(out)
  }

  
  override def recycle(nr:Int, nc:Int, nnz:Int):GMat = {
    if (nrows == nr && nc == ncols) {
      this
    } else if (realsize >= nr*nc) {
      new GMat(nr, nc, pdata, realsize)
    } else {
//      free
      GMat(nr, nc)
    }  
  }
  
  override def free() = {
    if (pdata == null) throw new RuntimeException("attempt to free a free'd matrix");
    cudaFree(pdata)
    pdata = null;
    this
  }
  
  override def finalize = {
//    if (pdata != null) free
  }
  
  override def getdiag():GMat = {
    if (nrows != ncols) throw new RuntimeException("getdiag requires a square matrix, but dims= %d %d" format (nrows, ncols))
    val out = GMat.newOrCheckGMat(nrows, 1, null, GUID, "getdiag".##)
    cudaMemcpy2D(out.pdata, Sizeof.FLOAT, pdata, (nrows+1)*Sizeof.FLOAT, Sizeof.FLOAT, nrows, cudaMemcpyDeviceToDevice)
    cudaStreamSynchronize(Mat.SyncMethod)
    val err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in getdiag " + cudaGetErrorString(err))
    }
    out
  }
  
    
  override def mkdiag():GMat = {
    if (math.min(nrows, ncols) != 1) throw new RuntimeException("mkdiag requires a vector argument, but dims= %d %d" format (nrows, ncols))
    val size = math.max(nrows, ncols)
    val out = GMat.newOrCheckGMat(size, size, null, GUID, "mkdiag".##)
    out.clear
    var err = cudaMemcpy2D(out.pdata, (nrows+1)*Sizeof.FLOAT, pdata, Sizeof.FLOAT, Sizeof.FLOAT, nrows, cudaMemcpyDeviceToDevice)
    cudaStreamSynchronize(Mat.SyncMethod)
    if (err == 0) err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in mkdiag " + cudaGetErrorString(err))
    }
    out
  }
  
  def blockGemm(transa:Int, transb:Int, nr:Int, nc:Int, k:Int, alpha:Float, aoff:Int, lda:Int, astep:Int, 
      b:GMat, boff:Int, ldb:Int, bstep:Int, beta:Float, c:GMat, coff:Int, ldc:Int, cstep:Int, nreps:Int):GMat = {

      Mat.nflops += 2L * nr * nc * k * nreps;
      if (lda > astep || ldb > bstep || ldc > cstep) { 
    	  CUMAT.myCublasSgemmStridedBatched(
    			  getHandle, transa, transb,	
    			  nr, nc, k, 
    			  alpha, 
    			  pdata.withByteOffset(1L * Sizeof.FLOAT * aoff), lda, astep, 
    			  b.pdata.withByteOffset(1L * Sizeof.FLOAT * boff), ldb, bstep, 
    			  beta, 
    			  c.pdata.withByteOffset(1L * Sizeof.FLOAT * coff), ldc, cstep,
    			  nreps);
      } else { 
    	  cublasSgemmStridedBatched(
    			  getHandle, transa, transb,	
    			  nr, nc, k, 
    			  Pointer.to(Array(alpha)),
    			  pdata.withByteOffset(1L * Sizeof.FLOAT * aoff), lda, astep, 
    			  b.pdata.withByteOffset(1L * Sizeof.FLOAT * boff), ldb, bstep, 
    			  Pointer.to(Array(beta)),
    			  c.pdata.withByteOffset(1L * Sizeof.FLOAT * coff), ldc, cstep,
    			  nreps);
      }
      cudaStreamSynchronize(Mat.SyncMethod)

      val err = cudaGetLastError()
      if (err != 0) {
    	  println("device is %d" format SciFunctions.getGPU)
    	  throw new RuntimeException("Cuda error in GMat blockGemm " + cudaGetErrorString(err))
      }
    c;
  }
  
  override def blockGemm(transa:Int, transb:Int, nr:Int, nc:Int, k:Int, alpha:Float, aoff:Int, lda:Int, astep:Int, 
      b:Mat, boff:Int, ldb:Int, bstep:Int, beta:Float, c:Mat, coff:Int, ldc:Int, cstep:Int, nreps:Int):GMat = {
  		blockGemm(transa, transb, nr, nc, k, alpha, aoff, lda, astep, b.asInstanceOf[GMat], boff, ldb, bstep, 
  		    beta, c.asInstanceOf[GMat], coff, ldc, cstep, nreps);
  }
  
  def blockGemm4D(transa:Int, transb:Int, nr:Int, nc:Int, k:Int, alpha:Float, 
		  aoff:Int, lda:Int, astep1:Int, astep2:Int,  
		  b:GMat, boff:Int, ldb:Int, bstep1:Int, bstep2:Int, beta:Float, 
		  c:GMat, coff:Int, ldc:Int, cstep1:Int, cstep2:Int, nreps1:Int, nreps2:Int):GMat = {

    Mat.nflops += 2L * nr * nc * k * nreps1 * nreps2;
    CUMAT.myCublasSgemmStridedBatched4D(
    		getHandle, transa, transb, nr, nc, k, alpha, 
    		pdata.withByteOffset(1L * Sizeof.FLOAT * aoff), lda, astep1, astep2, 
    		b.pdata.withByteOffset(1L * Sizeof.FLOAT * boff), ldb, bstep1, bstep2, beta, 
    		c.pdata.withByteOffset(1L * Sizeof.FLOAT * coff), ldc, cstep1, cstep2,
    		nreps1, nreps2);
    cudaStreamSynchronize(Mat.SyncMethod)

    val err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in GMat blockGemm4D " + cudaGetErrorString(err))
    }
    c;
  }
  
  override def blockGemm4D(transa:Int, transb:Int, nr:Int, nc:Int, k:Int, alpha:Float, aoff:Int, lda:Int, astep1:Int, astep2:Int, 
      b:Mat, boff:Int, ldb:Int, bstep1:Int, bstep2:Int, beta:Float, c:Mat, coff:Int, ldc:Int, cstep1:Int, cstep2:Int, nreps1:Int, nreps2:Int):GMat = {
  		blockGemm4D(transa, transb, nr, nc, k, alpha, aoff, lda, astep1, astep2, b.asInstanceOf[GMat], boff, ldb, bstep1, bstep2, 
  		    beta, c.asInstanceOf[GMat], coff, ldc, cstep1, cstep2, nreps1, nreps2);
  }
  
  override def cumsumByKey(fkeys:FMat, omat:Mat):GMat = {
    val keys = GMat(fkeys);
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumsumKey dimensions mismatch");
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, GUID, keys.GUID, "cumsumKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1 || ncols == 1) {
      val err = CUMAT.cumsumByKeyFF(pdata, keys.pdata, out.pdata, llength);
      if (err != 0) {
    		throw new RuntimeException("CUMAT.cumsumByKey error " + cudaGetErrorString(err))
      }
    } else {
      val tmp = GLMat(nrows, ncols);
      var err = CUMAT.embedmat2d(keys.pdata, tmp.pdata, nrows, ncols, 0);
      if (err == 0) err = CUMAT.cumsumByKeyFL(pdata, tmp.pdata, out.pdata, llength);
      if (err != 0) {
    		throw new RuntimeException("CUMAT.cumsumByKey error " + cudaGetErrorString(err))
      }      
      tmp.free;
    }
    out  
  }
  
  override def cumsumByKey(ikeys:IMat, omat:Mat):GMat = {
    val keys = GIMat(ikeys);
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumsumKey dimensions mismatch");
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, GUID, keys.GUID, "cumsumKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1 || ncols == 1) {
      val err = CUMAT.cumsumByKeyFI(pdata, keys.pdata, out.pdata, llength);
      if (err != 0) {
    		throw new RuntimeException("CUMAT.cumsumByKey error " + cudaGetErrorString(err))
      }
    } else {
    	val tmp = GLMat(nrows, ncols);
      var err = CUMAT.embedmat2d(keys.pdata, tmp.pdata, nrows, ncols, 0);
      if (err == 0) err = CUMAT.cumsumByKeyFL(pdata, tmp.pdata, out.pdata, llength);
      if (err != 0) {
    		throw new RuntimeException("CUMAT.cumsumByKey error " + cudaGetErrorString(err))
      }
      tmp.free;
    }
    out  
  }
  
  override def cumsumByKey(keys:FMat):GMat = cumsumByKey(keys, null);
    
  override def cumsumByKey(keys:IMat):GMat = cumsumByKey(keys, null);
  
  override def cummaxByKey(fkeys:FMat, omat:Mat):GMat = {
    val keys = GMat(fkeys);
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cummaxKey dimensions mismatch");
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, GUID, keys.GUID, "cummaxKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1 || ncols == 1) {
      val err = CUMAT.cummaxByKeyFF(pdata, keys.pdata, out.pdata, llength);
      if (err != 0) {
        throw new RuntimeException("CUMAT.cummaxByKey error " + cudaGetErrorString(err))
      }
    } else {
      val tmp = GLMat(nrows, ncols);
      var err = CUMAT.embedmat2d(keys.pdata, tmp.pdata, nrows, ncols, 0);
      if (err == 0) err = CUMAT.cummaxByKeyFL(pdata, tmp.pdata, out.pdata, llength);
      if (err != 0) {
        throw new RuntimeException("CUMAT.cummaxByKey error " + cudaGetErrorString(err))
      }      
      tmp.free;
    }
    out  
  }
  
  override def cummaxByKey(ikeys:IMat, omat:Mat):GMat = {
    val keys = GIMat(ikeys);
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cummaxKey dimensions mismatch");
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, GUID, keys.GUID, "cummaxKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1 || ncols == 1) {
      val err = CUMAT.cummaxByKeyFI(pdata, keys.pdata, out.pdata, llength);
      if (err != 0) {
        throw new RuntimeException("CUMAT.cummaxByKey error " + cudaGetErrorString(err))
      }
    } else {
      val tmp = GLMat(nrows, ncols);
      var err = CUMAT.embedmat2d(keys.pdata, tmp.pdata, nrows, ncols, 0);
      if (err == 0) err = CUMAT.cummaxByKeyFL(pdata, tmp.pdata, out.pdata, llength);
      if (err != 0) {
        throw new RuntimeException("CUMAT.cummaxByKey error " + cudaGetErrorString(err))
      }
      tmp.free;
    }
    out  
  }
  
  override def cummaxByKey(keys:FMat):GMat = cummaxByKey(keys, null);
    
  override def cummaxByKey(keys:IMat):GMat = cummaxByKey(keys, null);
  
   override def cumminByKey(fkeys:FMat, omat:Mat):GMat = {
     val keys = GMat(fkeys);
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumminKey dimensions mismatch");
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, GUID, keys.GUID, "cumminKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1 || ncols == 1) {
      val err = CUMAT.cumminByKeyFF(pdata, keys.pdata, out.pdata, llength);
      if (err != 0) {
        throw new RuntimeException("CUMAT.cumminByKey error " + cudaGetErrorString(err))
      }
    } else {
      val tmp = GLMat(nrows, ncols);
      var err = CUMAT.embedmat2d(keys.pdata, tmp.pdata, nrows, ncols, 0);
      if (err == 0) err = CUMAT.cumminByKeyFL(pdata, tmp.pdata, out.pdata, llength);
      if (err != 0) {
        throw new RuntimeException("CUMAT.cumminByKey error " + cudaGetErrorString(err))
      }      
      tmp.free;
    }
    out  
  }
  
  override def cumminByKey(ikeys:IMat, omat:Mat):GMat = {
    val keys = GIMat(ikeys);
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumminKey dimensions mismatch");
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, GUID, keys.GUID, "cumminKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1 || ncols == 1) {
      val err = CUMAT.cumminByKeyFI(pdata, keys.pdata, out.pdata, llength);
      if (err != 0) {
        throw new RuntimeException("CUMAT.cumminByKey error " + cudaGetErrorString(err))
      }
    } else {
      val tmp = GLMat(nrows, ncols);
      var err = CUMAT.embedmat2d(keys.pdata, tmp.pdata, nrows, ncols, 0);
      if (err == 0) err = CUMAT.cumminByKeyFL(pdata, tmp.pdata, out.pdata, llength);
      if (err != 0) {
        throw new RuntimeException("CUMAT.cumminByKey error " + cudaGetErrorString(err))
      }
      tmp.free;
    }
    out  
  }
  
  override def cumminByKey(keys:FMat):GMat = cumminByKey(keys, null);
    
  override def cumminByKey(keys:IMat):GMat = cumminByKey(keys, null);

  override def _reverse(omat:Mat):GMat = {
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, GUID,  "reverse".##);
    val err = CUMAT.reverse(pdata, out.pdata, llength);
    if (err != 0) {
    	throw new RuntimeException("CUMAT.reverse error " + cudaGetErrorString(err))
    }
    out
  }
  
  override def reverse:GMat = _reverse(null);
  
  override def reverse(omat:Mat):GMat = _reverse(omat);
  
  def reduce(inds:Array[Int], fctn:(GMat,Int)=>GMat, opname:String):GMat = {
    val alldims = MatFunctions.izeros(_dims.length,1);
    val xinds = new IMat(inds.length, 1, inds);
    val xdims = new IMat(_dims.length, 1, _dims);
    alldims(xinds) = 1;
    if (alldims.data.reduce(_+_) != inds.length) {
      throw new RuntimeException(opname+ " indices arent a legal subset of dims");
    }
    val restinds = MatFunctions.find(alldims == 0);
    if (restinds(0) == 0) {
    	val tmp = transpose((restinds on xinds).data);
    	val tmpF = new GMat(xdims(restinds).data.reduce(_*_), xdims(xinds).data.reduce(_*_), tmp.pdata, length);
    	tmpF.setGUID(ND.hash3(ND.hashInts(inds), GUID, ("reduce"+opname).##));
    	val reduced:GMat = fctn(tmpF, 2);
    	val pdims = xdims(restinds) on MatFunctions.iones(inds.length,1);
    	val out1 = new GMat(pdims.data, reduced.pdata, reduced.length);
    	out1.setGUID(ND.hash3(ND.hashInts(inds), GUID, ("reduce2"+opname).##));
    	out1.transpose(MatFunctions.invperm(restinds on xinds).data)
    } else {
    	val tmp = transpose((xinds on restinds).data);
    	val tmpF = new GMat(xdims(xinds).data.reduce(_*_), xdims(restinds).data.reduce(_*_), tmp.pdata, length);
    	tmpF.setGUID(ND.hash3(ND.hashInts(inds), GUID, ("reduce"+opname).##));
    	val reduced:GMat = fctn(tmpF, 1);
    	val newdims = MatFunctions.iones(inds.length,1) on xdims(restinds);
    	val out1 = new GMat(newdims.data, reduced.pdata, reduced.length);
    	out1.setGUID(ND.hash3(ND.hashInts(inds), GUID, ("reduce2"+opname).##));
    	out1.transpose(MatFunctions.invperm(xinds on restinds).data);
    }
  }
  
  def reduce(inds0:Array[Int], fctn:(GMat,Int)=>GMat, fred:(Pointer, Pointer, Int, Int, Int)=>Int, opname:String):GMat = {
    var i = 1;
    while (i < inds0.length) {
      if (inds0(i-1) >= inds0(i)) {
        throw new RuntimeException("GMat reduce bad index vector");
      }
      i += 1;
    }
    var inmat = this;
    var outmat = this;
    var inds = MatFunctions.irow(inds0);
    var nextinds = getNextInds(inds);
    var restinds = inds;
    while (nextinds.asInstanceOf[AnyRef] != null) {
    	restinds = if (restinds.length > nextinds.length) restinds.colslice(0, restinds.length-nextinds.length) else null;
    	val outdims = inmat.dims.copy;
    	outdims(nextinds) = 1;
    	var n = 1;
      for (i <- nextinds(0) to nextinds(nextinds.length-1)) n *= inmat.dims(i);
      var k = 1;
      for (i <- (nextinds(nextinds.length-1)+1) until inmat.dims.length) k *= inmat.dims(i);
    	if (nextinds(0) == 0) {
    	  val tmpin = inmat.reshapeView(n, k);
    	  val tmpout = fctn(tmpin,1);
    	  outmat = tmpout.reshapeView(outdims);
    	} else {
    		outmat = GMat.newOrCheckGMat(outdims, null, inmat.GUID, ND.hashInts(outdims.data), "GMat_reduce".##);
    		var m = 1;
    		for (i <- 0 until nextinds(0)) m *= inmat.dims(i);
    		fred(inmat.pdata, outmat.pdata, m, n, k);
    		Mat.nflops += inmat.length;
    	}
    	nextinds = getNextInds(restinds);
    	inmat = outmat;
    }
    outmat;
  }
  
  /*
   * Basic compute routines on pairs of GMats
   */
  override def unary_-() = {
    val minusOne = GMat.newOrCheckGMat(1,1,null,-1,"minusOne".##);
    minusOne.set(-1f);
    gOp(minusOne, null, op_mul)
  }
  
  def + (a : GMat) = gOp(a, null, op_add)
  def - (a : GMat) = gOp(a, null, op_sub)
  def *@ (a : GMat) = gOp(a, null, op_mul)
  def ∘  (a : GMat) = gOp(a, null, op_mul)
  def /  (a : GMat) = gOp(a, null, op_div)
  def ^  (a : GMat) = gOp(a, null, op_pow)
  def * (a : GMat) = GMult(a, null)
  def * (a : GSMat) = GSMult(a, null)
  def *^ (a : GMat) = GMultT(a, null)
  def *^ (a : GSMat) = GSMultT(a, null)
  def xT (a : GMat) = GMultT(a, null)
  def xT (a : GSMat) = GSMultT(a, null)
  def ^* (a : GMat) = GTMult(a, null)
  def *+^ (a : GMat) = GMST(a, null)
  def Tx (a : GMat) = GTMult(a, null)
  def kron(a: GMat):GMat = kron(a, null)
  def ⊗  (b : GMat) = kron(b, null)
  def ∙  (a : GMat) = dot(a)
  def ∙→ (a : GMat) = dotr(a)
  
  def > (b : GMat) = gOp(b, null, op_gt)
  def < (b : GMat) = gOp(b, null, op_lt)
  def == (b : GMat) = gOp(b, null, op_eq)
  def === (b : GMat) = gOp(b, null, op_eq)
  def >= (b : GMat) = gOp(b, null, op_ge)
  def <= (b : GMat) = gOp(b, null, op_le)
  def != (b : GMat) = gOp(b, null, op_ne)
  
  def max (b : GMat) = gOp(b, null, op_max)
  def min (b : GMat) = gOp(b, null, op_min)

  override def sum(ind:Int):GMat = reduceOp(null, ind, 0f, op_add);
  override def prod(ind:Int):GMat = reduceOp(null, ind, 1f, op_mul);
  override def maxi(ind:Int):GMat = reduceOp(null, ind, Float.MinValue, op_max);
  override def amin(ind:Int):GMat = reduceOp(null, ind, Float.MaxValue, op_min);
  override def amax(ind:Int):GMat = reduceOp(null, ind, Float.MinValue, op_max);
  override def mini(ind:Int):GMat = reduceOp(null, ind, Float.MaxValue, op_min); 
  override def mean(ind:Int):GMat = SciFunctions._mean(this, ind).asInstanceOf[GMat];
  override def variance(ind:Int):GMat = SciFunctions._variance(this, ind).asInstanceOf[GMat];
  
  override def sum(inds:Array[Int]):FMat = reduce(inds, (a:GMat, dir:Int) => GFunctions.sum(a,dir,null), CUMAT.sumTensor, "sum");
  override def prod(inds:Array[Int]):FMat = reduce(inds, (a:GMat, dir:Int) => GFunctions.prod(a,dir,null), CUMAT.prodTensor, "prod");
  override def maxi(inds:Array[Int]):FMat = reduce(inds, (a:GMat, dir:Int) => GFunctions.maxi(a,dir,null), CUMAT.maxTensor, "maxi")
  override def mini(inds:Array[Int]):FMat = reduce(inds, (a:GMat, dir:Int) => GFunctions.mini(a,dir,null), CUMAT.minTensor, "mini")
  override def amax(inds:Array[Int]):FMat = reduce(inds, (a:GMat, dir:Int) => GFunctions.maxi(a,dir,null), CUMAT.maxTensor, "amax")
  override def amin(inds:Array[Int]):FMat = reduce(inds, (a:GMat, dir:Int) => GFunctions.mini(a,dir,null), CUMAT.minTensor,"amin")
  override def mean(inds:Array[Int]):FMat = reduce(inds, (a:GMat, dir:Int) => SciFunctions.mean(a,dir), "mean")
  override def variance(inds:Array[Int]):FMat = reduce(inds, (a:GMat, dir:Int) => SciFunctions.variance(a,dir), "variance")

  override def sum(inds:IMat):FMat = reduce(inds.data, (a:GMat, dir:Int) => GFunctions.sum(a,dir,null), CUMAT.sumTensor, "sum");
  override def prod(inds:IMat):FMat = reduce(inds.data, (a:GMat, dir:Int) => GFunctions.prod(a,dir,null), CUMAT.prodTensor, "prod");
  override def maxi(inds:IMat):FMat = reduce(inds.data, (a:GMat, dir:Int) => GFunctions.maxi(a,dir,null), CUMAT.maxTensor, "maxi")
  override def mini(inds:IMat):FMat = reduce(inds.data, (a:GMat, dir:Int) => GFunctions.mini(a,dir,null), CUMAT.minTensor, "mini")
  override def amax(inds:IMat):FMat = reduce(inds.data, (a:GMat, dir:Int) => GFunctions.maxi(a,dir,null), CUMAT.maxTensor, "amax")
  override def amin(inds:IMat):FMat = reduce(inds.data, (a:GMat, dir:Int) => GFunctions.mini(a,dir,null), CUMAT.minTensor,"amin")
  override def mean(inds:IMat):FMat = reduce(inds.data, (a:GMat, dir:Int) => SciFunctions.mean(a,dir), "mean")
  override def variance(inds:IMat):FMat = reduce(inds.data, (a:GMat, dir:Int) => SciFunctions.variance(a,dir), "variance")

  override def * (a : FMat) = GMult(GMat(a), null)
  override def * (a : SMat) = GSMult(GSMat(a), null)
  override def *^ (a : FMat) = GMultT(GMat(a), null)
  override def *^ (a : SMat) = GSMultT(GSMat(a), null)
  override def xT (a : FMat) = GMultT(GMat(a), null)
  override def xT (a : SMat) = GSMultT(GSMat(a), null)
  override def ^* (a : FMat) = GTMult(GMat(a), null)
  def *+^ (a : FMat) = GMST(GMat(a), null)
  override def Tx (a : FMat) = GTMult(GMat(a), null)
  override def kron(a: FMat) = kron(GMat(a), null)
  override def ⊗  (a : FMat) = kron(GMat(a), null)
  override def + (a : FMat) = gOp(GMat(a), null, op_add)
  override def - (a : FMat) = gOp(GMat(a), null, op_sub)
  override def *@ (a : FMat) = gOp(GMat(a), null, op_mul)
  override def ∘  (a : FMat) = gOp(GMat(a), null, op_mul)
  override def /  (a : FMat) = gOp(GMat(a), null, op_div)
  override def ^  (a : FMat) = gOp(GMat(a), null, op_pow)
  override def ∙  (a : FMat) = dot(a)
  override def ∙→ (a : FMat) = dotr(a)
  
  override def > (a : FMat) = gOp(GMat(a), null, op_gt)
  override def < (a : FMat) = gOp(GMat(a), null, op_lt)
  override def == (a : FMat) = gOp(GMat(a), null, op_eq)
  override def === (a : FMat) = gOp(GMat(a), null, op_eq)
  override def >= (a : FMat) = gOp(GMat(a), null, op_ge)
  override def <= (a : FMat) = gOp(GMat(a), null, op_le)
  override def != (a : FMat) = gOp(GMat(a), null, op_ne)
  
  override def max (a : FMat) = gOp(GMat(a), null, op_max)
  override def min (a : FMat) = gOp(GMat(a), null, op_min)
  
  override def + (a : Float) = gOp(GMat.elem(a), null, op_add)
  override def - (a : Float) = gOp(GMat.elem(a), null, op_sub)
  override def *@ (a : Float) = gOp(GMat.elem(a), null, op_mul)
  override def * (a : Float) = gOp(GMat.elem(a), null, op_mul)
  override def ∘  (a : Float) = gOp(GMat.elem(a), null, op_mul)
  override def /  (a : Float) = gOp(GMat.elem(a), null, op_div)
  override def ^  (a : Float) = gOp(GMat.elem(a), null, op_pow)
  
  override def < (b : Float) = gOp(GMat.elem(b), null, op_lt);
  override def > (b : Float) = gOp(GMat.elem(b), null, op_gt);
  override def <= (b : Float) = gOp(GMat.elem(b), null, op_le);
  override def >= (b : Float) = gOp(GMat.elem(b), null, op_ge);
  override def == (b : Float) = gOp(GMat.elem(b), null, op_eq);
  override def != (b : Float) = gOp(GMat.elem(b), null, op_ne);
  
  override def max (b : Float) = gOp(GMat.elem(b), null, op_max)
  override def min (b : Float) = gOp(GMat.elem(b), null, op_min)
  
  
  override def + (a : Double) = gOp(GMat(a.toFloat), null, op_add)
  override def - (a : Double) = gOp(GMat(a.toFloat), null, op_sub)
  override def *@ (a : Double) = gOp(GMat(a.toFloat), null, op_mul)
  override def * (a : Double) = gOp(GMat(a.toFloat), null, op_mul)
  override def ∘  (a : Double) = gOp(GMat(a.toFloat), null, op_mul)
  override def /  (a : Double) = gOp(GMat(a.toFloat), null, op_div)
  override def ^  (a : Double) = gOp(GMat(a.toFloat), null, op_pow)
  
  override def < (b : Double) = gOp(GMat(b.toFloat), null, op_lt)
  override def > (b : Double) = gOp(GMat(b.toFloat), null, op_gt)
  override def <= (b : Double) = gOp(GMat(b.toFloat), null, op_le)
  override def >= (b : Double) = gOp(GMat(b.toFloat), null, op_ge)
  override def == (b : Double) = gOp(GMat(b.toFloat), null, op_eq)
  override def != (b : Double) = gOp(GMat(b.toFloat), null, op_ne)

  override def max (b : Double) = gOp(GMat.elem(b), null, op_max)
  override def min (b : Double) = gOp(GMat.elem(b), null, op_min)
  
  
  override def + (a : Int) = gOp(GMat(a.toFloat), null, op_add)
  override def - (a : Int) = gOp(GMat(a.toFloat), null, op_sub)
  override def *@ (a : Int) = gOp(GMat(a.toFloat), null, op_mul)
  override def * (a : Int) = gOp(GMat(a.toFloat), null, op_mul)
  override def ∘  (a : Int) = gOp(GMat(a.toFloat), null, op_mul)
  override def /  (a : Int) = gOp(GMat(a.toFloat), null, op_div)
  override def ^  (a : Int) = gOp(GMat(a.toFloat), null, op_pow)
  
  override def < (b : Int) = gOp(GMat(b.toFloat), null, op_lt)
  override def > (b : Int) = gOp(GMat(b.toFloat), null, op_gt)
  override def <= (b : Int) = gOp(GMat(b.toFloat), null, op_le)
  override def >= (b : Int) = gOp(GMat(b.toFloat), null, op_ge)
  override def == (b : Int) = gOp(GMat(b.toFloat), null, op_eq)
  override def != (b : Int) = gOp(GMat(b.toFloat), null, op_ne)
  
  override def max (b : Int) = gOp(GMat.elem(b), null, op_max)
  override def min (b : Int) = gOp(GMat.elem(b), null, op_min)

  
  override def + (a : Long) = gOp(GMat(a.toFloat), null, op_add)
  override def - (a : Long) = gOp(GMat(a.toFloat), null, op_sub)
  override def *@ (a : Long) = gOp(GMat(a.toFloat), null, op_mul)
  override def * (a : Long) = gOp(GMat(a.toFloat), null, op_mul)
  override def ∘  (a : Long) = gOp(GMat(a.toFloat), null, op_mul)
  override def /  (a : Long) = gOp(GMat(a.toFloat), null, op_div)
  override def ^  (a : Long) = gOp(GMat(a.toFloat), null, op_pow)
  
  override def < (b : Long) = gOp(GMat(b.toFloat), null, op_lt)
  override def > (b : Long) = gOp(GMat(b.toFloat), null, op_gt)
  override def <= (b : Long) = gOp(GMat(b.toFloat), null, op_le)
  override def >= (b : Long) = gOp(GMat(b.toFloat), null, op_ge)
  override def == (b : Long) = gOp(GMat(b.toFloat), null, op_eq)
  override def != (b : Long) = gOp(GMat(b.toFloat), null, op_ne)
  
  override def max (b : Long) = gOp(GMat.elem(b), null, op_max)
  override def min (b : Long) = gOp(GMat.elem(b), null, op_min)

  
  def on(a : GMat) = vertcat(a, null)
  def \ (a : GMat) = horzcat(a, null)
   
 /*
  * Specialize to IMats to help the type system. 
  */
  override def *   (b : IMat) = Mop_Times.op(this, b, null) 
  override def *^  (b : IMat) = Mop_TimesT.op(this, b, null)
  override def xT  (b : IMat) = Mop_TimesT.op(this, b, null)
  override def Tx  (b : IMat) = Mop_TTimes.op(this, b, null)
  override def ^*  (b : IMat) = Mop_TTimes.op(this, b, null)
  override def +   (b : IMat) = Mop_Plus.op(this, b, null)
  override def -   (b : IMat) = Mop_Minus.op(this, b, null)
  override def *@  (b : IMat) = Mop_ETimes.op(this, b, null)
  override def ∘   (b : IMat) = Mop_ETimes.op(this, b, null)
  override def /<  (b : IMat) = Mop_Div.op(this, b, null)
  override def \\  (b : IMat) = Mop_RSolve.op(this, b, null)
  override def ◁   (b : IMat) = Mop_Div.op(this, b, null)
  override def ▷   (b : IMat) = Mop_RSolve.op(this, b, null)
  override def /   (b : IMat) = Mop_EDiv.op(this, b, null)  
  override def ^   (b : IMat) = Mop_Pow.op(this, b, null) 
  override def ∙   (b : IMat) = Mop_Dot.op(this, b, null)
  override def ∙→  (b : IMat) = Mop_Dotr.op(this, b, null)
  override def dot (b : IMat) = Mop_Dot.op(this, b, null)
  override def dotr(b : IMat) = Mop_Dotr.op(this, b, null)
  override def \   (b : IMat) = Mop_HCat.op(this, b, null)
  override def on  (b : IMat) = Mop_VCat.op(this, b, null)

  override def >   (b : IMat) = Mop_GT.op(this, b, null)
  override def <   (b : IMat) = Mop_LT.op(this, b, null)
  override def ==  (b : IMat) = Mop_EQ.op(this, b, null)
  override def === (b : IMat) = Mop_EQ.op(this, b, null)
  override def >=  (b : IMat) = Mop_GE.op(this, b, null)
  override def <=  (b : IMat) = Mop_LE.op(this, b, null)
  override def !=  (b : IMat) = Mop_NE.op(this, b, null)
   
 /*
  * Specialize to DMats to help the type system. 
  */ 
  override def *   (b : DMat) = Mop_Times.op(this, b, null) 
  override def *^  (b : DMat) = Mop_TimesT.op(this, b, null)
  override def xT  (b : DMat) = Mop_TimesT.op(this, b, null)
  override def Tx  (b : DMat) = Mop_TTimes.op(this, b, null)
  override def ^*  (b : DMat) = Mop_TTimes.op(this, b, null)
  override def +   (b : DMat) = Mop_Plus.op(this, b, null)
  override def -   (b : DMat) = Mop_Minus.op(this, b, null)
  override def *@  (b : DMat) = Mop_ETimes.op(this, b, null)
  override def ∘   (b : DMat) = Mop_ETimes.op(this, b, null)
  override def /<  (b : DMat) = Mop_Div.op(this, b, null)
  override def \\  (b : DMat) = Mop_RSolve.op(this, b, null)
  override def ◁   (b : DMat) = Mop_Div.op(this, b, null)
  override def ▷   (b : DMat) = Mop_RSolve.op(this, b, null)
  override def /   (b : DMat) = Mop_EDiv.op(this, b, null)  
  override def ^   (b : DMat) = Mop_Pow.op(this, b, null) 
  override def ∙   (b : DMat) = Mop_Dot.op(this, b, null)
  override def ∙→  (b : DMat) = Mop_Dotr.op(this, b, null)
  override def dot (b : DMat) = Mop_Dot.op(this, b, null)
  override def dotr(b : DMat) = Mop_Dotr.op(this, b, null)
  override def \   (b : DMat) = Mop_HCat.op(this, b, null)
  override def on  (b : DMat) = Mop_VCat.op(this, b, null)
  
  override def >   (b : DMat) = Mop_GT.op(this, b, null)
  override def <   (b : DMat) = Mop_LT.op(this, b, null)
  override def ==  (b : DMat) = Mop_EQ.op(this, b, null)
  override def === (b : DMat) = Mop_EQ.op(this, b, null)
  override def >=  (b : DMat) = Mop_GE.op(this, b, null)
  override def <=  (b : DMat) = Mop_LE.op(this, b, null)
  override def !=  (b : DMat) = Mop_NE.op(this, b, null)
 
  
 /*
  * Operators whose second arg is generic. 
  */ 
  override def *  (b : Mat) = Mop_Times.op(this, b, null)
  override def *^ (b : Mat) = Mop_TimesT.op(this, b, null)
  override def xT (b : Mat) = Mop_TimesT.op(this, b, null)
  override def Tx (b : Mat) = Mop_TTimes.op(this, b, null)
  override def ^* (b : Mat) = Mop_TTimes.op(this, b, null)
  override def +  (b : Mat) = Mop_Plus.op(this, b, null)
  override def -  (b : Mat) = Mop_Minus.op(this, b, null)
  override def *@ (b : Mat) = Mop_ETimes.op(this, b, null)
  override def ∘  (b : Mat) = Mop_ETimes.op(this, b, null)
  override def /  (b : Mat) = Mop_EDiv.op(this, b, null)
  override def /< (b : Mat) = Mop_Div.op(this, b, null)
  override def \\ (b : Mat) = Mop_RSolve.op(this, b, null)
  override def ◁  (b : Mat) = Mop_Div.op(this, b, null)
  override def ▷  (b : Mat) = Mop_RSolve.op(this, b, null)
  override def ^  (b : Mat) = Mop_Pow.op(this, b, null) 
  override def ∙  (b : Mat) = Mop_Dot.op(this, b, null)
  override def ∙→ (b : Mat) = Mop_Dotr.op(this, b, null)
  override def dot  (b : Mat) = Mop_Dot.op(this, b, null)
  override def dotr (b : Mat) = Mop_Dotr.op(this, b, null)
  override def \  (b : Mat) = Mop_HCat.op(this, b, null)
  override def on (b : Mat) = Mop_VCat.op(this, b, null)
  
  override def >   (b : Mat) = Mop_GT.op(this, b, null)
  override def <   (b : Mat) = Mop_LT.op(this, b, null)
  override def >=  (b : Mat) = Mop_GE.op(this, b, null)
  override def <=  (b : Mat) = Mop_LE.op(this, b, null)
  override def ==  (b : Mat) = Mop_EQ.op(this, b, null)
  override def === (b : Mat) = Mop_EQ.op(this, b, null) 
  override def !=  (b : Mat) = Mop_NE.op(this, b, null)
 
  /*
   * Tilde operator
   */
  def ~ (b: GMat) = new GPair(this, b)
  def ~ (b: GSMat) = new GSPair(this, b)
  override def ~ (b: FMat):FPair = new GPair(this, GMat(b));
  
  override def ~ (b: Mat):Pair = 
    b match {
    case t:TMat => new GTPair(this,t)
    case bb:GMat => new GPair(this, bb)
    case bb:GSMat => new GSPair(this, bb)
  }
 
  /*
   * @@ operator for DDS
   */  
  def @@ (b : GSMat) = new GDSPair(this, b)
  def ^* (b : GDSPair) = MatFunctions.DDS(this, b.left, b.right, null)
  def Tx (b : GDSPair) = MatFunctions.DDS(this, b.left, b.right, null)
  override def ^* (b0 : DSPair) = {val b = b0.asInstanceOf[GDSPair]; MatFunctions.DDS(this, b.left, b.right, null)}
  override def Tx (b0 : DSPair) = {val b = b0.asInstanceOf[GDSPair]; MatFunctions.DDS(this, b.left, b.right, null)}

}

class GTPair(val omat:GMat,val mat:TMat) extends Pair(omat, mat) {
    override def * (a:Mat) = a match {
        case g:GMat => mat.tMult(g,omat)
        case g:GSMat => mat.tMult(g,omat)
    }
}

/*
 * Result of a@@b for DDS
 */
class GDSPair(val left:GMat, val right:GSMat) extends DSPair {}

/*
 * GPair is the result of a~b
 */
class GPair(omat:Mat, override val mat:GMat) extends FPair(omat, mat) {
	import GMat.BinOp._
	
	override def t = {
    val out = GMat.newOrCheckGMat(mat.ncols, mat.nrows, omat, mat.GUID, "pt".##)
    CUMAT.transpose(mat.pdata, mat.nrows, out.pdata, mat.ncols, mat.nrows, mat.ncols)
    out
  }
/*	
  def *  (a : GMat) = mat.GMult(a, omat)
  def *  (a : GSMat) = mat.GSMult(a, omat)
  def *^ (a : GMat) = mat.GMultT(a, omat)
  def *^ (a : GSMat) = mat.GSMultT(a, omat)
  def xT (a : GMat) = mat.GMultT(a, omat)
  def xT (a : GSMat) = mat.GSMultT(a, omat)
  def ^* (a : GMat) = mat.GTMult(a, omat)
  def *+^ (a : GMat) = mat.GMST(a, omat)
  def Tx (a : GMat) = mat.GTMult(a, omat)
  def kron(a: GMat):GMat = mat.kron(a, omat)
  def ⊗  (b : GMat) = mat.kron(b, omat)
	def +  (a : GMat) = mat.gOp(a, omat, op_add)
	def -  (a : GMat) = mat.gOp(a, omat, op_sub)
	def *@ (a : GMat) = mat.gOp(a, omat, op_mul)
	def ∘  (a : GMat) = mat.gOp(a, omat, op_mul)
	def /  (a : GMat) = mat.gOp(a, omat, op_div)
	def ^  (a : GMat) = mat.gOp(a, omat, op_pow)
	def >  (a : GMat) = mat.gOp(a, omat, op_gt)
	def <  (a : GMat) = mat.gOp(a, omat, op_lt)
	def == (a : GMat) = mat.gOp(a, omat, op_eq)
	def === (a : GMat) = mat.gOp(a, omat, op_eq)
	def >= (a : GMat) = mat.gOp(a, omat, op_ge)
	def <= (a : GMat) = mat.gOp(a, omat, op_le)
	def != (a : GMat) = mat.gOp(a, omat, op_ne)
	
  def max (a : GMat) = mat.gOp(a, omat, op_max)
	def min (a : GMat) = mat.gOp(a, omat, op_min)
	
	def dot (a :GMat) = mat.dot(a, omat) 
	def dotr (a :GMat) = mat.dotr(a, omat) 
	def ∙ (a :GMat) = mat.dot(a, omat)
	def ∙→ (a :GMat) = mat.dotr(a, omat)
	def on(a : GMat) = mat.vertcat(a, omat)
	def \ (a : GMat) = mat.horzcat(a, omat)
	*/
  
	def checkOne(a:Seq[Int], name:String):Int = {
    if (a.length > 1) throw new RuntimeException("GMat %s only takes one argument" format name);
    a(0);
  }
	
	
  override def *  (a : FMat) = mat.GMult(GMat(a), omat);
  override def *  (a : SMat) = mat.GSMult(GSMat(a), omat);
  override def *^ (a : FMat) = mat.GMultT(GMat(a), omat)
  override def *^ (a : SMat) = mat.GSMultT(GSMat(a), omat)
  override def xT (a : FMat) = mat.GMultT(GMat(a), omat)
  override def xT (a : SMat) = mat.GSMultT(GSMat(a), omat)
  override def ^* (a : FMat) = mat.GTMult(GMat(a), omat)
  def *+^ (a : FMat) = mat.GMST(GMat(a), omat)
  override def Tx (a : FMat) = mat.GTMult(GMat(a), omat)
  def kron(a: FMat):FMat = mat.kron(GMat(a), omat)
  override def ⊗  (b : FMat) = mat.kron(b, omat)
	override def +  (a : FMat) = mat.gOp(GMat(a), omat, op_add)
	override def -  (a : FMat) = mat.gOp(GMat(a), omat, op_sub)
	override def *@ (a : FMat) = mat.gOp(GMat(a), omat, op_mul)
	override def ∘  (a : FMat) = mat.gOp(GMat(a), omat, op_mul)
	override def /  (a : FMat) = mat.gOp(GMat(a), omat, op_div)
	override def ^  (a : FMat) = mat.gOp(GMat(a), omat, op_pow)
	override def >  (a : FMat) = mat.gOp(GMat(a), omat, op_gt)
	override def <  (a : FMat) = mat.gOp(GMat(a), omat, op_lt)
	override def == (a : FMat) = mat.gOp(GMat(a), omat, op_eq)
	override def === (a : FMat) = mat.gOp(GMat(a), omat, op_eq)
	override def >= (a : FMat) = mat.gOp(GMat(a), omat, op_ge)
	override def <= (a : FMat) = mat.gOp(GMat(a), omat, op_le)
	override def != (a : FMat) = mat.gOp(GMat(a), omat, op_ne)
	
  override def max (a : FMat) = mat.gOp(GMat(a), omat, op_max)
	override def min (a : FMat) = mat.gOp(GMat(a), omat, op_min)
	
	override def dot (a :FMat) = mat.dot(GMat(a), omat) 
	override def dotr (a :FMat) = mat.dotr(GMat(a), omat) 
	override def ∙ (a :FMat) = mat.dot(GMat(a), omat)
	override def ∙→ (a :FMat) = mat.dotr(GMat(a), omat)
	def on(a : FMat) = mat.vertcat(GMat(a), omat)
	def \ (a : FMat) = mat.horzcat(GMat(a), omat)
  
  override def * (b : Float) = mat.gOp(GMat(b), omat, op_mul)
  override def *@ (b : Float) = mat.gOp(GMat(b), omat, op_mul)
  override def ∘ (b : Float) = mat.gOp(GMat(b), omat, op_mul)
  override def + (b : Float) = mat.gOp(GMat(b), omat, op_add)
  override def - (b : Float) = mat.gOp(GMat(b), omat, op_sub)
  override def / (b : Float) = mat.gOp(GMat(b), omat, op_div)
  override def ^ (b : Float) = mat.gOp(GMat(b), omat, op_pow)
  override def >  (b : Float) = mat.gOp(GMat(b), omat, op_gt)
	override def <  (b : Float) = mat.gOp(GMat(b), omat, op_lt)
  override def == (b : Float) = mat.gOp(GMat(b), omat, op_eq)
  override def != (b : Float) = mat.gOp(GMat(b), omat, op_ne)
  override def >= (b : Float) = mat.gOp(GMat(b), omat, op_ge)
	override def <= (b : Float) = mat.gOp(GMat(b), omat, op_le)
	override def max (b : Float) = mat.gOp(GMat(b), omat, op_max)
	override def min (b : Float) = mat.gOp(GMat(b), omat, op_min)
  
  override def * (b : Double) = mat.gOp(GMat(b.toFloat), omat, op_mul)
  override def *@ (b : Double) = mat.gOp(GMat(b.toFloat), omat, op_mul)
  override def ∘ (b : Double) = mat.gOp(GMat(b.toFloat), omat, op_mul)
  override def + (b : Double) = mat.gOp(GMat(b.toFloat), omat, op_add)
  override def - (b : Double) = mat.gOp(GMat(b.toFloat), omat, op_sub)
  override def / (b : Double) = mat.gOp(GMat(b.toFloat), omat, op_div)
  override def ^ (b : Double) = mat.gOp(GMat(b.toFloat), omat, op_pow)
  override def >  (b : Double) = mat.gOp(GMat(b.toFloat), omat, op_gt)
  override def <  (b : Double) = mat.gOp(GMat(b.toFloat), omat, op_lt)
  override def == (b : Double) = mat.gOp(GMat(b.toFloat), omat, op_eq)
  override def != (b : Double) = mat.gOp(GMat(b.toFloat), omat, op_ne)
  override def >= (b : Double) = mat.gOp(GMat(b.toFloat), omat, op_ge)
  override def <= (b : Double) = mat.gOp(GMat(b.toFloat), omat, op_le)
 override def max (b : Double) = mat.gOp(GMat(b), omat, op_max)
	override def min (b : Double) = mat.gOp(GMat(b), omat, op_min)
	
  
	override def * (b : Int) = mat.gOp(GMat(b.toFloat), omat, op_mul)
	override def *@ (b : Int) = mat.gOp(GMat(b.toFloat), omat, op_mul)
  override def ∘ (b : Int) = mat.gOp(GMat(b.toFloat), omat, op_mul)
  override def + (b : Int) = mat.gOp(GMat(b.toFloat), omat, op_add)
  override def - (b : Int) = mat.gOp(GMat(b.toFloat), omat, op_sub)
  override def / (b : Int) = mat.gOp(GMat(b.toFloat), omat, op_div)
  override def ^ (b : Int) = mat.gOp(GMat(b.toFloat), omat, op_pow)
  override def >  (b : Int) = mat.gOp(GMat(b.toFloat), omat, op_gt)
	override def <  (b : Int) = mat.gOp(GMat(b.toFloat), omat, op_lt)
  override def == (b : Int) = mat.gOp(GMat(b.toFloat), omat, op_eq)
  override def != (b : Int) = mat.gOp(GMat(b.toFloat), omat, op_ne)
  override def >= (b : Int) = mat.gOp(GMat(b.toFloat), omat, op_ge)
	override def <= (b : Int) = mat.gOp(GMat(b.toFloat), omat, op_le)
	override def max (b : Int) = mat.gOp(GMat(b.toFloat), omat, op_max)
	override def min (b : Int) = mat.gOp(GMat(b.toFloat), omat, op_min)
  
  
  override def * (b : Long) = mat.gOp(GMat(b.toFloat), omat, op_mul)
  override def *@ (b : Long) = mat.gOp(GMat(b.toFloat), omat, op_mul)
  override def ∘ (b : Long) = mat.gOp(GMat(b.toFloat), omat, op_mul)
  override def + (b : Long) = mat.gOp(GMat(b.toFloat), omat, op_add)
  override def - (b : Long) = mat.gOp(GMat(b.toFloat), omat, op_sub)
  override def / (b : Long) = mat.gOp(GMat(b.toFloat), omat, op_div)
  override def ^ (b : Long) = mat.gOp(GMat(b.toFloat), omat, op_pow)
  override def >  (b : Long) = mat.gOp(GMat(b.toFloat), omat, op_gt)
  override def <  (b : Long) = mat.gOp(GMat(b.toFloat), omat, op_lt)
  override def == (b : Long) = mat.gOp(GMat(b.toFloat), omat, op_eq)
  override def != (b : Long) = mat.gOp(GMat(b.toFloat), omat, op_ne)
  override def >= (b : Long) = mat.gOp(GMat(b.toFloat), omat, op_ge)
  override def <= (b : Long) = mat.gOp(GMat(b.toFloat), omat, op_le)
	override def max (b : Long) = mat.gOp(GMat(b.toFloat), omat, op_max)
	override def min (b : Long) = mat.gOp(GMat(b.toFloat), omat, op_min)


  def ^* (b : GDSPair) = MatFunctions.DDS(mat, b.left, b.right, omat)
  def Tx (b : GDSPair) = MatFunctions.DDS(mat, b.left, b.right, omat)
  /*
   * Specialize to IMat
   */
  override def *   (b : IMat) = Mop_Times.op(mat, b, omat) 
  override def *^  (b : IMat) = Mop_TimesT.op(mat, b, omat)
  override def xT  (b : IMat) = Mop_TimesT.op(mat, b, omat)
  override def Tx  (b : IMat) = Mop_TTimes.op(mat, b, omat)
  override def ^*  (b : IMat) = Mop_TTimes.op(mat, b, omat)
  override def +   (b : IMat) = Mop_Plus.op(mat, b, omat)
  override def -   (b : IMat) = Mop_Minus.op(mat, b, omat)
  override def *@  (b : IMat) = Mop_ETimes.op(mat, b, omat)
  override def ∘   (b : IMat) = Mop_ETimes.op(mat, b, omat)
  override def /   (b : IMat) = Mop_EDiv.op(mat, b, omat)  
  override def ^   (b : IMat) = Mop_Pow.op(mat, b, omat) 
  override def ∙   (b : IMat) = Mop_Dot.op(mat, b, omat)
  override def ∙→  (b : IMat) = Mop_Dotr.op(mat, b, omat)
  override def dot (b : IMat) = Mop_Dot.op(mat, b, omat)
  override def dotr(b : IMat) = Mop_Dotr.op(mat, b, omat)
  override def \   (b : IMat) = Mop_HCat.op(mat, b, omat)
  override def on  (b : IMat) = Mop_VCat.op(mat, b, omat)

  override def >   (b : IMat) = Mop_GT.op(mat, b, omat)
  override def <   (b : IMat) = Mop_LT.op(mat, b, omat)
  override def ==  (b : IMat) = Mop_EQ.op(mat, b, omat)
  override def === (b : IMat) = Mop_EQ.op(mat, b, omat)
  override def >=  (b : IMat) = Mop_GE.op(mat, b, omat)
  override def <=  (b : IMat) = Mop_LE.op(mat, b, omat)
  override def !=  (b : IMat) = Mop_NE.op(mat, b, omat)
  
  /*
   * Specialize to DMat
   */
  override def *   (b : DMat) = Mop_Times.op(mat, b, omat) 
  override def *^  (b : DMat) = Mop_TimesT.op(mat, b, omat)
  override def xT  (b : DMat) = Mop_TimesT.op(mat, b, omat)
  override def Tx  (b : DMat) = Mop_TTimes.op(mat, b, omat)
  override def ^*  (b : DMat) = Mop_TTimes.op(mat, b, omat)
  override def +   (b : DMat) = Mop_Plus.op(mat, b, omat)
  override def -   (b : DMat) = Mop_Minus.op(mat, b, omat)
  override def *@  (b : DMat) = Mop_ETimes.op(mat, b, omat)
  override def ∘   (b : DMat) = Mop_ETimes.op(mat, b, omat)
  override def /   (b : DMat) = Mop_EDiv.op(mat, b, omat)  
  override def ^   (b : DMat) = Mop_Pow.op(mat, b, omat) 
  override def ∙   (b : DMat) = Mop_Dot.op(mat, b, omat)
  override def ∙→  (b : DMat) = Mop_Dotr.op(mat, b, omat)
  override def dot (b : DMat) = Mop_Dot.op(mat, b, omat)
  override def dotr(b : DMat) = Mop_Dotr.op(mat, b, omat)
  override def \   (b : DMat) = Mop_HCat.op(mat, b, omat)
  override def on  (b : DMat) = Mop_VCat.op(mat, b, omat)

  override def >   (b : DMat) = Mop_GT.op(mat, b, omat)
  override def <   (b : DMat) = Mop_LT.op(mat, b, omat)
  override def ==  (b : DMat) = Mop_EQ.op(mat, b, omat)
  override def === (b : DMat) = Mop_EQ.op(mat, b, omat)
  override def >=  (b : DMat) = Mop_GE.op(mat, b, omat)
  override def <=  (b : DMat) = Mop_LE.op(mat, b, omat)
  override def !=  (b : DMat) = Mop_NE.op(mat, b, omat)
  
  /*
   * Generics
   */
  override def *  (b : Mat):Mat = Mop_Times.op(mat, b, omat)
  override def xT (b : Mat):Mat = Mop_TimesT.op(mat, b, omat)
  override def *^ (b : Mat):Mat = Mop_TimesT.op(mat, b, omat)
  override def Tx (b : Mat):Mat = Mop_TTimes.op(mat, b, omat)
  override def ^* (b : Mat):Mat = Mop_TTimes.op(mat, b, omat)
  override def +  (b : Mat):Mat = Mop_Plus.op(mat, b, omat)
  override def -  (b : Mat):Mat = Mop_Minus.op(mat, b, omat)
  override def *@ (b : Mat):Mat = Mop_ETimes.op(mat, b, omat)
  override def ∘  (b : Mat):Mat = Mop_ETimes.op(mat, b, omat)
  override def /  (b : Mat):Mat = Mop_EDiv.op(mat, b, omat)
  override def ^  (b : Mat):Mat = Mop_Pow.op(mat, b, omat) 
  override def /< (b : Mat):Mat = Mop_Div.op(mat, b, omat)
  override def \\ (b : Mat):Mat = Mop_RSolve.op(mat, b, omat)
  override def ◁  (b : Mat):Mat = Mop_Div.op(mat, b, omat)
  override def ▷  (b : Mat):Mat = Mop_RSolve.op(mat, b, omat)
  override def ∙   (b : Mat) = Mop_Dot.op(mat, b, omat)
  override def ∙→  (b : Mat) = Mop_Dotr.op(mat, b, omat)
  override def dot (b : Mat) = Mop_Dot.op(mat, b, omat)
  override def dotr(b : Mat) = Mop_Dotr.op(mat, b, omat)
  override def \  (b : Mat):Mat = Mop_HCat.op(mat, b, omat)
  override def on (b : Mat):Mat = Mop_VCat.op(mat, b, omat)
  
  override def >   (b : Mat):Mat = Mop_GT.op(mat, b, omat)
  override def <   (b : Mat):Mat = Mop_LT.op(mat, b, omat)
  override def >=  (b : Mat):Mat = Mop_GE.op(mat, b, omat)
  override def <=  (b : Mat):Mat = Mop_LE.op(mat, b, omat)
  override def ==  (b : Mat):Mat = Mop_EQ.op(mat, b, omat)
  override def === (b : Mat):Mat = Mop_EQ.op(mat, b, omat) 
  override def !=  (b : Mat):Mat = Mop_NE.op(mat, b, omat)
}


object GMat {
  
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
	  val op_ifpos=14
  }  
  
  object TransF {
    val abs=0
    val exp=1 
    val expm1=2
    val sqrt=3
    val ln=4
    val log10=5 
    val log1p=6
    val cos=7
    val sin=8
    val tan=9
    val cosh=10 
    val sinh=11 
    val tanh=12 
    val acos=13 
    val asin=14
    val atan=15 
    val acosh=16 
    val asinh=17 
    val atanh=18
    val erf=19
    val erfinv=20 
    val erfc=21
    val erfcinv=22 
    val gammaln=23
    val gamma=24
    val ceil=25
    val floor=26
    val round=27
    val trunc=28
    val sign=29
    val exppsi=34
    val normcdf=35
    val normcdfinv=36
    val logistic=37    
  }
  
  var multYieldSize = 1e6f;
  
  object TransF2 {
    val atan2=0
    val pow=1 
  }  
  
   val nullPointer = new Pointer
  
  def zeros(nr:Int, nc:Int) = {
    val out = GMat(nr, nc);
    out.clear;
    out;
  }
  
  def zeros(dims:IMat) = {
    val out = make(dims);
    out.clear;
    out
  }
  
  def ones(nr:Int, nc:Int) = {
    val out = GMat(nr, nc);
    out.set(1f);
    out
  }
  
  def ones(dims:IMat) = {
    val out = make(dims);
    out.set(1f);
    out
  }
  
  val pONE = Pointer.to(Array(1f));
  
  val pZERO = Pointer.to(Array(0f));
 
  
  def apply(nr:Int, nc:Int):GMat = {
    val retv = new GMat(nr, nc, new Pointer(), 1L*nr*nc)  
    if (Mat.debugMem) {
      println("GMat %d %d, %d %f" format (nr, nc, SciFunctions.getGPU, SciFunctions.GPUmem._1))
      if (nr*nc > Mat.debugMemThreshold) throw new RuntimeException("GMat alloc too large");
    }
    var err = if (1L*nr*nc*Sizeof.FLOAT > Mat.hostAllocSize) {
      cudaMallocHost(retv.pdata, 1L*nr*nc*Sizeof.FLOAT);
    } else {
      cudaMalloc(retv.pdata, 1L*nr*nc*Sizeof.FLOAT);
    }
    cudaStreamSynchronize(Mat.SyncMethod);
    if (err == 0) err = cudaGetLastError();
    if (err != 0) throw new RuntimeException("CUDA alloc failed " + cudaGetErrorString(err));
    retv        
  }   
  
  def make(dims:Array[Int]):GMat = {
  	var err = cudaGetLastError();
  	if (err != 0) throw new RuntimeException("Weird previous error in GMat.make " + cudaGetErrorString(err));
	  val len = dims.reduce(_*_);
    val retv = new GMat(dims, new Pointer, len);
    if (Mat.debugMem) {
      println("GMat %d, %d %f" format (len, SciFunctions.getGPU, SciFunctions.GPUmem._1))
      if (len > Mat.debugMemThreshold) throw new RuntimeException("GMat alloc too large");
    }
    err = if (1L*len*Sizeof.FLOAT > Mat.hostAllocSize) {
      cudaMallocHost(retv.pdata, 1L*len*Sizeof.FLOAT);
    } else {
      cudaMalloc(retv.pdata, 1L*len*Sizeof.FLOAT);
    }
    cudaStreamSynchronize(Mat.SyncMethod);
    if (err == 0) err = cudaGetLastError();
    if (err != 0) throw new RuntimeException("CUDA alloc failed " + cudaGetErrorString(err));
    retv       
  }
  
  def make(dims:IMat):GMat = make(dims.data);
  
  def apply(a:FMat):GMat = {
    a match {
      case g:GMat => g;
      case _ => {
    	  val rsize = a.nrows*a.ncols
    			  val retv = GMat.newOrCheckGMat(a.dims, null, a.GUID, "GMat_FMat".##)
    			  cudaMemcpy(retv.pdata, Pointer.to(a.data), 1L*rsize*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyHostToDevice)
    			  cudaStreamSynchronize(Mat.SyncMethod);
    			  val err = cudaGetLastError()
    			  if (err != 0) {
    				  println("device is %d" format SciFunctions.getGPU)
    				  throw new RuntimeException("CUDA error in GMat() " + cudaGetErrorString(err))
    			  }
    	  retv
      }
    }
  }

  def apply(a:DMat):GMat = {
    a match {
      case g:GDMat_ => {
    	val rsize = a.length
    	val retv = GMat.newOrCheckGMat(a.dims, null, a.GUID, "GMat_FMat".##)
        var err = CUMAT.doubleToFloat(g.pdata, retv.pdata, rsize)
    	cudaStreamSynchronize(Mat.SyncMethod);
    	if (err == 0) err = cudaGetLastError()
    	if (err != 0) {
    	  println("device is %d" format SciFunctions.getGPU)
    	  throw new RuntimeException("CUDA error in GMat() " + cudaGetErrorString(err))
    	}
    	retv
      }
      case _ => GMat(FMat(a))
    }
  }

  def apply(a:GIMat):GMat = {
 
    val rsize = a.nrows*a.ncols
    val retv = GMat.newOrCheckGMat(a.dims, null, a.GUID, "GMat_GIMat".##)
    var err = CUMAT.intToFloat(a.pdata, retv.pdata, a.length)
    cudaStreamSynchronize(Mat.SyncMethod);
    if (err == 0) err = cudaGetLastError()
    if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("GMat(GIMat) error " + cudaGetErrorString(err))
    }
    retv
  }
  
  def apply(a:Mat):GMat = a match {
    case aa:GMat => aa
    case aa:GIMat => GMat(aa)
    case aa:FMat => GMat(aa)
    case aa:DMat => GMat(aa)
    case aa:IMat => GMat(FMat(aa))
  }
  
  def apply(a:Float):GMat = {
    val out = GMat.newOrCheckGMat(1, 1, null, a.##, "GMat_Float".##)
    out.set(a)
    out
  }
  
  def apply(a:Double):GMat = {
    val out = GMat.newOrCheckGMat(1, 1, null, a.##, "GMat_Float".##)
    out.set(a.toFloat)
    out
  }

  def elem(a:Float):GMat = {
    val out = GMat.newOrCheckGMat(1, 1, null, a.##, "Gelem".##);
    out.set(a)
    out
  }
  
  def elem(a:Double):GMat = {
    val out = GMat.newOrCheckGMat(1, 1, null, a.##, "Gelem".##);
    out.set(a.toFloat)
    out
  }
  
  def toFMat(a:GMat):FMat = a.toFMat(null)  
  
  def fromFMat(a:FMat, b:GMat):GMat = {
    val bb = GMat.newOrCheckGMat(a.dims, b, a.GUID, "GMat_fromFMat".##)
    cudaMemcpy(bb.pdata, Pointer.to(a.data), a.length*1L*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyHostToDevice)
    cudaStreamSynchronize(Mat.SyncMethod);
    var err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("CUDA error in fromFMat " + cudaGetErrorString(err))
    }
    bb
  }

 
  def GPUtoGPUarraycopy(a:Pointer, aoffset:Int,  b:Pointer, boffset:Int, len:Int, msg:String ) = {
	  cudaMemcpy(b.withByteOffset(1L*boffset*Sizeof.FLOAT), a.withByteOffset(1L*aoffset*Sizeof.FLOAT), 1L*len*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyDeviceToDevice);
    cudaStreamSynchronize(Mat.SyncMethod);
    val err = cudaGetLastError;
    if (err != 0) throw new RuntimeException(msg +" error in memcpy "+ cudaGetErrorString(err));
  }
  
  def GPUtoCPUarraycopy(a:Pointer, aoffset:Int,  b:Array[Float], boffset:Int, len:Int, msg:String ) = {
    cudaMemcpy(Pointer.to(b).withByteOffset(1L*boffset*Sizeof.FLOAT), a.withByteOffset(1L*aoffset*Sizeof.FLOAT), 1L*len*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyDeviceToHost);
    cudaStreamSynchronize(Mat.SyncMethod);
    val err = cudaGetLastError;
    if (err != 0) throw new RuntimeException(msg +" error in memcpy "+ cudaGetErrorString(err));
  }
  
  def CPUtoGPUarraycopy(a:Array[Float], aoffset:Int,  b:Pointer, boffset:Int, len:Int, msg:String ) = {
    cudaMemcpy(b.withByteOffset(1L*boffset*Sizeof.FLOAT), Pointer.to(a).withByteOffset(1L*aoffset*Sizeof.FLOAT), 1L*len*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyHostToDevice);
    cudaStreamSynchronize(Mat.SyncMethod);
    val err = cudaGetLastError;
    if (err != 0) throw new RuntimeException(msg +" error in memcpy "+ cudaGetErrorString(err));
  } 
  
  def newOrCheckGMat(nr:Int, nc:Int, outmat:Mat):GMat = {
    if (outmat.asInstanceOf[AnyRef] == null || (outmat.nrows == 0 && outmat.ncols == 0)) {
      GMat(nr, nc)
    } else {
      outmat match {
        case omat:GMat => if (omat.nrows != nr || omat.ncols != nc) {
        omat.recycle(nr, nc, 0)
      } else {
      	omat
      }
      }
    }
  }  
  
  def newOrCheckGMat(dims:Array[Int], out:Mat):GMat = {
    if (out.asInstanceOf[AnyRef] != null && ND.compareDims(out.dims.data, dims)) {
      out.asInstanceOf[GMat]
    } else {
      GMat.make(dims)
    }
  }
  
  def newOrCheckGMat(dims:IMat, out:Mat):GMat = newOrCheckGMat(dims.data, out);
    
  def newOrCheckGMat(nr:Int, nc:Int, outmat:Mat, matGuid:Long, opHash:Int):GMat = {
    val m = if (outmat.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
      newOrCheckGMat(nr, nc, outmat)
    } else {
      val key = (matGuid, opHash.toLong, SciFunctions.getGPU)
      val res = Mat.cache3(key)
      if (res != null) {
        newOrCheckGMat(nr, nc, res)
      } else {
        val omat = newOrCheckGMat(nr, nc, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
    if (m.myGPU != SciFunctions.getGPU) {
    	throw new RuntimeException("newOrCheckGMat1 problem with mat %d" format m.GUID)
    }
    m
  }
  
  def newOrCheckGMat(dims:Array[Int], out:Mat, matGuid:Long, opHash:Int):GMat = {
    if (out.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
       newOrCheckGMat(dims, out)
    } else {
      val key = (matGuid, opHash.toLong, SciFunctions.getGPU)
      val res = Mat.cache3(key)
      if (res != null) {
        newOrCheckGMat(dims, res)
      } else {
        val omat = newOrCheckGMat(dims, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckGMat(dims:IMat, out:Mat, g1:Long, opHash:Int):GMat = newOrCheckGMat(dims.data, out, g1, opHash);
  
  def newOrCheckGMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int):GMat = {
    val m = if (outmat.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
      newOrCheckGMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, opHash.toLong, SciFunctions.getGPU)
      val res = Mat.cache4(key)
      if (res != null) {
      	newOrCheckGMat(nr, nc, res)
      } else {
        val omat = newOrCheckGMat(nr, nc, null)
        Mat.cache4put(key, omat)
        omat
      }
    }
    if (m.myGPU != SciFunctions.getGPU) {
    	throw new RuntimeException("newOrCheckGMat2 problem with mat %d" format m.GUID)
    }
    m
  }
  
  def newOrCheckGMat(dims:Array[Int], out:Mat, guid1:Long, guid2:Long, opHash:Int):GMat = {
    if (out.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
      newOrCheckGMat(dims, out)
    } else {
      val key = (guid1, guid2, opHash.toLong, SciFunctions.getGPU)
      val res = Mat.cache4(key)
      if (res != null) {
        newOrCheckGMat(dims, res)
      } else {
        val omat = newOrCheckGMat(dims, null)
        Mat.cache4put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckGMat(dims0:IMat, out:Mat, g1:Long, g2:Long, opHash:Int):GMat = {
    newOrCheckGMat(dims0.data, out, g1, g2, opHash);
  }
    
  def newOrCheckGMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):GMat = {
    val m = if (outmat.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
    	newOrCheckGMat(nr, nc, outmat)
    } else {
    	val key = (guid1, guid2, guid3, opHash.toLong, SciFunctions.getGPU)
        val res = Mat.cache5(key)
    	if (res != null) {
    		newOrCheckGMat(nr, nc, res)
    	} else {
    		val omat = newOrCheckGMat(nr, nc, null)
    		Mat.cache5put(key, omat)
    		omat
    	}
    }
    if (m.myGPU != SciFunctions.getGPU) {
    	throw new RuntimeException("newOrCheckGMat3 problem with mat %d" format m.GUID)
    }
    m
  }
  
  def newOrCheckGMat(dims:Array[Int], out:Mat, g1:Long, g2:Long, g3:Long, opHash:Int):GMat = {
    if (out.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
      newOrCheckGMat(dims, out)
    } else {
      val key = (g1, g2, g3, opHash.toLong, SciFunctions.getGPU)
      val res = Mat.cache5(key)
      if (res != null) {
        newOrCheckGMat(dims, res)
      } else {
        val omat = newOrCheckGMat(dims, null)
        Mat.cache5put(key, omat)
        omat
      }
    }
  }

  def newOrCheckGMat(dims:IMat, out:Mat, g1:Long, g2:Long, g3:Long, opHash:Int):GMat = newOrCheckGMat(dims.data, out, g1, g2, g3, opHash);

  def newOrCheckGMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, guid4:Long, opHash:Int):GMat = {
    val m = if (outmat.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
    	newOrCheckGMat(nr, nc, outmat)
    } else {
    	val key = (guid1, guid2, guid3, guid4, opHash.toLong, SciFunctions.getGPU)
        val res = Mat.cache6(key)
    	if (res != null) {
    		newOrCheckGMat(nr, nc, res)
    	} else {
    		val omat = newOrCheckGMat(nr, nc, null)
    		Mat.cache6put(key, omat)
    		omat
    	}
    }
    if (m.myGPU != SciFunctions.getGPU) {
    	throw new RuntimeException("newOrCheckGMat3 problem with mat %d" format m.GUID)
    }
    m
  }

  def newOrCheckGMat(dims:Array[Int], out:Mat, g1:Long, g2:Long, g3:Long, g4:Long, opHash:Int):GMat = {
    if (out.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
      newOrCheckGMat(dims, out)
    } else {
      val key = (g1, g2, g3, g4, opHash.toLong, SciFunctions.getGPU)
      val res = Mat.cache6(key)
      if (res != null) {
        newOrCheckGMat(dims, res)
      } else {
        val omat = newOrCheckGMat(dims, null)
        Mat.cache6put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckGMat(dims:IMat, out:Mat, g1:Long, g2:Long, g3:Long, g4:Long, opHash:Int):GMat = newOrCheckGMat(dims.data, out, g1, g2, g3, g4, opHash);

}







