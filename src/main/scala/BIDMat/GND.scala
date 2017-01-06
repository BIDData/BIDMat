//-*-coding:utf-8-*-
// N-dimensional arrays with floating point contents
package BIDMat
import MatFunctions._
import jcuda._
import jcuda.runtime._
import jcuda.runtime.JCuda._
import jcuda.runtime.cudaMemcpyKind._
import jcuda.jcublas._
import jcuda.jcublas.JCublas._
import edu.berkeley.bid.CUMAT
import GMat._
import java.util.Arrays
import java.util.concurrent.atomic._
import scala.concurrent.future
import scala.concurrent.ExecutionContext.Implicits.global
import scala.util.hashing.MurmurHash3


case class GND(dims0:Array[Int], val data:Pointer) extends ND(dims0) { 
  
	override def mytype = "GND"
	  
	override def dv:Double =
    if (nrows > 1 || ncols > 1) {
      throw new RuntimeException("Matrix should be 1x1 to extract value")
    } else {
      toFND(null).data(0)
    }
  
  override def fv:Float =
    if (nrows > 1 || ncols > 1) {
      throw new RuntimeException("Matrix should be 1x1 to extract value")
    } else {
      toFND(null).data(0)
    }
  
  override def contents() = {
    val out = new GMat(length, 1, data, length);
    out.setGUID(MurmurHash3.mix(MurmurHash3.mix(length, 1), (GUID*7897889).toInt));
    out
  }

  def applyf(indx:Int):Float = { 
    if (indx >= 0 && indx < length) { 
    	val tmp = new Array[Float](1);
    	cudaMemcpy(Pointer.to(tmp), data.withByteOffset(1L*indx*Sizeof.FLOAT), Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyDeviceToHost);
    	var err = cudaGetLastError;
    	if (err != 0) throw new RuntimeException("GND apply() error " + cudaGetErrorString(err));
    	tmp(0);
    } else { 
      throw new RuntimeException("GND index out of range")
    }
  }
  
    
  def apply(i:Int):Float = applyf(i);

  def apply(i1:Int, i2:Int):Float = apply(Array(i1, i2))
  def apply(i1:Int, i2:Int, i3:Int):Float = apply(Array(i1, i2, i3))
  def apply(i1:Int, i2:Int, i3:Int, i4:Int):Float = apply(Array(i1, i2, i3, i4))
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int):Float = apply(Array(i1, i2, i3, i4, i5))
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int):Float = apply(Array(i1, i2, i3, i4, i5, i6))
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int):Float = apply(Array(i1, i2, i3, i4, i5, i6, i7))
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, i8:Int):Float = apply(Array(i1, i2, i3, i4, i5, i6, i7, i8))
  
  def apply(inds:List[Int]):Float = apply(inds.toArray)
  
  def apply(inds:Array[Int]):Float = {
  	val indx = ND.linearize(inds, dims);
  	val tmp = new Array[Float](1);
    GPUtoCPUarraycopy(data, indx, tmp, 0, 1, "GND apply");
  	tmp(0);
  }
  
  def apply(inds:GIMat):GMat = {
  	inds match {
  	case aa:MatrixWildcard => {
  	  val out = GMat.newOrCheckGMat(length, 1, null, GUID, inds.GUID, "apply(?)".##);
  	  GPUtoGPUarraycopy(data, 0, out.data, 0, length, "GND apply(?)");
  		out
  	}
  	case _ => {
      val out = GMat.newOrCheckGMat(inds.nrows, inds.ncols, null, GUID, inds.GUID, "applyI".##);
      val err = CUMAT.copyFromInds(data, out.data, inds.data, inds.llength);
      if (err != 0) {
        throw new RuntimeException("GND apply error " + cudaGetErrorString(err))
      }
  		out;
  	}
  	}
  }
  
  def apply(inds0:List[GIMat]):GND = apply(inds0.toArray);
  
  def safePointer(ind:GIMat):Pointer = {
    if (ind.asInstanceOf[AnyRef] == null) {
      GMat.nullPointer;
    } else {
      ind.data;
    }
  }
  
  def apply(inds:Array[GIMat]):GND = {  
    val newdims = new Array[Int](_dims.length)
    val newinds = new Array[GIMat](_dims.length)
    var j = 0
    for (i <- 0 until _dims.length) {
      inds(i) match {
        case aa:MatrixWildcard => {
          newdims(i) = _dims(i); 
        }
        case _ => {
          newdims(i) = inds(i).length;
          newinds(i) = inds(i)
        }
      }
    }
    val out = GND.newOrCheckGND(newdims, null, GUID, ND.hashGUIDs(inds), "apply".##);
    inds.length match {
      case 2 => {
    	  val omat = out.toGMatView(newdims(0), newdims(1));
    	  toGMatView(dims(0), dims(1)).applyx(inds(0), inds(1), omat);
      }
      case 3 => {
        val err = CUMAT.copyFromInds3D(data, dims(0), dims(1), out.data, newdims(0), newdims(1), 
            safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1), safePointer(newinds(2)), newdims(2));
        if (err != 0) throw new RuntimeException("GND apply(I, J, K) error" + cudaGetErrorString(err));
      }
      case 4 => {
      	val err = CUMAT.copyFromInds4D(data, dims(0), dims(1), dims(2), out.data, newdims(0), newdims(1), newdims(2),
      			safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1), safePointer(newinds(2)), newdims(2), safePointer(newinds(3)), newdims(3));
      	if (err != 0) throw new RuntimeException("GND apply(I, J, K, L) error" + cudaGetErrorString(err));
      }   
      case _ => throw new RuntimeException("GND slice access with more than 4 indices not supported");
    }
    out;
  }
  
  def apply(i1:GIMat, i2:GIMat):GND = apply(Array(i1, i2));
  def apply(i1:GIMat, i2:GIMat, i3:GIMat):GND = apply(Array(i1, i2, i3));
  def apply(i1:GIMat, i2:GIMat, i3:GIMat, i4:GIMat):GND = apply(Array(i1, i2, i3, i4));
  
  def apply (a1:IMat, a2:GIMat):GND = apply(GIMat(a1), a2);
  def apply(a1:GIMat, a2:IMat):GND = apply(a1, GIMat(a2));
  def apply(a1:IMat, a2:IMat):GND = apply(GIMat(a1), GIMat(a2));
  
  def apply(a1:IMat, a2:GIMat, a3:GIMat):GND = apply(GIMat(a1), a2, a3);
  def apply(a1:GIMat, a2:IMat, a3:GIMat):GND = apply(a1, GIMat(a2), a3);
  def apply(a1:IMat, a2:IMat, a3:GIMat):GND = apply(GIMat(a1), GIMat(a2), a3);
  def apply(a1:GIMat, a2:GIMat, a3:IMat):GND = apply(a1, a2, GIMat(a3));
  def apply(a1:IMat, a2:GIMat, a3:IMat):GND = apply(GIMat(a1), a2, GIMat(a3));
  def apply(a1:GIMat, a2:IMat, a3:IMat):GND = apply(a1, GIMat(a2), GIMat(a3));
  def apply(a1:IMat, a2:IMat, a3:IMat):GND = apply(GIMat(a1), GIMat(a2), GIMat(a3));

  def apply(a1:IMat, a2:GIMat, a3:GIMat, a4:GIMat):GND = apply(GIMat(a1), a2, a3, a4);
  def apply(a1:GIMat, a2:IMat, a3:GIMat, a4:GIMat):GND = apply(a1, GIMat(a2), a3, a4);
  def apply(a1:IMat, a2:IMat, a3:GIMat, a4:GIMat):GND = apply(GIMat(a1), GIMat(a2), a3, a4);
  def apply(a1:GIMat, a2:GIMat, a3:IMat, a4:GIMat):GND = apply(a1, a2, GIMat(a3), a4);
  def apply(a1:IMat, a2:GIMat, a3:IMat, a4:GIMat):GND = apply(GIMat(a1), a2, GIMat(a3), a4);
  def apply(a1:GIMat, a2:IMat, a3:IMat, a4:GIMat):GND = apply(a1, GIMat(a2), GIMat(a3), a4);
  def apply(a1:IMat, a2:IMat, a3:IMat, a4:GIMat):GND = apply(GIMat(a1), GIMat(a2), GIMat(a3), a4);
  def apply(a1:GIMat, a2:GIMat, a3:GIMat, a4:IMat):GND = apply(a1, a2, a3, GIMat(a4));
  def apply(a1:IMat, a2:GIMat, a3:GIMat, a4:IMat):GND = apply(GIMat(a1), a2, a3, GIMat(a4));
  def apply(a1:GIMat, a2:IMat, a3:GIMat, a4:IMat):GND = apply(a1, GIMat(a2), a3, GIMat(a4));
  def apply(a1:IMat, a2:IMat, a3:GIMat, a4:IMat):GND = apply(GIMat(a1), GIMat(a2), a3, GIMat(a4));
  def apply(a1:GIMat, a2:GIMat, a3:IMat, a4:IMat):GND = apply(a1, a2, GIMat(a3), GIMat(a4));
  def apply(a1:IMat, a2:GIMat, a3:IMat, a4:IMat):GND = apply(GIMat(a1), a2, GIMat(a3), GIMat(a4));
  def apply(a1:GIMat, a2:IMat, a3:IMat, a4:IMat):GND = apply(a1, GIMat(a2), GIMat(a3), GIMat(a4));
  def apply(a1:IMat, a2:IMat, a3:IMat, a4:IMat):GND = apply(GIMat(a1), GIMat(a2), GIMat(a3), GIMat(a4));
  
  override def apply(inds:Mat):Mat = {
    inds match {
      case ii:GIMat => apply(ii);
      case ii:IMat => apply(GIMat(ii));
    }
  }
  
  override def apply(i1:Mat, i2:Mat):ND = {
    (i1, i2) match {
      case (a1:GIMat, a2:GIMat) => apply(a1, a2);
      case (a1:IMat, a2:GIMat) => apply(GIMat(a1), a2);
      case (a1:GIMat, a2:IMat) => apply(a1, GIMat(a2));
      case (a1:IMat, a2:IMat) => apply(GIMat(a1), GIMat(a2));
    }
  }
  
  override def apply(i1:Mat, i2:Mat, i3:Mat):ND = {
    (i1, i2, i3) match {
      case (a1:GIMat, a2:GIMat, a3:GIMat) => apply(a1, a2, a3);
      case (a1:IMat, a2:GIMat, a3:GIMat) => apply(GIMat(a1), a2, a3);
      case (a1:GIMat, a2:IMat, a3:GIMat) => apply(a1, GIMat(a2), a3);
      case (a1:IMat, a2:IMat, a3:GIMat) => apply(GIMat(a1), GIMat(a2), a3);
      case (a1:GIMat, a2:GIMat, a3:IMat) => apply(a1, a2, GIMat(a3));
      case (a1:IMat, a2:GIMat, a3:IMat) => apply(GIMat(a1), a2, GIMat(a3));
      case (a1:GIMat, a2:IMat, a3:IMat) => apply(a1, GIMat(a2), GIMat(a3));
      case (a1:IMat, a2:IMat, a3:IMat) => apply(GIMat(a1), GIMat(a2), GIMat(a3));
    }
  }
  
  override def apply(i1:Mat, i2:Mat, i3:Mat, i4:Mat):ND = {
    (i1, i2, i3, i4) match {
      case (a1:GIMat, a2:GIMat, a3:GIMat, a4:GIMat) => apply(a1, a2, a3, a4);
      case (a1:IMat, a2:GIMat, a3:GIMat, a4:GIMat) => apply(GIMat(a1), a2, a3, a4);
      case (a1:GIMat, a2:IMat, a3:GIMat, a4:GIMat) => apply(a1, GIMat(a2), a3, a4);
      case (a1:IMat, a2:IMat, a3:GIMat, a4:GIMat) => apply(GIMat(a1), GIMat(a2), a3, a4);
      case (a1:GIMat, a2:GIMat, a3:IMat, a4:GIMat) => apply(a1, a2, GIMat(a3), a4);
      case (a1:IMat, a2:GIMat, a3:IMat, a4:GIMat) => apply(GIMat(a1), a2, GIMat(a3), a4);
      case (a1:GIMat, a2:IMat, a3:IMat, a4:GIMat) => apply(a1, GIMat(a2), GIMat(a3), a4);
      case (a1:IMat, a2:IMat, a3:IMat, a4:GIMat) => apply(GIMat(a1), GIMat(a2), GIMat(a3), a4);
      case (a1:GIMat, a2:GIMat, a3:GIMat, a4:IMat) => apply(a1, a2, a3, GIMat(a4));
      case (a1:IMat, a2:GIMat, a3:GIMat, a4:IMat) => apply(GIMat(a1), a2, a3, GIMat(a4));
      case (a1:GIMat, a2:IMat, a3:GIMat, a4:IMat) => apply(a1, GIMat(a2), a3, GIMat(a4));
      case (a1:IMat, a2:IMat, a3:GIMat, a4:IMat) => apply(GIMat(a1), GIMat(a2), a3, GIMat(a4));
      case (a1:GIMat, a2:GIMat, a3:IMat, a4:IMat) => apply(a1, a2, GIMat(a3), GIMat(a4));
      case (a1:IMat, a2:GIMat, a3:IMat, a4:IMat) => apply(GIMat(a1), a2, GIMat(a3), GIMat(a4));
      case (a1:GIMat, a2:IMat, a3:IMat, a4:IMat) => apply(a1, GIMat(a2), GIMat(a3), GIMat(a4));
      case (a1:IMat, a2:IMat, a3:IMat, a4:IMat) => apply(GIMat(a1), GIMat(a2), GIMat(a3), GIMat(a4));
    }
  }

  
  def reshape(newdims:Int*):GND = reshape(newdims.toArray)
  
  def reshape(newdims:Array[Int]):GND = {
    if (newdims.reduce(_*_) == length) {
      val out = GND.newOrCheckGND(newdims, null, GUID, ND.hashInts(newdims), "reshape".##)
      GPUtoGPUarraycopy(data, 0, out.data, 0, length, "GND reshape");
      out
    } else {
      throw new RuntimeException("GND reshape total length doesnt match")
    }
  }
  
  def reshapeView(newdims:Int*):GND = reshapeView(newdims.toArray)
  
  def reshapeView(newdims:Array[Int]):GND = {
    if (newdims.reduce(_*_) == length) {
      val out = GND(newdims, data)
      out
    } else {
      throw new RuntimeException("GND reshapeView total length doesnt match")
    }
  }
  
  def toGMat(nr:Int, nc:Int):GMat = toGMat(nr, nc, null);
  
  def toGMat(nr:Int, nc:Int, omat:Mat):GMat = {
    if (nr*nc != length) {
      throw new RuntimeException("GND and output GMat dims dont match")
    } else {
    	val out = GMat.newOrCheckGMat(nr, nc, omat, GUID, nr, nc, "toGMat".##);
      GPUtoGPUarraycopy(data, 0, out.data, 0, length, "GND toGMat");
      out
    }
  }
  
   def toFND(omat:FND):FND = {
      val out = FND.newOrCheckFND(_dims, omat, GUID, "toFND".##);
      GPUtoCPUarraycopy(data, 0, out.data, 0, length, "GND toFND");
      out
   }

  def toGMatView(nr:Int, nc:Int):GMat = {
    if (nr*nc != length) {
      throw new RuntimeException("GND output GMat dims dont match")
    } else {
      val out = new GMat(nr, nc, data, nr * nc);
      out.setGUID(edu.berkeley.bid.MurmurHash3.MurmurHash3_x64_64(Array(GUID, nr), "toGMatView".##));
      out;
    }
  }

  def update(indx:Int, v:Float) = { 
    if (indx < 0 || indx >= length) {
    	throw new RuntimeException("GND index out of range")
    } else { 
      val tmp = Array[Float](v);
      CPUtoGPUarraycopy(tmp, 0, data, indx, 1, "GND update");
    }
    this
  }
  
  def updatex(inds:GIMat, vv:GMat):GND = {
  	inds match {
  	case aa:MatrixWildcard => {
  		if (vv.length == length) {
  			GPUtoGPUarraycopy(vv.data, 0, data, 0, length, "GND udpate");
  			this;
  	  } else {
  	    throw new RuntimeException("GND update(?) RHS dimension doesnt match")
  	  }
  	}
  	case _ => {
  		if (inds.length != vv.length) {
  			throw new RuntimeException("GND updatex error: I and v have unequal lengths " + inds.length + " and " + vv.length + ", respectively.")
  		}
  		val err = CUMAT.copyToInds(data, vv.data, inds.data, inds.llength);
  		if (err != 0) {
  			throw new RuntimeException("GND update error " + cudaGetErrorString(err))
  		}
  		this;
  	}
  	}
  }
  
  def update(inds:GIMat, vv:GMat):GND = updatex(inds:GIMat, vv:GMat):GND;
  
  def update(inds:List[Int], v:Float):GND = update(inds.toArray, v)
  
  def update(inds:Array[Int], v:Float):GND = {
  	val indx = ND.linearize(inds, dims); 
  	val tmp = Array[Float](v);
  	CPUtoGPUarraycopy(tmp, 0, data, indx, 1, "GND update");
  	this
  }
  
  def update(i1:Int, i2:Int, vv:Float):GND = update(Array(i1, i2), vv)
  def update(i1:Int, i2:Int, i3:Int, vv:Float):GND = update(Array(i1, i2, i3), vv)
  def update(i1:Int, i2:Int, i3:Int, i4:Int, vv:Float):GND = update(Array(i1, i2, i3, i4), vv)
  
  def update(inds:Array[GIMat], vv:GND):GND = {
    if (inds.length != _dims.length) {
      throw new RuntimeException("GND update wrong number of dims")
    }
    val newdims = new Array[Int](_dims.length)
    val newinds = new Array[GIMat](_dims.length)
    var j = 0
    for (i <- 0 until _dims.length) {
      inds(i) match {
        case aa:MatrixWildcard => {
          newdims(i) = _dims(i); 
        }
        case _ => {
          newdims(i) = inds(i).length;
          newinds(i) = inds(i)
        }
      }
    }
    ND.checkDims("GND update:", ND.trimDims(newdims), ND.trimDims(vv._dims));
    inds.length match {
    case 2 => {
    	val omat = vv.toGMatView(newdims(0), newdims(1));
    	toGMatView(dims(0), dims(1)).updatex(inds(0), inds(1), omat);
    }
    case 3 => {
    	val err = CUMAT.copyToInds3D(vv.data, vv.dims(0), vv.dims(1), data, dims(0), dims(1), 
    			safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1), safePointer(newinds(2)), newdims(2));
    	if (err != 0) throw new RuntimeException("GND apply(I, J, K) error" + cudaGetErrorString(err));
    }
    case 4 => {
    	val err = CUMAT.copyToInds4D(vv.data, vv.dims(0), vv.dims(1), vv.dims(2), data, dims(0), dims(1), dims(2),
    			safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1), safePointer(newinds(2)), newdims(2), safePointer(newinds(3)), newdims(3));
    	if (err != 0) throw new RuntimeException("GND apply(I, J, K, L) error" + cudaGetErrorString(err));
    }
    case _ => throw new RuntimeException("GND slice access with more than 4 indices not supported");
    }
    this
  }
  
  def update(i1:GIMat, i2:GIMat, vv:GND):GND = update(Array(i1, i2), vv);
  def update(i1:GIMat, i2:GIMat, i3:GIMat, vv:GND):GND = update(Array(i1, i2, i3), vv);
  def update(i1:GIMat, i2:GIMat, i3:GIMat, i4:GIMat, vv:GND):GND = update(Array(i1, i2, i3, i4), vv);
  
  def update(i1:GIMat, i2:GIMat, vv:ND):GND = update(Array(i1, i2), vv.asInstanceOf[GND]);
  def update(i1:GIMat, i2:GIMat, i3:GIMat, vv:ND):GND = update(Array(i1, i2, i3), vv.asInstanceOf[GND]);
  def update(i1:GIMat, i2:GIMat, i3:GIMat, i4:GIMat, vv:ND):GND = update(Array(i1, i2, i3, i4), vv.asInstanceOf[GND]);
  
  def update(a1:IMat, a2:GIMat, uu:ND):GND = update(GIMat(a1), a2, uu);
  def update(a1:GIMat, a2:IMat, uu:ND):GND = update(a1, GIMat(a2), uu);
  def update(a1:IMat, a2:IMat, uu:ND):GND = update(GIMat(a1), GIMat(a2), uu);

  def update(a1:IMat, a2:GIMat, a3:GIMat, uu:ND):GND = update(GIMat(a1), a2, a3, uu);
  def update(a1:GIMat, a2:IMat, a3:GIMat, uu:ND):GND = update(a1, GIMat(a2), a3, uu);
  def update(a1:IMat, a2:IMat, a3:GIMat, uu:ND):GND = update(GIMat(a1), GIMat(a2), a3, uu);
  def update(a1:GIMat, a2:GIMat, a3:IMat, uu:ND):GND = update(a1, a2, GIMat(a3), uu);
  def update(a1:IMat, a2:GIMat, a3:IMat, uu:ND):GND = update(GIMat(a1), a2, GIMat(a3), uu);
  def update(a1:GIMat, a2:IMat, a3:IMat, uu:ND):GND = update(a1, GIMat(a2), GIMat(a3), uu);
  def update(a1:IMat, a2:IMat, a3:IMat, uu:ND):GND = update(GIMat(a1), GIMat(a2), GIMat(a3), uu);

  def update(a1:IMat, a2:GIMat, a3:GIMat, a4:GIMat, uu:ND):GND = update(GIMat(a1), a2, a3, a4, uu);
  def update(a1:GIMat, a2:IMat, a3:GIMat, a4:GIMat, uu:ND):GND = update(a1, GIMat(a2), a3, a4, uu);
  def update(a1:IMat, a2:IMat, a3:GIMat, a4:GIMat, uu:ND):GND = update(GIMat(a1), GIMat(a2), a3, a4, uu);
  def update(a1:GIMat, a2:GIMat, a3:IMat, a4:GIMat, uu:ND):GND = update(a1, a2, GIMat(a3), a4, uu);
  def update(a1:IMat, a2:GIMat, a3:IMat, a4:GIMat, uu:ND):GND = update(GIMat(a1), a2, GIMat(a3), a4, uu);
  def update(a1:GIMat, a2:IMat, a3:IMat, a4:GIMat, uu:ND):GND = update(a1, GIMat(a2), GIMat(a3), a4, uu);
  def update(a1:IMat, a2:IMat, a3:IMat, a4:GIMat, uu:ND):GND = update(GIMat(a1), GIMat(a2), GIMat(a3), a4, uu);
  def update(a1:GIMat, a2:GIMat, a3:GIMat, a4:IMat, uu:ND):GND = update(a1, a2, a3, GIMat(a4), uu);
  def update(a1:IMat, a2:GIMat, a3:GIMat, a4:IMat, uu:ND):GND = update(GIMat(a1), a2, a3, GIMat(a4), uu);
  def update(a1:GIMat, a2:IMat, a3:GIMat, a4:IMat, uu:ND):GND = update(a1, GIMat(a2), a3, GIMat(a4), uu);
  def update(a1:IMat, a2:IMat, a3:GIMat, a4:IMat, uu:ND):GND = update(GIMat(a1), GIMat(a2), a3, GIMat(a4), uu);
  def update(a1:GIMat, a2:GIMat, a3:IMat, a4:IMat, uu:ND):GND = update(a1, a2, GIMat(a3), GIMat(a4), uu);
  def update(a1:IMat, a2:GIMat, a3:IMat, a4:IMat, uu:ND):GND = update(GIMat(a1), a2, GIMat(a3), GIMat(a4), uu);
  def update(a1:GIMat, a2:IMat, a3:IMat, a4:IMat, uu:ND):GND = update(a1, GIMat(a2), GIMat(a3), GIMat(a4), uu);
  def update(a1:IMat, a2:IMat, a3:IMat, a4:IMat, uu:ND):GND = update(GIMat(a1), GIMat(a2), GIMat(a3), GIMat(a4), uu);
  
  def update(i1:Mat, vv:Mat):GND = {
    (i1, vv) match {
      case (a1:GIMat, uu:GMat) => update(a1, uu);
      case (a1:GIMat, uu:FMat) => update(a1, GMat(uu));
      case (a1:IMat, uu:FMat) => update(GIMat(a1), GMat(uu));
      case (a1:IMat, uu:GMat) => update(GIMat(a1), uu);
    }
  }
  
  def update(i1:Mat, i2:Mat, vv:ND):GND = {
    (i1, i2, vv) match {
      case (a1:GIMat, a2:GIMat, uu:ND) => update(a1, a2, uu);
      case (a1:IMat, a2:GIMat, uu:ND) => update(GIMat(a1), a2, uu);
      case (a1:GIMat, a2:IMat, uu:ND) => update(a1, GIMat(a2), uu);
      case (a1:IMat, a2:IMat, uu:ND) => update(GIMat(a1), GIMat(a2), uu);
    }
  }
    
  def update(i1:Mat, i2:Mat, i3:Mat, vv:ND):GND = {
    (i1, i2, i3, vv) match {
      case (a1:GIMat, a2:GIMat, a3:GIMat, uu:ND) => update(a1, a2, a3, uu);
      case (a1:IMat, a2:GIMat, a3:GIMat, uu:ND) => update(GIMat(a1), a2, a3, uu);
      case (a1:GIMat, a2:IMat, a3:GIMat, uu:ND) => update(a1, GIMat(a2), a3, uu);
      case (a1:IMat, a2:IMat, a3:GIMat, uu:ND) => update(GIMat(a1), GIMat(a2), a3, uu);
      case (a1:GIMat, a2:GIMat, a3:IMat, uu:ND) => update(a1, a2, GIMat(a3), uu);
      case (a1:IMat, a2:GIMat, a3:IMat, uu:ND) => update(GIMat(a1), a2, GIMat(a3), uu);
      case (a1:GIMat, a2:IMat, a3:IMat, uu:ND) => update(a1, GIMat(a2), GIMat(a3), uu);
      case (a1:IMat, a2:IMat, a3:IMat, uu:ND) => update(GIMat(a1), GIMat(a2), GIMat(a3), uu);
    }
  }
  
  def update(i1:Mat, i2:Mat, i3:Mat, i4:Mat, vv:ND):GND = {
    (i1, i2, i3, i4, vv) match {
      case (a1:GIMat, a2:GIMat, a3:GIMat, a4:GIMat, uu:ND) => update(a1, a2, a3, a4, uu);
      case (a1:IMat, a2:GIMat, a3:GIMat, a4:GIMat, uu:ND) => update(GIMat(a1), a2, a3, a4, uu);
      case (a1:GIMat, a2:IMat, a3:GIMat, a4:GIMat, uu:ND) => update(a1, GIMat(a2), a3, a4, uu);
      case (a1:IMat, a2:IMat, a3:GIMat, a4:GIMat, uu:ND) => update(GIMat(a1), GIMat(a2), a3, a4, uu);
      case (a1:GIMat, a2:GIMat, a3:IMat, a4:GIMat, uu:ND) => update(a1, a2, GIMat(a3), a4, uu);
      case (a1:IMat, a2:GIMat, a3:IMat, a4:GIMat, uu:ND) => update(GIMat(a1), a2, GIMat(a3), a4, uu);
      case (a1:GIMat, a2:IMat, a3:IMat, a4:GIMat, uu:ND) => update(a1, GIMat(a2), GIMat(a3), a4, uu);
      case (a1:IMat, a2:IMat, a3:IMat, a4:GIMat, uu:ND) => update(GIMat(a1), GIMat(a2), GIMat(a3), a4, uu);
      case (a1:GIMat, a2:GIMat, a3:GIMat, a4:IMat, uu:ND) => update(a1, a2, a3, GIMat(a4), uu);
      case (a1:IMat, a2:GIMat, a3:GIMat, a4:IMat, uu:ND) => update(GIMat(a1), a2, a3, GIMat(a4), uu);
      case (a1:GIMat, a2:IMat, a3:GIMat, a4:IMat, uu:ND) => update(a1, GIMat(a2), a3, GIMat(a4), uu);
      case (a1:IMat, a2:IMat, a3:GIMat, a4:IMat, uu:ND) => update(GIMat(a1), GIMat(a2), a3, GIMat(a4), uu);
      case (a1:GIMat, a2:GIMat, a3:IMat, a4:IMat, uu:ND) => update(a1, a2, GIMat(a3), GIMat(a4), uu);
      case (a1:IMat, a2:GIMat, a3:IMat, a4:IMat, uu:ND) => update(GIMat(a1), a2, GIMat(a3), GIMat(a4), uu);
      case (a1:GIMat, a2:IMat, a3:IMat, a4:IMat, uu:ND) => update(a1, GIMat(a2), GIMat(a3), GIMat(a4), uu);
      case (a1:IMat, a2:IMat, a3:IMat, a4:IMat, uu:ND) => update(GIMat(a1), GIMat(a2), GIMat(a3), GIMat(a4), uu);
    }
  }
  
  def update(inds:Array[GIMat], v:Float):GND = {
    val newdims = new Array[Int](dims.length)
    for (i <- 0 until dims.length) {
      newdims(i) = inds(i) match {case aa:MatrixWildcard => _dims(i); case _ => inds(i).length}
    }
//    updateHelper(inds, v, 0, inds.length-1)
    this
  }
  
  def update(i1:GIMat, i2:GIMat, vv:Float):GND = update(Array(i1, i2), vv)
  def update(i1:GIMat, i2:GIMat, i3:GIMat, vv:Float):GND = update(Array(i1, i2, i3), vv)
  def update(i1:GIMat, i2:GIMat, i3:GIMat, i4:GIMat, vv:Float):GND = update(Array(i1, i2, i3, i4), vv)
  
  def update(i1:Mat, i2:Mat, vv:Float):ND = {
    (i1, i2) match {
      case (a1:GIMat, a2:GIMat) => update(a1, a2, vv);
    }
  }
  
  def update(i1:Mat, vv:Float):ND = {
    (i1) match {
      case (a1:GIMat) => update(a1, vv);
    }
  }

  def update(i1:Mat, i2:Mat, i3:Mat, vv:Float):ND = {
    (i1, i2, i3) match {
      case (a1:GIMat, a2:GIMat, a3:GIMat) => update(a1, a2, a3, vv);
    }
  }
    
  def update(i1:Mat, i2:Mat, i3:Mat, i4:Mat, vv:Float):ND = {
    (i1, i2, i3, i4) match {
      case (a1:GIMat, a2:GIMat, a3:GIMat, a4:GIMat) => update(a1, a2, a3, a4, vv);
    }
  }
  
  val asMat:GMat = {
    val out = new GMat(nrows, ncols, data, nrows * ncols);
    out.setGUID(edu.berkeley.bid.MurmurHash3.MurmurHash3_x64_64(Array(GUID), 0x45239234));
    out
  }
    
 override def colslice(a:Int, b:Int, out:ND):GND = {
    val dims0 = dims;
    dims0(dims.length - 1) = b - a;
    val outx = GND.newOrCheckGND(dims0, out, GUID, a, b, "colslice".##);
    val omat = asMat.colslice(a, b, outx.asMat);
    outx;
  }
  
  override def colslice(a:Int, b:Int):GND = colslice(a, b, null)
  
  override def colslice(a:Int, b:Int, out:ND, c:Int):GND = {
    asMat.colslice(a, b, out.asMat, c);
    out.asInstanceOf[GND];
  }
  
  def copyTo(a:GND):GND = {
    GMat.GPUtoGPUarraycopy(data, 0, a.data, 0, length, "GND copyTo");
    a
  }
  
  def copyTo(a:FND):FND = {
    GMat.GPUtoCPUarraycopy(data, 0, a.data, 0, length, "GND copyTo");
    a
  }
  
  def copyTo(a:ND):ND = {
    a match {
      case aa:FND => copyTo(aa);
      case aa:GND => copyTo(aa);
    }
  }
  
  def zeros(nr:Int, nc:Int):GMat = GMat.zeros(nr, nc)
   
  def zeros(dims0:IMat):GND = GND.zeros(dims0);
 
  def zeros:GND = GND.zeros(dims);
  
  override def ones(dims0:IMat):GND = GND.ones(dims0);

  def clear:GND = {
	  cudaMemset(data, 0, Sizeof.FLOAT*length);
    cudaDeviceSynchronize;
    val err = cudaGetLastError
    if (err != 0) throw new RuntimeException("GND clear error");
    this
  }
  
  def set(v:Float):GND = {
    asMat.set(v);
    this
  }
  
  def transpose(dims:Array[Int]):GND = transpose(irow(dims))

  def transpose(perm:IMat):GND = { 
    val nd = _dims.length
    if (perm.length != nd) { 
      throw new RuntimeException("GND transpose bad permutation ")
    }
    val xdims = irow(_dims);
    val iperm = invperm(perm);
    val pdims = xdims(perm).data;
    var out = GND.newOrCheckGND(pdims, null, GUID, ND.hashInts(pdims), "transpose".##);
    var out2 = GND.newOrCheckGND(pdims, null, GUID, ND.hashInts(pdims), "transpose1".##);
    cudaMemcpy(out.data, data, 4L*length, cudaMemcpyDeviceToDevice);
    for (i <- (nd - 1) until 0 by -1) { 
      if (iperm(i) != i) { 
        val (d1, d2, d3) = ND.getDims(i, iperm, xdims)
        if (d1 > 1 && d2 > 1) { 
 //         println("spermute %d %d %d" format (d1,d2,d3))
          CUMAT.spermute(d1, d2, d3, out.data, out2.data)
          val tmp = out2
          out2 = out
          out = tmp
        }
        ND.rotate(i, iperm, xdims)
      } 
    }
    out
  }
  
  def transpose(i1:Int, i2:Int):GND = transpose(Array(i1, i2))
  def transpose(i1:Int, i2:Int, i3:Int):GND = transpose(Array(i1, i2, i3))
  def transpose(i1:Int, i2:Int, i3:Int, i4:Int):GND = transpose(Array(i1, i2, i3, i4))
  def transpose(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int):GND = transpose(Array(i1, i2, i3, i4, i5))
  def transpose(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int):GND = transpose(Array(i1, i2, i3, i4, i5, i6))
  def transpose(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int):GND = transpose(Array(i1, i2, i3, i4, i5, i6, i7))
  def transpose(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, i8:Int):GND = transpose(Array(i1, i2, i3, i4, i5, i6, i7, i8))
  
  
  
  
  def * (b : GND):GND = {
	  val dims0 = dims(0->(dims.length-1));
	  val dims1 = b.dims(1->b.dims.length);
	  val x = toGMatView(SciFunctions.prod(dims0).v, dims(dims.length-1)); 
	  val y = b.toGMatView(b.dims(0), SciFunctions.prod(dims1).v);
	  val zz = GND.newOrCheckGND((dims0 \ dims1).data, null, GUID, b.GUID, "*".##);
	  val z = zz.toGMatView(x.nrows, y.ncols);
	  z ~ x * y;      
	  zz
  }
  
  def unary_-():GND = {
    val zz = GND.newOrCheckGND(dims.data, null, GUID, "-".##);
    val a = toGMatView(length, 1);
    val b = zz.toGMatView(length, 1);
    val minusOne = GMat.newOrCheckGMat(1,1,null,-1,"minusOne".##);
    minusOne.set(-1f);
    a.gOp(minusOne, null, GMat.BinOp.op_mul);
    zz;
  }
   
  def + (mat:GND):GND = {val (a, b, c, d) = GND.asGMats(this, mat, null, "+"); c ~ a + b; d}
  def - (mat:GND):GND = {val (a, b, c, d) = GND.asGMats(this, mat, null, "-"); c ~ a - b; d}
  def *@ (mat:GND):GND = {val (a, b, c, d) = GND.asGMats(this, mat, null, "*@"); c ~ a *@ b; d}
  def / (mat:GND):GND = {val (a, b, c, d) = GND.asGMats(this, mat, null, "/"); c ~ a / b; d}
  def ^ (mat:GND):GND = {val (a, b, c, d) = GND.asGMats(this, mat, null, "^"); c ~ a ^ b; d}
  
  def > (mat:GND):GND = {val (a, b, c, d) = GND.asGMats(this, mat, null, ">"); c ~ a > b; d}
  def < (mat:GND):GND = {val (a, b, c, d) = GND.asGMats(this, mat, null, "<"); c ~ a < b; d}
  def >= (mat:GND):GND = {val (a, b, c, d) = GND.asGMats(this, mat, null, ">="); c ~ a >= b; d}
  def <= (mat:GND):GND = {val (a, b, c, d) = GND.asGMats(this, mat, null, "<="); c ~ a <= b; d}
  def != (mat:GND):GND = {val (a, b, c, d) = GND.asGMats(this, mat, null, "!="); c ~ a != b; d}
  def == (mat:GND):GND = {val (a, b, c, d) = GND.asGMats(this, mat, null, "=="); c ~ a == b; d}
  def === (mat:GND):GND = {val (a, b, c, d) = GND.asGMats(this, mat, null, "==="); c ~ a === b; d}
     
  def max (mat:GND):GND = {val (a, b, c, d) = GND.asGMats(this, mat, null, "max"); SciFunctions.max(a, b, c); d}
  def max (mat:GND, omat:ND):GND = {val (a, b, c, d) = GND.asGMats(this, mat, omat, "max"); SciFunctions.max(a, b, c); d}
  def min (mat:GND):GND = {val (a, b, c, d) = GND.asGMats(this, mat, null, "max"); SciFunctions.min(a, b, c); d}
  def min (mat:GND, omat:ND):GND = {val (a, b, c, d) = GND.asGMats(this, mat, omat, "max"); SciFunctions.min(a, b, c); d}
  
  
  override def + (b:Float):GND = {val (a, c, d) = GND.asGMats(this, null, "+"); c ~ a + b; d}
  override def - (b:Float):GND = {val (a, c, d) = GND.asGMats(this, null, "-"); c ~ a - b; d}
  override def *@ (b:Float):GND = {val (a, c, d) = GND.asGMats(this, null, "*@"); c ~ a *@ b; d}
  override def ∘ (b:Float):GND = {val (a, c, d) = GND.asGMats(this, null, "∘"); c ~ a *@ b; d}
  override def * (b:Float):GND = {val (a, c, d) = GND.asGMats(this, null, "*"); c ~ a *@ b; d}
  override def / (b:Float):GND = {val (a, c, d) = GND.asGMats(this, null, "/"); c ~ a / b; d}
  override def ^ (b:Float):GND = {val (a, c, d) = GND.asGMats(this, null, "^"); c ~ a ^ b; d}
  
  override def > (b:Float):GND = {val (a, c, d) = GND.asGMats(this, null, ">"); c ~ a > b; d}
  override def < (b:Float):GND = {val (a, c, d) = GND.asGMats(this, null, "<"); c ~ a < b; d}
  override def >= (b:Float):GND = {val (a, c, d) = GND.asGMats(this, null, ">="); c ~ a >= b; d}
  override def <= (b:Float):GND = {val (a, c, d) = GND.asGMats(this, null, "<="); c ~ a <= b; d}
  override def != (b:Float):GND = {val (a, c, d) = GND.asGMats(this, null, "!="); c ~ a != b; d}
  override def == (b:Float):GND = {val (a, c, d) = GND.asGMats(this, null, "=="); c ~ a == b; d}
  override def === (b:Float):GND = {val (a, c, d) = GND.asGMats(this, null, "==="); c ~ a === b; d}  
    
  override def max(b:Float):GND = {val (a, c, d) = GND.asGMats(this, null, "max"); SciFunctions.max(a, b, c); d}
  override def min(b:Float):GND = {val (a, c, d) = GND.asGMats(this, null, "min"); SciFunctions.min(a, b, c); d}
  def max(b:Float, omat:ND):GND = {val (a, c, d) = GND.asGMats(this, omat, "max"); SciFunctions.max(a, b, c); d}
  def min(b:Float, omat:ND):GND = {val (a, c, d) = GND.asGMats(this, omat, "min"); SciFunctions.min(a, b, c); d}
  
  
  
  override def + (b:Double):GND = {val (a, c, d) = GND.asGMats(this, null, "+"); c ~ a + b; d}
  override def - (b:Double):GND = {val (a, c, d) = GND.asGMats(this, null, "-"); c ~ a - b; d}
  override def *@ (b:Double):GND = {val (a, c, d) = GND.asGMats(this, null, "*@"); c ~ a *@ b; d}
  override def ∘ (b:Double):GND = {val (a, c, d) = GND.asGMats(this, null, "∘"); c ~ a *@ b; d}
  override def / (b:Double):GND = {val (a, c, d) = GND.asGMats(this, null, "/"); c ~ a / b; d}
  override def ^ (b:Double):GND = {val (a, c, d) = GND.asGMats(this, null, "^"); c ~ a ^ b; d}
  override def * (b:Double):GND = {val (a, c, d) = GND.asGMats(this, null, "*"); c ~ a *@ b; d}
  
  override def > (b:Double):GND = {val (a, c, d) = GND.asGMats(this, null, ">"); c ~ a > b; d}
  override def < (b:Double):GND = {val (a, c, d) = GND.asGMats(this, null, "<"); c ~ a < b; d}
  override def >= (b:Double):GND = {val (a, c, d) = GND.asGMats(this, null, ">="); c ~ a >= b; d}
  override def <= (b:Double):GND = {val (a, c, d) = GND.asGMats(this, null, "<="); c ~ a <= b; d}
  override def != (b:Double):GND = {val (a, c, d) = GND.asGMats(this, null, "!="); c ~ a != b; d}
  override def == (b:Double):GND = {val (a, c, d) = GND.asGMats(this, null, "=="); c ~ a == b; d}
  override def === (b:Double):GND = {val (a, c, d) = GND.asGMats(this, null, "==="); c ~ a === b; d}
  
  
  
  override def + (b:Int):GND = {val (a, c, d) = GND.asGMats(this, null, "+"); c ~ a + b; d}
  override def - (b:Int):GND = {val (a, c, d) = GND.asGMats(this, null, "-"); c ~ a - b; d}
  override def *@ (b:Int):GND = {val (a, c, d) = GND.asGMats(this, null, "*@"); c ~ a *@ b; d}
  override def ∘ (b:Int):GND = {val (a, c, d) = GND.asGMats(this, null, "∘"); c ~ a *@ b; d}
  override def / (b:Int):GND = {val (a, c, d) = GND.asGMats(this, null, "/"); c ~ a / b; d}
  override def ^ (b:Int):GND = {val (a, c, d) = GND.asGMats(this, null, "^"); c ~ a ^ b; d}
  override def * (b:Int):GND = {val (a, c, d) = GND.asGMats(this, null, "*"); c ~ a *@ b; d}
  
  override def > (b:Int):GND = {val (a, c, d) = GND.asGMats(this, null, ">"); c ~ a > b; d}
  override def < (b:Int):GND = {val (a, c, d) = GND.asGMats(this, null, "<"); c ~ a < b; d}
  override def >= (b:Int):GND = {val (a, c, d) = GND.asGMats(this, null, ">="); c ~ a >= b; d}
  override def <= (b:Int):GND = {val (a, c, d) = GND.asGMats(this, null, "<="); c ~ a <= b; d}
  override def != (b:Int):GND = {val (a, c, d) = GND.asGMats(this, null, "!="); c ~ a != b; d}
  override def == (b:Int):GND = {val (a, c, d) = GND.asGMats(this, null, "=="); c ~ a == b; d}
  override def === (b:Int):GND = {val (a, c, d) = GND.asGMats(this, null, "==="); c ~ a === b; d}
  
  
  override def + (b:Long):GND = {val (a, c, d) = GND.asGMats(this, null, "+"); c ~ a + b; d}
  override def - (b:Long):GND = {val (a, c, d) = GND.asGMats(this, null, "-"); c ~ a - b; d}
  override def *@ (b:Long):GND = {val (a, c, d) = GND.asGMats(this, null, "*@"); c ~ a *@ b; d}
  override def ∘ (b:Long):GND = {val (a, c, d) = GND.asGMats(this, null, "∘"); c ~ a *@ b; d}
  override def / (b:Long):GND = {val (a, c, d) = GND.asGMats(this, null, "/"); c ~ a / b; d}
  override def ^ (b:Long):GND = {val (a, c, d) = GND.asGMats(this, null, "^"); c ~ a ^ b; d}
  override def * (b:Long):GND = {val (a, c, d) = GND.asGMats(this, null, "*"); c ~ a *@ b; d}
  
  override def > (b:Long):GND = {val (a, c, d) = GND.asGMats(this, null, ">"); c ~ a > b; d}
  override def < (b:Long):GND = {val (a, c, d) = GND.asGMats(this, null, "<"); c ~ a < b; d}
  override def >= (b:Long):GND = {val (a, c, d) = GND.asGMats(this, null, ">="); c ~ a >= b; d}
  override def <= (b:Long):GND = {val (a, c, d) = GND.asGMats(this, null, "<="); c ~ a <= b; d}
  override def != (b:Long):GND = {val (a, c, d) = GND.asGMats(this, null, "!="); c ~ a != b; d}
  override def == (b:Long):GND = {val (a, c, d) = GND.asGMats(this, null, "=="); c ~ a == b; d}
  override def === (b:Long):GND = {val (a, c, d) = GND.asGMats(this, null, "==="); c ~ a === b; d}
  
  
  override def + (b : ND):ND = this + b.asInstanceOf[GND];
  override def - (b : ND):ND = this - b.asInstanceOf[GND]; 
  override def * (b : ND):ND = this * b.asInstanceOf[GND];
  override def *@ (b : ND):ND = this *@ b.asInstanceOf[GND];
  override def ∘  (b : ND):ND = this *@ b.asInstanceOf[GND];
  override def /  (b : ND):ND = this / b.asInstanceOf[GND];
  override def ^ (b : ND):ND = this ^ b.asInstanceOf[GND];
  
  override def > (b : ND):ND = this > b.asInstanceOf[GND];
  override def < (b : ND):ND = this < b.asInstanceOf[GND];
  override def >= (b : ND):ND = this >= b.asInstanceOf[GND];
  override def <= (b : ND):ND = this <= b.asInstanceOf[GND];
  override def == (b : ND):ND = this == b.asInstanceOf[GND];
  override def === (b : ND):ND = this === b.asInstanceOf[GND];
  override def != (b : ND):ND = this != b.asInstanceOf[GND];
  
  override def \ (b : ND):ND = this \ b.asInstanceOf[GND];
  override def on (b : ND):ND = this on b.asInstanceOf[GND];
  
  def reduce(inds:Array[Int], fctn:(GMat)=>GMat, opname:String):GND = {
    val alldims = izeros(_dims.length,1);
    val xinds = new IMat(inds.length, 1, inds);
    val xdims = new IMat(_dims.length, 1, _dims);
    alldims(xinds) = 1;
    if (SciFunctions.sum(alldims).v != inds.length) {
      throw new RuntimeException(opname+ " indices arent a legal subset of dims");
    }
    val restdims = find(alldims == 0);
    val tmp = transpose((xinds on restdims).data);
    val dimprodx = SciFunctions.prod(xdims(xinds)).v;
    val dimprody = SciFunctions.prod(xdims(restdims)).v;
    val tmpF = new GMat(dimprodx, dimprody, tmp.data, dimprodx * dimprody);
    val tmpSum:GMat = fctn(tmpF)
    val out1 = new GND((iones(inds.length,1) on xdims(restdims)).data, tmpSum.data)
    out1.transpose(invperm(xinds on restdims).data)
  }
  
  def sum(inds:Array[Int]):GND = reduce(inds, SciFunctions.sum, "sum")
  def prod(inds:Array[Int]):GND = reduce(inds, SciFunctions.prod, "prod")
  def mean(inds:Array[Int]):GND = reduce(inds, SciFunctions.mean, "mean")
  def variance(inds:Array[Int]):GND = reduce(inds, SciFunctions.variance, "variance")
  def maxi(inds:Array[Int]):GND = reduce(inds, SciFunctions.maxi, "maxi")
  def mini(inds:Array[Int]):GND = reduce(inds, SciFunctions.mini, "mini")
  
  override def sum(inds:Int*):GND = sum(inds.toArray)
  override def prod(inds:Int*):GND = prod(inds.toArray)
  override def mean(inds:Int*):GND = mean(inds.toArray)
  override def variance(inds:Int*):GND = variance(inds.toArray)
  override def maxi(inds:Int*):GND = maxi(inds.toArray)
  override def mini(inds:Int*):GND = mini(inds.toArray)  
    
  def ~ (b : GND):GNDPair = new GNDPair(this, b)
  def ~ (b : ND):GNDPair = new GNDPair(this, b.asInstanceOf[GND])

}

class GNDPair(val omat:ND, val amat:GND) extends NDPair {
  
	def * (b : GND):GND = {
			val dims0 = amat.dims(0->(amat.dims.length-1));
			val dims1 = b.dims(1->b.dims.length);
			val x = amat.toGMatView(SciFunctions.prod(dims0).v, amat.dims(amat.dims.length-1)); 
			val y = b.toGMatView(b.dims(0), SciFunctions.prod(dims1).v);
			val odims = dims0 \ dims1;
			ND.checkDims("GND *", odims, omat.dims);
			val zz = GND.newOrCheckGND(odims.data, omat, amat.GUID, b.GUID, "*".##);
			val z = zz.toGMatView(x.nrows, y.ncols);
			z ~ x * y;      
			zz
	}
	
	def *^ (b : GND):GND = {
     omat match {
       case gg:GFilter => {gg}
     }
  }
  
  def *^ (b : ND):ND = {
     (omat, b) match {
       case (ff:GFilter, bb:GND) => {ff}
     }
  }
  
  def + (bmat:GND):GND = {val (a, b, c, d) = GND.asGMats(amat, bmat, omat, "+"); c ~ a + b; d}
  def - (bmat:GND):GND = {val (a, b, c, d) = GND.asGMats(amat, bmat, omat, "-"); c ~ a - b; d}
  def *@ (bmat:GND):GND = {val (a, b, c, d) = GND.asGMats(amat, bmat, omat, "*@"); c ~ a *@ b; d}
  def ∘ (bmat:GND):GND = {val (a, b, c, d) = GND.asGMats(amat, bmat, omat, "∘"); c ~ a *@ b; d}
  def / (bmat:GND):GND = {val (a, b, c, d) = GND.asGMats(amat, bmat, omat, "/"); c ~ a / b; d}
  def ^ (bmat:GND):GND = {val (a, b, c, d) = GND.asGMats(amat, bmat, omat, "^"); c ~ a ^ b; d}
  
  def > (bmat:GND):GND = {val (a, b, c, d) = GND.asGMats(amat, bmat, omat, ">"); c ~ a > b; d}
  def < (bmat:GND):GND = {val (a, b, c, d) = GND.asGMats(amat, bmat, omat, "<"); c ~ a < b; d}
  def >= (bmat:GND):GND = {val (a, b, c, d) = GND.asGMats(amat, bmat, omat, ">="); c ~ a >= b; d}
  def <= (bmat:GND):GND = {val (a, b, c, d) = GND.asGMats(amat, bmat, omat, "<="); c ~ a <= b; d}
  def != (bmat:GND):GND = {val (a, b, c, d) = GND.asGMats(amat, bmat, omat, "!="); c ~ a != b; d}
  def == (bmat:GND):GND = {val (a, b, c, d) = GND.asGMats(amat, bmat, omat, "=="); c ~ a == b; d}
  def === (bmat:GND):GND = {val (a, b, c, d) = GND.asGMats(amat, bmat, omat, "==="); c ~ a === b; d}
  
  
  def + (b:Float):GND = {val (a, c, d) = GND.asGMats(amat, omat, "+"); c ~ a + GMat(b); d}
  def - (b:Float):GND = {val (a, c, d) = GND.asGMats(amat, omat, "-"); c ~ a - GMat(b); d}
  def *@ (b:Float):GND = {val (a, c, d) = GND.asGMats(amat, omat, "*@"); c ~ a *@ GMat(b); d}
  def ∘ (b:Float):GND = {val (a, c, d) = GND.asGMats(amat, omat, "∘"); c ~ a *@ GMat(b); d}
  def * (b:Float):GND = {val (a, c, d) = GND.asGMats(amat, omat, "*"); c ~ a * GMat(b); d}
  def / (b:Float):GND = {val (a, c, d) = GND.asGMats(amat, omat, "/"); c ~ a / GMat(b); d}  
  def ^ (b:Float):GND = {val (a, c, d) = GND.asGMats(amat, omat, "^"); c ~ a ^ GMat(b); d}
  
  def > (b:Float):GND = {val (a, c, d) = GND.asGMats(amat, omat, ">"); c ~ a > GMat(b); d}
  def < (b:Float):GND = {val (a, c, d) = GND.asGMats(amat, omat, "<"); c ~ a < GMat(b); d}
  def >= (b:Float):GND = {val (a, c, d) = GND.asGMats(amat, omat, ">="); c ~ a >= GMat(b); d}
  def <= (b:Float):GND = {val (a, c, d) = GND.asGMats(amat, omat, "<="); c ~ a <= GMat(b); d}
  def != (b:Float):GND = {val (a, c, d) = GND.asGMats(amat, omat, "!="); c ~ a != GMat(b); d}
  def == (b:Float):GND = {val (a, c, d) = GND.asGMats(amat, omat, "=="); c ~ a == GMat(b); d}  
  def === (b:Float):GND = {val (a, c, d) = GND.asGMats(amat, omat, "==="); c ~ a === GMat(b); d}
  
  
  def + (b:Double):GND = {val (a, c, d) = GND.asGMats(amat, omat, "+"); c ~ a + GMat(b.toFloat); d}
  def - (b:Double):GND = {val (a, c, d) = GND.asGMats(amat, omat, "-"); c ~ a - GMat(b.toFloat); d}
  def *@ (b:Double):GND = {val (a, c, d) = GND.asGMats(amat, omat, "*@"); c ~ a *@ GMat(b.toFloat); d}
  def ∘ (b:Double):GND = {val (a, c, d) = GND.asGMats(amat, omat, "∘"); c ~ a *@ GMat(b.toFloat); d}
  def * (b:Double):GND = {val (a, c, d) = GND.asGMats(amat, omat, "*"); c ~ a * GMat(b.toFloat); d}
  def / (b:Double):GND = {val (a, c, d) = GND.asGMats(amat, omat, "/"); c ~ a / GMat(b.toFloat); d}  
  def ^ (b:Double):GND = {val (a, c, d) = GND.asGMats(amat, omat, "^"); c ~ a ^ GMat(b.toFloat); d}
  
  def > (b:Double):GND = {val (a, c, d) = GND.asGMats(amat, omat, ">"); c ~ a > GMat(b.toFloat); d}
  def < (b:Double):GND = {val (a, c, d) = GND.asGMats(amat, omat, "<"); c ~ a < GMat(b.toFloat); d}
  def >= (b:Double):GND = {val (a, c, d) = GND.asGMats(amat, omat, ">="); c ~ a >= GMat(b.toFloat); d}
  def <= (b:Double):GND = {val (a, c, d) = GND.asGMats(amat, omat, "<="); c ~ a <= GMat(b.toFloat); d}
  def != (b:Double):GND = {val (a, c, d) = GND.asGMats(amat, omat, "!="); c ~ a != GMat(b.toFloat); d}
  def == (b:Double):GND = {val (a, c, d) = GND.asGMats(amat, omat, "=="); c ~ a == GMat(b.toFloat); d}  
  def === (b:Double):GND = {val (a, c, d) = GND.asGMats(amat, omat, "==="); c ~ a === GMat(b.toFloat); d}
  
  
  def + (b:Int):GND = {val (a, c, d) = GND.asGMats(amat, omat, "+"); c ~ a + GMat(b.toFloat); d}
  def - (b:Int):GND = {val (a, c, d) = GND.asGMats(amat, omat, "-"); c ~ a - GMat(b.toFloat); d}
  def *@ (b:Int):GND = {val (a, c, d) = GND.asGMats(amat, omat, "*@"); c ~ a *@ GMat(b.toFloat); d}
  def ∘ (b:Int):GND = {val (a, c, d) = GND.asGMats(amat, omat, "∘"); c ~ a *@ GMat(b.toFloat); d}
  def * (b:Int):GND = {val (a, c, d) = GND.asGMats(amat, omat, "*"); c ~ a * GMat(b.toFloat); d}
  def / (b:Int):GND = {val (a, c, d) = GND.asGMats(amat, omat, "/"); c ~ a / GMat(b.toFloat); d}  
  def ^ (b:Int):GND = {val (a, c, d) = GND.asGMats(amat, omat, "^"); c ~ a ^ GMat(b.toFloat); d}
  
  def > (b:Int):GND = {val (a, c, d) = GND.asGMats(amat, omat, ">"); c ~ a > GMat(b.toFloat); d}
  def < (b:Int):GND = {val (a, c, d) = GND.asGMats(amat, omat, "<"); c ~ a < GMat(b.toFloat); d}
  def >= (b:Int):GND = {val (a, c, d) = GND.asGMats(amat, omat, ">="); c ~ a >= GMat(b.toFloat); d}
  def <= (b:Int):GND = {val (a, c, d) = GND.asGMats(amat, omat, "<="); c ~ a <= GMat(b.toFloat); d}
  def != (b:Int):GND = {val (a, c, d) = GND.asGMats(amat, omat, "!="); c ~ a != GMat(b.toFloat); d}
  def == (b:Int):GND = {val (a, c, d) = GND.asGMats(amat, omat, "=="); c ~ a == GMat(b.toFloat); d}  
  def === (b:Int):GND = {val (a, c, d) = GND.asGMats(amat, omat, "==="); c ~ a === GMat(b.toFloat); d}
  
  
  def + (b:Long):GND = {val (a, c, d) = GND.asGMats(amat, omat, "+"); c ~ a + GMat(b.toFloat); d}
  def - (b:Long):GND = {val (a, c, d) = GND.asGMats(amat, omat, "-"); c ~ a - GMat(b.toFloat); d}
  def *@ (b:Long):GND = {val (a, c, d) = GND.asGMats(amat, omat, "*@"); c ~ a *@ GMat(b.toFloat); d}
  def ∘ (b:Long):GND = {val (a, c, d) = GND.asGMats(amat, omat, "∘"); c ~ a *@ GMat(b.toFloat); d}
  def * (b:Long):GND = {val (a, c, d) = GND.asGMats(amat, omat, "*"); c ~ a * GMat(b.toFloat); d}
  def / (b:Long):GND = {val (a, c, d) = GND.asGMats(amat, omat, "/"); c ~ a / GMat(b.toFloat); d}  
  def ^ (b:Long):GND = {val (a, c, d) = GND.asGMats(amat, omat, "^"); c ~ a ^ GMat(b.toFloat); d}
  
  def > (b:Long):GND = {val (a, c, d) = GND.asGMats(amat, omat, ">"); c ~ a > GMat(b.toFloat); d}
  def < (b:Long):GND = {val (a, c, d) = GND.asGMats(amat, omat, "<"); c ~ a < GMat(b.toFloat); d}
  def >= (b:Long):GND = {val (a, c, d) = GND.asGMats(amat, omat, ">="); c ~ a >= GMat(b.toFloat); d}
  def <= (b:Long):GND = {val (a, c, d) = GND.asGMats(amat, omat, "<="); c ~ a <= GMat(b.toFloat); d}
  def != (b:Long):GND = {val (a, c, d) = GND.asGMats(amat, omat, "!="); c ~ a != GMat(b.toFloat); d}
  def == (b:Long):GND = {val (a, c, d) = GND.asGMats(amat, omat, "=="); c ~ a == GMat(b.toFloat); d}  
  def === (b:Long):GND = {val (a, c, d) = GND.asGMats(amat, omat, "==="); c ~ a === GMat(b.toFloat); d}
  
  
  def * (bmat:ND):GND = this * bmat.asInstanceOf[GND];
  def + (bmat:ND):GND = {val (a, b, c, d) = GND.asGMats(amat, bmat.asInstanceOf[GND], omat, "+"); c ~ a + b; d}
  def - (bmat:ND):GND = {val (a, b, c, d) = GND.asGMats(amat, bmat.asInstanceOf[GND], omat, "-"); c ~ a - b; d}
  def *@ (bmat:ND):GND = {val (a, b, c, d) = GND.asGMats(amat, bmat.asInstanceOf[GND], omat, "*@"); c ~ a *@ b; d}
  def ∘ (bmat:ND):GND = {val (a, b, c, d) = GND.asGMats(amat, bmat.asInstanceOf[GND], omat, "∘"); c ~ a *@ b; d}
  def / (bmat:ND):GND = {val (a, b, c, d) = GND.asGMats(amat, bmat.asInstanceOf[GND], omat, "/"); c ~ a / b; d}
  def ^ (bmat:ND):GND = {val (a, b, c, d) = GND.asGMats(amat, bmat.asInstanceOf[GND], omat, "^"); c ~ a ^ b; d}
  
  def > (bmat:ND):GND = {val (a, b, c, d) = GND.asGMats(amat, bmat.asInstanceOf[GND], omat, ">"); c ~ a > b; d}
  def < (bmat:ND):GND = {val (a, b, c, d) = GND.asGMats(amat, bmat.asInstanceOf[GND], omat, "<"); c ~ a < b; d}
  def >= (bmat:ND):GND = {val (a, b, c, d) = GND.asGMats(amat, bmat.asInstanceOf[GND], omat, ">="); c ~ a >= b; d}
  def <= (bmat:ND):GND = {val (a, b, c, d) = GND.asGMats(amat, bmat.asInstanceOf[GND], omat, "<="); c ~ a <= b; d}
  def != (bmat:ND):GND = {val (a, b, c, d) = GND.asGMats(amat, bmat.asInstanceOf[GND], omat, "!="); c ~ a != b; d}
  def == (bmat:ND):GND = {val (a, b, c, d) = GND.asGMats(amat, bmat.asInstanceOf[GND], omat, "=="); c ~ a == b; d}
  def === (bmat:ND):GND = {val (a, b, c, d) = GND.asGMats(amat, bmat.asInstanceOf[GND], omat, "==="); c ~ a === b; d}
}

object GND {
  
  def scalar(v:Float, nd:Int):GND = {
    val newdims = new Array[Int](nd);
    Arrays.fill(newdims,1);
    val out = GND(newdims);
    val tmp = Array[Float](v);
    CPUtoGPUarraycopy(tmp, 0, out.data, 0, 1, "GND scalar");
    out
  }
  
  def apply(dims:Array[Int]):GND = {
    val retv = new GND(dims, new Pointer);
    val len = dims.reduce(_*_);
    if (Mat.debugMem) {
      println("GND %d, %d %f" format (len, SciFunctions.getGPU, SciFunctions.GPUmem._1))
      if (len > Mat.debugMemThreshold) throw new RuntimeException("GND alloc too large");
    }
    var err = if (1L*len*Sizeof.FLOAT > Mat.hostAllocSize) {
      cudaMallocHost(retv.data, 1L*len*Sizeof.FLOAT);
    } else {
      cudaMalloc(retv.data, 1L*len*Sizeof.FLOAT);
    }
    cudaDeviceSynchronize;
    if (err == 0) err = cudaGetLastError();
    if (err != 0) throw new RuntimeException("CUDA alloc failed " + cudaGetErrorString(err));
    retv       
  }
  
  def apply(dims0:IMat):GND = apply(dims0.data);
  
  def apply(dims:Int*):GND = apply(dims.toArray);
  
  def apply(f:GMat):GND = {
    val out = GND.newOrCheckGND(Array(f.nrows, f.ncols), null, f.GUID, "apply".##);
    GPUtoGPUarraycopy(f.data, 0, out.data, 0, f.length, "GND apply");
    out
  }
  
  def apply(f:FND):GND = {
    val out = GND.newOrCheckGND(f.dims.data, null, f.GUID, "apply".##);
    CPUtoGPUarraycopy(f.data, 0, out.data, 0, f.length, "GND apply");
    out
  }
  
  def apply(f:FMat):GND = {
    val out = GND.newOrCheckGND(Array(f.nrows, f.ncols), null, f.GUID, "apply".##);
    CPUtoGPUarraycopy(f.data, 0, out.data, 0, f.length, "GND apply");
    out
  }
  
  def zeros(dims0:IMat) = {
    val out = GND(dims0)
    cudaMemset(out.data, 0, Sizeof.FLOAT*out.length)
    cudaDeviceSynchronize()
    val err = cudaGetLastError()
    if (err != 0) throw new RuntimeException("GPU "+SciFunctions.getGPU+": Cuda error in GND.zeros " + cudaGetErrorString(err));
    out
  }
  
  def ones(dims0:IMat) = {
    val out = GND(dims0)
    CUMAT.setval(out.data, 1f, out.length)
    cudaDeviceSynchronize()
    val err = cudaGetLastError()
    if (err != 0) throw new RuntimeException("GPU "+SciFunctions.getGPU+": Cuda error in GND.ones " + cudaGetErrorString(err));
    out
  }
  
  def applyGNDfun(in:GND, omat:ND, opn:Int, kflops:Long):GND = {
    val out = GND.newOrCheckGND(in.dims, omat, in.GUID, opn)
    CUMAT.applygfun(in.data, out.data, in.length, opn)
    jcuda.runtime.JCuda.cudaDeviceSynchronize()
    Mat.nflops += kflops*in.length
    out
  }

  def applyGNDfun(in:GND, opn:Int, kflops:Long):GND = {
    val out = GND.newOrCheckGND(in.dims, null, in.GUID, opn)
    CUMAT.applygfun(in.data, out.data, in.length, opn)
    jcuda.runtime.JCuda.cudaDeviceSynchronize()
    Mat.nflops += kflops*in.length
    out
  }
  
  def applyGNDfun2(a:GND, b:GND, omat:ND, opn:Int, kflops:Long):GND = {   
  	ND.checkDims("applyGNDfun2", a.dims, b.dims);
  	val out = GND.newOrCheckGND(a.dims, omat, a.GUID, b.GUID, opn);
  	CUMAT.applygfun2(a.data, b.data, out.data, a.nrows*a.ncols, opn);
  	jcuda.runtime.JCuda.cudaDeviceSynchronize();
  	Mat.nflops += kflops*a.length;
  	out;
  }
  
  def applyGNDfun2(a:GND, b:GND, opn:Int, kflops:Long):GND = {
    ND.checkDims("applyGNDfun2", a.dims, b.dims);
    val out = GND.newOrCheckGND(a.dims, null, a.GUID, b.GUID, opn);
    CUMAT.applygfun2(a.data, b.data, out.data, a.nrows*a.ncols, opn);
    jcuda.runtime.JCuda.cudaDeviceSynchronize();
    Mat.nflops += kflops*a.length;
    out;
  }
  
  def rand(out:GND):GND = {
    GMat.rand(out.asMat);
    out;
  }
   
  def normrnd(mu:Float, sig:Float, out:GND):GND = {
    GMat.normrnd(mu, sig, out.asMat);
    out;
  }
  
  def gamrnd(a:GND, b:GND, out:GND):GND = { 
    GMat.gamrnd(a.asMat, b.asMat, out.asMat);
    out;
  }

  
  def asGMats(mat1:GND, mat2:GND, omat:ND, opname:String):(GMat, GMat, GMat, GND) = {
    if (mat1._dims.length != mat2._dims.length) {
      throw new RuntimeException("Operator "+opname+" inconsistent number of dims in operands")
    }
    val (nr1, nc1, nr2, nc2) = ND.compatibleDims(mat1._dims, mat2._dims, opname);
    val a = new GMat(nr1, nc1, mat1.data, nr1*nc1);
    val b = new GMat(nr2, nc2, mat2.data, nr2*nc2);
    val nr3 = math.max(nr1, nr2)
    val nc3 = math.max(nc1, nc2)
    val xdims = if (mat1.length > mat2.length) mat1._dims else mat2._dims
    val d = GND.newOrCheckGND(xdims, omat, mat1.GUID, mat2.GUID, opname.##)
    val c = new GMat(nr3, nc3, d.data, nr3*nc3);
    (a, b, c, d)
  }
  
  def asGMats(mat1:GND, omat:ND, opname:String):(GMat, GMat, GND) = {
    val d = GND.newOrCheckGND(mat1._dims, omat, mat1.GUID, opname.##)
    val a = new GMat(mat1.length, 1, mat1.data, mat1.length);
    val c = new GMat(mat1.length, 1, d.data, mat1.length);
    (a, c, d)
  }
  
  def newOrCheckGND(dims:Array[Int], out:ND):GND = {
    if (out.asInstanceOf[AnyRef] != null && ND.checkDims("GND newOrCheckGND: ", out.dims.data, dims)) {
    	out.asInstanceOf[GND]
    } else {
      GND(dims)
    }
  }
  
  def newOrCheckGND(dims:IMat, out:ND):GND = newOrCheckGND(dims.data, out);
    
  def newOrCheckGND(dims:Array[Int], out:ND, matGuid:Long, opHash:Int):GND = {
    if (out.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
      newOrCheckGND(dims, out)
    } else {
      val key = (matGuid, opHash)
      val res = ND.cache2(key)
      if (res != null) {
      	newOrCheckGND(dims, res)
      } else {
        val omat = newOrCheckGND(dims, null)
        ND.cache2put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckGND(dims:IMat, out:ND, g1:Long, opHash:Int):GND = newOrCheckGND(dims.data, out, g1, opHash);
  
  def newOrCheckGND(dims:Array[Int], out:ND, guid1:Long, guid2:Long, opHash:Int):GND = {
    if (out.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
      newOrCheckGND(dims, out)
    } else {
      val key = (guid1, guid2, opHash)
      val res = ND.cache3(key)
      if (res != null) {
      	newOrCheckGND(dims, res)
      } else {
        val omat = newOrCheckGND(dims, null)
        ND.cache3put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckGND(dims:IMat, out:ND, g1:Long, g2:Long, opHash:Int):GND = newOrCheckGND(dims.data, out, g1, g2, opHash);
  
  def newOrCheckGND(dims:Array[Int], out:ND, g1:Long, g2:Long, g3:Long, opHash:Int):GND = {
    if (out.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
      newOrCheckGND(dims, out)
    } else {
      val key = (g1, g2, g3, opHash)
      val res = ND.cache4(key)
      if (res != null) {
      	newOrCheckGND(dims, res)
      } else {
        val omat = newOrCheckGND(dims, null)
        ND.cache4put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckGND(dims:IMat, out:ND, g1:Long, g2:Long, g3:Long, opHash:Int):GND = newOrCheckGND(dims.data, out, g1, g2, g3, opHash);
}






