package BIDMat
import scala.collection.mutable.HashMap;
import java.lang.ref._;
import jcuda.NativePointerObject;
import edu.berkeley.bid.UTILS;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.Handler;
import java.util.logging.FileHandler;
import java.util.logging.SimpleFormatter;

@SerialVersionUID(100L)
class Mat(val _dims:Array[Int]) extends ND with Serializable {
  
  def this(nr:Int, nc:Int) = this(Array(nr, nc));
  
  def mytype = "Mat";
  
  var _idims:IMat = null
   
  def dims = {
    if (_idims.asInstanceOf[IMat] == null) {
      _idims = new IMat(1, _dims.length, _dims);
    }
    _idims;
  }
  
  val ndims = _dims.length;
  val ncols = _dims(ndims-1)
  val nrows = if (ndims == 1) 1 else _dims.slice(0,ndims-1).reduce(_*_);
  val length = nrows * ncols; 
  val llength = 1L*nrows*ncols;
  def nnz = length;
 
  private var _GUID = Mat.myrand.nextLong
  def GUID = _GUID;
  def setGUID(v:Long) {
    _GUID = v;
  }
  
  def dv:Double = {
    throw new RuntimeException("dv not implemented for "+this.mytype)
  }
  
  def fv:Float = {
    throw new RuntimeException("fv not implemented for "+this.mytype)
  }
  
  def notImplemented0(s:String):Mat = { 
    throw new RuntimeException(s+" not implemented for "+this.mytype)
  } 
  
  def notImplemented1(s:String,that:Float):Mat = { 
    throw new RuntimeException(s+" not implemented for "+this.mytype+" and Float")
  }
  
  def notImplemented1(s:String,that:Int):Mat = { 
    throw new RuntimeException(s+" not implemented for "+this.mytype+" and Int")
  }
  
  def notImplemented1(s:String,that:Double):Mat = { 
    throw new RuntimeException(s+" not implemented for "+this.mytype+" and Double")
  }
  
  def notImplemented1(s:String,that:Mat):Mat = { 
    throw new RuntimeException(s+" not implemented for "+this.mytype+" and Double")
  }
    
  def t = notImplemented0("t")  
  
  def apply(i1:IMat):Mat = notImplemented0("1D slice");
  def apply(i1:IMat, i2:IMat):Mat = notImplemented0("2D slice");
  def apply(i1:IMat, i2:IMat, i3:IMat):Mat = notImplemented0("3D slice");
  def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat):Mat = notImplemented0("4D slice");
  def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat):Mat = notImplemented0("5D slice");
  def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat):Mat = notImplemented0("6D slice");
  def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat):Mat = notImplemented0("7D slice");
  def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat):Mat = notImplemented0("8D slice");
  def apply(i1:IMat, i2:Int):Mat = notImplemented0("2D slice");
  def apply(i1:Int, i2:IMat):Mat = notImplemented0("2D slice");
    
  def update(i1:IMat, v:Mat):Mat = notImplemented0("1D update");
  def update(i1:IMat, i2:IMat, v:Mat):Mat = notImplemented0("2D update");
  def update(i1:IMat, i2:IMat, i3:IMat, v:Mat):Mat = notImplemented0("3D update");
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, v:Mat):Mat = notImplemented0("4D update");
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, v:Mat):Mat = notImplemented0("5D update");
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, v:Mat):Mat = notImplemented0("6D update");
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, v:Mat):Mat = notImplemented0("7D update");
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat, v:Mat):Mat = notImplemented0("8D update");
  def update(i1:IMat, i2:Int, v:Mat):Mat = notImplemented0("2D update");
  def update(i1:Int, i2:IMat, v:Mat):Mat = notImplemented0("2D update");
  
  def update(i1:IMat, v:Float):Mat = notImplemented0("1D update");
  def update(i1:IMat, i2:IMat, v:Float):Mat = notImplemented0("2D update");
  def update(i1:IMat, i2:IMat, i3:IMat, v:Float):Mat = notImplemented0("3D update");
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, v:Float):Mat = notImplemented0("4D update");
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, v:Float):Mat = notImplemented0("5D update");
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, v:Float):Mat = notImplemented0("6D update");
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, v:Float):Mat = notImplemented0("7D update");
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat, v:Float):Mat = notImplemented0("8D update");  
  
  def update(i1:Int, v:Float):Mat = notImplemented0("1D update");
  def update(i1:Int, i2:Int, v:Float):Mat = notImplemented0("2D update");
  def update(i1:Int, i2:IMat, v:Float):Mat = notImplemented0("2D update");
  def update(i1:IMat, i2:Int, v:Float):Mat = notImplemented0("2D update");
  
  def update(i1:IMat, v:Double):Mat = notImplemented0("1D update");
  def update(i1:IMat, i2:IMat, v:Double):Mat = notImplemented0("2D update");
  def update(i1:IMat, i2:IMat, i3:IMat, v:Double):Mat = notImplemented0("3D update");
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, v:Double):Mat = notImplemented0("4D update");
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, v:Double):Mat = notImplemented0("5D update");
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, v:Double):Mat = notImplemented0("6D update");
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, v:Double):Mat = notImplemented0("7D update");
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat, v:Double):Mat = notImplemented0("8D update");  
  
  def update(i1:Int, v:Double):Mat = notImplemented0("1D update");
  def update(i1:Int, i2:Int, v:Double):Mat = notImplemented0("2D update");
  def update(i1:Int, i2:IMat, v:Double):Mat = notImplemented0("2D update");
  def update(i1:IMat, i2:Int, v:Double):Mat = notImplemented0("2D update");
  
  def update(i1:IMat, v:Int):Mat = notImplemented0("1D update");
  def update(i1:IMat, i2:IMat, v:Int):Mat = notImplemented0("2D update");
  def update(i1:IMat, i2:IMat, i3:IMat, v:Int):Mat = notImplemented0("3D update");
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, v:Int):Mat = notImplemented0("4D update");
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, v:Int):Mat = notImplemented0("5D update");
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, v:Int):Mat = notImplemented0("6D update");
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, v:Int):Mat = notImplemented0("7D update");
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat, v:Int):Mat = notImplemented0("8D update");  
  
  def update(i1:Int, v:Int):Mat = notImplemented0("1D update");
  def update(i1:Int, i2:Int, v:Int):Mat = notImplemented0("2D update");
  def update(i1:Int, i2:IMat, v:Int):Mat = notImplemented0("2D update");
  def update(i1:IMat, i2:Int, v:Int):Mat = notImplemented0("2D update");
  
  def update(i1:IMat, v:Long):Mat = notImplemented0("1D update");
  def update(i1:IMat, i2:IMat, v:Long):Mat = notImplemented0("2D update");
  def update(i1:IMat, i2:IMat, i3:IMat, v:Long):Mat = notImplemented0("3D update");
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, v:Long):Mat = notImplemented0("4D update");
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, v:Long):Mat = notImplemented0("5D update");
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, v:Long):Mat = notImplemented0("6D update");
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, v:Long):Mat = notImplemented0("7D update");
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat, v:Long):Mat = notImplemented0("8D update");  
  
  def update(i1:Int, v:Long):Mat = notImplemented0("1D update");
  def update(i1:Int, i2:Int, v:Long):Mat = notImplemented0("2D update");
  def update(i1:Int, i2:IMat, v:Long):Mat = notImplemented0("2D update");
  def update(i1:IMat, i2:Int, v:Long):Mat = notImplemented0("2D update"); 
  
  def copyTo(a:Mat) = notImplemented0("copyTo");
  def copy = notImplemented0("copy");
  def newcopy = notImplemented0("newcopy");
  def set(v:Float) = notImplemented0("set");
  def set(v:Double) = notImplemented0("set");

  def zeros(nr:Int, nc:Int) = notImplemented0("zeros");
  def zeros(dims:IMat) = notImplemented0("zeros");
  def zero = notImplemented0("zero");
  def zerosLike = notImplemented0("zerosLike");
  def zeros(nr:Int, nc:Int, nnz:Int):Mat = zeros(nr, nc);
   
  def ones(nr:Int, nc:Int) = notImplemented0("ones");
  def ones(dims:IMat) = notImplemented0("ones");
  def one = notImplemented0("one");
  def onesLike = notImplemented0("onesLike");
  def ones(nr:Int, nc:Int, nnz:Int):Mat = ones(nr, nc);

  def izeros(nr:Int, nc:Int) = notImplemented0("izeros");
  def izeros(dims:IMat) = notImplemented0("izeros");
  def iones(nr:Int, nc:Int) = notImplemented0("iones");
  def iones(dims:IMat) = notImplemented0("iones");
  
  def clear = notImplemented0("clear");
  def clearUpper(i:Int) = notImplemented0("clearUpper");
  def clearLower(i:Int) = notImplemented0("clearLower"); 
  def clearUpper = notImplemented0("clearUpper");
  def clearLower = notImplemented0("clearLower");
  def free = notImplemented0("free");
  def view(nr:Int, nc:Int):Mat = notImplemented0("view");
  def view(nr:Int, nc:Int, setGUID:Boolean):Mat = notImplemented0("view");
  
  def recycle(nr:Int, nc:Int, nnz:Int):Mat = notImplemented0("recycle");
  def contents:Mat = notImplemented0("contents");
  
  def colslice(a:Int, b:Int):Mat = notImplemented0("colslice");
  def colslice(a:Int, b:Int, m:Mat):Mat = notImplemented0("colslice");
  def colslice(a:Int, b:Int, m:Mat, c:Int):Mat = notImplemented0("colslice");
  def colslice(a:Int, b:Int, m:Mat, c:Int, pb:Boolean):Mat = notImplemented0("colslice");

  def rowslice(a:Int, b:Int, out:Mat):Mat = notImplemented0("rowslice");
  def rowslice(a:Int, b:Int, out:Mat, c:Int):Mat = notImplemented0("rowslice");
  def rowslice(a:Int, b:Int):Mat = notImplemented0("rowslice");
  
  def reshape(newdims:Int*):Mat = notImplemented0("reshape");
  def reshape(newdims:Array[Int]):Mat = notImplemented0("reshape");
  def reshape(newdims:IMat):Mat = notImplemented0("reshape");

  def reshapeView(newdims:Int*):Mat  = notImplemented0("reshapeView");
  def reshapeView(newdims:Array[Int]):Mat = notImplemented0("reshapeView");
  def reshapeView(newdims:IMat):Mat = notImplemented0("reshapeView");
  
  def transpose(dims:Array[Int]):Mat = notImplemented0("transpose");
  def transpose(dims:IMat):Mat = notImplemented0("transpose");
  def transpose(d1:Int, d2:Int):Mat = notImplemented0("transpose");
  def transpose(d1:Int, d2:Int, d3:Int):Mat = notImplemented0("transpose");
  def transpose(d1:Int, d2:Int, d3:Int, d4:Int):Mat = notImplemented0("transpose");
  def transpose(d1:Int, d2:Int, d3:Int, d4:Int, d5:Int):Mat = notImplemented0("transpose");
  def transpose(d1:Int, d2:Int, d3:Int, d4:Int, d5:Int, d6:Int):Mat = notImplemented0("transpose"); 
  def transpose(d1:Int, d2:Int, d3:Int, d4:Int, d5:Int, d6:Int, d7:Int):Mat = notImplemented0("transpose"); 
  def transpose(d1:Int, d2:Int, d3:Int, d4:Int, d5:Int, d6:Int, d7:Int, d8:Int):Mat = notImplemented0("transpose"); 
   
  
  def tileMult(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:Mat, broff:Int, bcoff:Int, c:Mat, croff:Int, ccoff:Int):Mat =
    notImplemented0("tileMult");
  def tileMultNT(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:Mat, broff:Int, bcoff:Int, c:Mat, croff:Int, ccoff:Int):Mat =
    notImplemented0("tileMultT");
  def tileMultTN(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:Mat, broff:Int, bcoff:Int, c:Mat, croff:Int, ccoff:Int):Mat =
  		notImplemented0("tileMultT");
  def tileMultTT(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:Mat, broff:Int, bcoff:Int, c:Mat, croff:Int, ccoff:Int):Mat =
    notImplemented0("tileMultT");
  
  def tileCopy(fromrow:Int, fromcol:Int, to:Mat, torow:Int, tocol:Int, height:Int, width:Int):Mat = {
  		notImplemented0("tileCopy");
  }
  
  def vecAdd(fromi:Int, to:Mat, toi:Int, height:Int):Mat = {
  		notImplemented0("vecAdd");
  }

  def blockGemm(transa:Int, transb:Int, nr:Int, nc:Int, k:Int, reps:Int, aoff:Int, lda:Int, astep:Int, 
      b:Mat, boff:Int, ldb:Int, bstep:Int, c:Mat, coff:Int, ldc:Int, cstep:Int, addC:Boolean):Mat = notImplemented0("blockGemm");
  
  def mult(a:Mat, b:Mat, at:Boolean, bt:Boolean):Mat= notImplemented1("mult", a);
  def mult(a:Mat, b:Mat):Mat= notImplemented1("mult", a);

  def madd(a:Mat, b:Mat, at:Boolean, bt:Boolean):Mat = notImplemented1("madd", a);
  def madd(a:Mat, b:Mat):Mat = notImplemented1("madd", a);
  
  def blockmult(a:Mat, b:Mat, ngroups:Int, at:Boolean, bt:Boolean):Mat= notImplemented1("groupmult", a);
  def blockmult(a:Mat, b:Mat, ngroups:Int):Mat= notImplemented1("groupmult", a);

  def blockmadd(a:Mat, b:Mat, ngroups:Int, at:Boolean, bt:Boolean):Mat = notImplemented1("groupmadd", a);
  def blockmadd(a:Mat, b:Mat, ngroups:Int):Mat = notImplemented1("groupmadd", a);
  
  def unary_-():Mat = notImplemented1("-", this)
  
  def <-- (b : Mat):Mat = b.copyTo(this)
  
  def +  (b : Mat):Mat = notImplemented1("+", b)
  def -  (b : Mat):Mat = notImplemented1("-", b)
  def *  (b : Mat):Mat = notImplemented1("*", b)
  def *^ (b : Mat):Mat = notImplemented1("*^", b)
  def xT (b : Mat):Mat = notImplemented1("*", b)
  def Tx (b : Mat):Mat = notImplemented1("*", b)
  def ^* (b : Mat):Mat = notImplemented1("^*", b)
  def ** (b : Mat):Mat = notImplemented1("**", b)
  def ⊗  (b : Mat):Mat = notImplemented1("⊗", b)       // unicode 8855, 0x2297
  def /< (b : Mat):Mat = notImplemented1("/<", b)
  def ∘  (b : Mat):Mat = notImplemented1("∘", b)        // unicode 8728, 0x2218 
  def *@ (b : Mat):Mat = notImplemented1("*@", b)
  def /  (b : Mat):Mat = notImplemented1("/", b)
  def \\ (b : Mat):Mat = notImplemented1("\\\\", b)
  def ^  (b : Mat):Mat = notImplemented1("^", b) 
  def ◁  (b : Mat):Mat = notImplemented1("◁", b)        // unicode 9665, 0x25C1 
  def ▷  (b : Mat):Mat = notImplemented1("▷", b)        // unicode 9666, 0x25C2
  def dot (b : Mat):Mat = notImplemented1("dot", b)
  def dotr (b : Mat):Mat = notImplemented1("dotr", b) 
  def ∙ (b : Mat):Mat = notImplemented1("dot", b)       // unicode 8729, 0x2219 
  def ∙→ (b : Mat):Mat = notImplemented1("dotr", b)     // unicode (8729, 8594) (0x2219, 0x2192)
    
  def >  (b : Mat):Mat = notImplemented1(">", b)
  def <  (b : Mat):Mat = notImplemented1("<", b)
  def >= (b : Mat):Mat = notImplemented1(">=", b)
  def <= (b : Mat):Mat = notImplemented1("<=", b)
  def == (b : Mat):Mat = notImplemented1("==", b)
  def === (b : Mat):Mat = notImplemented1("===", b)
  def != (b : Mat):Mat = notImplemented1("!=", b)
  
  def \ (b : Mat):Mat = notImplemented1("\\", b)
  def on (b : Mat):Mat = notImplemented1("on", b)
  
  def max(b : Mat):Mat = notImplemented1("max", b)
  def min(b : Mat):Mat = notImplemented1("min", b)
  
   def sum(b : IMat):Mat = notImplemented0("sum")
  def prod(b : IMat):Mat = notImplemented0("prod")
  def maxi(b : IMat):Mat = notImplemented0("maxi")
  def mini(b : IMat):Mat = notImplemented0("mini")
  def amax(b : IMat):Mat = notImplemented0("maxi")
  def amin(b : IMat):Mat = notImplemented0("mini")
  def mean(b : IMat):Mat = notImplemented0("mean")
  def variance(b : IMat):Mat = notImplemented0("variance")
  
  def sum(i1:Int):Mat = notImplemented0("sum");
  def sum(i1:Int, i2:Int):Mat = notImplemented0("sum");
  def sum(i1:Int, i2:Int, i3:Int):Mat = notImplemented0("sum");  
  def sum(i1:Int, i2:Int, i3:Int, i4:Int):Mat = notImplemented0("sum");   
  def sum(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int):Mat = notImplemented0("sum");   
  def sum(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int):Mat = notImplemented0("sum");  
  
  def prod(i1:Int):Mat = notImplemented0("prod");
  def prod(i1:Int, i2:Int):Mat = notImplemented0("prod");
  def prod(i1:Int, i2:Int, i3:Int):Mat = notImplemented0("prod");  
  def prod(i1:Int, i2:Int, i3:Int, i4:Int):Mat = notImplemented0("prod");   
  def prod(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int):Mat = notImplemented0("prod");   
  def prod(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int):Mat = notImplemented0("prod");  
  
  def maxi(i1:Int):Mat = notImplemented0("maxi");
  def amax(i1:Int):Mat = notImplemented0("amax"); 
    
  def mini(i1:Int):Mat = notImplemented0("mini");
  def amin(i1:Int):Mat = notImplemented0("amin");
  
  def mean(i1:Int):Mat = notImplemented0("mean");  
  
  def variance(i1:Int):Mat = notImplemented0("variance");
  
  def *  (b : Float):Mat = notImplemented1("*", b)
  def +  (b : Float):Mat = notImplemented1("+", b)
  def -  (b : Float):Mat = notImplemented1("-", b)
  def *@ (b : Float):Mat = notImplemented1("*@", b)
  def ∘  (b : Float):Mat = notImplemented1("∘", b)
  def / (b : Float):Mat = notImplemented1("/", b)
  def ^  (b : Float):Mat = notImplemented1("^", b)
  
  def >  (b : Float):Mat = notImplemented1(">", b)
  def <  (b : Float):Mat = notImplemented1("<", b)
  def >= (b : Float):Mat = notImplemented1(">=", b)
  def <= (b : Float):Mat = notImplemented1("<=", b)
  def == (b : Float):Mat = notImplemented1("==", b)
  def === (b : Float):Mat = notImplemented1("===", b)
  def != (b : Float):Mat = notImplemented1("!=", b)
  
  def max(b : Float):Mat = notImplemented1("max", b)
  def min(b : Float):Mat = notImplemented1("min", b)
  
  def *  (b : Int):Mat = notImplemented1("*", b)
  def +  (b : Int):Mat = notImplemented1("+", b)
  def -  (b : Int):Mat = notImplemented1("-", b)
  def *@ (b : Int):Mat = notImplemented1("*@", b)
  def ∘  (b : Int):Mat = notImplemented1("∘", b)
  def / (b : Int):Mat = notImplemented1("/", b)
  def ^  (b : Int):Mat = notImplemented1("^", b)
  
  def >  (b : Int):Mat = notImplemented1(">", b)
  def <  (b : Int):Mat = notImplemented1("<", b)
  def >= (b : Int):Mat = notImplemented1(">=", b)
  def <= (b : Int):Mat = notImplemented1("<=", b)
  def == (b : Int):Mat = notImplemented1("==", b)
  def === (b : Int):Mat = notImplemented1("===", b)
  def != (b : Int):Mat = notImplemented1("!=", b)
  
  def max(b : Int):Mat = notImplemented1("max", b)
  def min(b : Int):Mat = notImplemented1("min", b)
  
  def + (b : Long):Mat = notImplemented1("+", b)
  def - (b : Long):Mat = notImplemented1("-", b)
  def * (b : Long):Mat = notImplemented1("*", b)
  def *@ (b : Long):Mat = notImplemented1("*@", b)
  def ∘  (b : Long):Mat = notImplemented1("∘", b)
  def /  (b : Long):Mat = notImplemented1("/", b)
  def ^ (b : Long):Mat   = notImplemented1("^", b)
  
  def > (b : Long):Mat = notImplemented1(">", b)
  def < (b : Long):Mat = notImplemented1("<", b)
  def >= (b : Long):Mat = notImplemented1(">=", b)
  def <= (b : Long):Mat = notImplemented1("<=", b)
  def == (b : Long):Mat = notImplemented1("==", b)
  def === (b : Long):Mat = notImplemented1("===", b)
  def != (b : Long):Mat = notImplemented1("!=", b)
  
  def max(b : Long):Mat = notImplemented1("max", b)
  def min(b : Long):Mat = notImplemented1("min", b)
  
  def *  (b : Double):Mat = notImplemented1("*", b)
  def +  (b : Double):Mat = notImplemented1("+", b)
  def -  (b : Double):Mat = notImplemented1("-", b)
  def *@ (b : Double):Mat = notImplemented1("*@", b)
  def ∘  (b : Double):Mat = notImplemented1("∘", b)
  def / (b : Double):Mat = notImplemented1("/", b)
  def ^  (b : Double):Mat = notImplemented1("^", b)
  
  def >  (b : Double):Mat = notImplemented1(">", b)
  def <  (b : Double):Mat = notImplemented1("<", b)
  def >= (b : Double):Mat = notImplemented1(">=", b)
  def <= (b : Double):Mat = notImplemented1("<=", b)
  def == (b : Double):Mat = notImplemented1("==", b)
  def === (b : Double):Mat = notImplemented1("===", b)
  def != (b : Double):Mat = notImplemented1("!=", b)
  
  def max(b : Double):Mat = notImplemented1("max", b)
  def min(b : Double):Mat = notImplemented1("min", b)
  
  
  def ddot (b : Mat):Double = {notImplemented1("ddot", b); 0}
  def ∙∙ (b : Mat):Double = {notImplemented1("ddot", b); 0}
  
  
  def ^* (b : DSPair):Mat = notImplemented0("^*")
  def Tx (b : DSPair):Mat = notImplemented0("Tx")
  def @@ (b : Mat):DSPair = (this, b) match {
    case (aa:GMat, bb:GSMat) => new GDSPair(aa, bb)
    case (aa:FMat, bb:SMat) => new FDSPair(aa, bb) 
  }
  
  def ~ (b : Mat):Pair = b match {
  case bb:GMat => new GPair(this, bb)
  case bb:GIMat => new GIPair(this, bb)
  case bb:GDMat => new GDPair(this, bb)
  case bb:GLMat => new GLPair(this, bb)
  case bb:FMat => new FPair(this, bb)
  case bb:DMat => new DPair(this, bb)
  case bb:IMat => new IPair(this, bb)
  case bb:SMat => new SPair(this, bb)
  case bb:SDMat => new SDPair(this, bb)
  case bb:CMat => new CPair(this, bb)

  }
  
}

abstract class DSPair extends Serializable {}

abstract class Pair(omat:Mat, mat:Mat) extends Serializable {
  
  def notImplemented0(s:String):Mat = { 
    throw new RuntimeException("operator "+s+" not implemented for "+this)
  }
  
  def t = notImplemented0("t")
  
  def notImplemented1(s:String,that:Mat):Mat = { 
    throw new RuntimeException(s+" not implemented for "+mat.mytype+" and Int")
  }
  
  def +  (b : Mat):Mat = notImplemented1("+", b)
  def -  (b : Mat):Mat = notImplemented1("-", b)
  def *  (b : Mat):Mat = notImplemented1("*", b)
  def *^ (b : Mat):Mat = notImplemented1("*^", b)
  def xT (b : Mat):Mat = notImplemented1("*", b)
  def Tx (b : Mat):Mat = notImplemented1("*", b)
  def ^* (b : Mat):Mat = notImplemented1("^*", b)
  def ** (b : Mat):Mat = notImplemented1("**", b)
  def ⊗  (b : Mat):Mat = notImplemented1("⊗", b)       // unicode 8855, 0x2297
  def /< (b : Mat):Mat = notImplemented1("/<", b)
  def ∘  (b : Mat):Mat = notImplemented1("∘", b)        // unicode 8728, 0x2218 
  def *@ (b : Mat):Mat = notImplemented1("*@", b)
  def /  (b : Mat):Mat = notImplemented1("/", b)
  def \\ (b : Mat):Mat = notImplemented1("\\\\", b)
  def ^  (b : Mat):Mat = notImplemented1("^", b) 
  def ◁  (b : Mat):Mat = notImplemented1("◁", b)        // unicode 9665, 0x25C1 
  def ▷  (b : Mat):Mat = notImplemented1("▷", b)        // unicode 9666, 0x25C2
  def dot (b : Mat):Mat = notImplemented1("dot", b)
  def dotr (b : Mat):Mat = notImplemented1("dotr", b) 
  def ∙ (b : Mat):Mat = notImplemented1("dot", b)       // unicode 8729, 0x2219 
  def ∙→ (b : Mat):Mat = notImplemented1("dotr", b)     // unicode (8729, 8594) (0x2219, 0x2192)
    
  def >  (b : Mat):Mat = notImplemented1(">", b)
  def <  (b : Mat):Mat = notImplemented1("<", b)
  def >= (b : Mat):Mat = notImplemented1(">=", b)
  def <= (b : Mat):Mat = notImplemented1("<=", b)
  def == (b : Mat):Mat = notImplemented1("==", b)
  def === (b : Mat):Mat = notImplemented1("===", b)
  def != (b : Mat):Mat = notImplemented1("!=", b)
  
  def max(b : Mat):Mat = notImplemented1("max", b)
  def min(b : Mat):Mat = notImplemented1("min", b)
  
  def \ (b : Mat):Mat = notImplemented1("max", b)
  def on(b : Mat):Mat = notImplemented1("min", b)
  
  def + (b : Float):Mat = notImplemented0("+")
  def - (b : Float):Mat = notImplemented0("-")
  def * (b : Float):Mat = notImplemented0("*")
  def xT (b : Float):Mat = notImplemented0("xT")
  def *^ (b : Float):Mat = notImplemented0("*^")
  def Tx (b : Float):Mat = notImplemented0("Tx")
  def ^* (b : Float):Mat = notImplemented0("*^")
  def /< (b : Float):Mat = notImplemented0("/<")
  def *@ (b : Float):Mat = notImplemented0("*@")
  def ∘  (b : Float):Mat = notImplemented0("∘")
  def /  (b : Float):Mat = notImplemented0("/")
  def \\ (b : Float):Mat = notImplemented0("\\\\")
  def ^ (b : Float):Mat = notImplemented0("^") 
  def ◁ (b : Float):Mat = notImplemented0("◁")
  def ▷ (b : Float):Mat = notImplemented0("▷")
  def dot (b : Float):Mat = notImplemented0("dot")
  def dotr (b : Float):Mat = notImplemented0("dotr")
  def ∙ (b : Float):Mat = notImplemented0("dot")
  def ∙→ (b : Float):Mat = notImplemented0("dotr")
  
  def > (b : Float):Mat = notImplemented0(">")
  def < (b : Float):Mat = notImplemented0("<")
  def >= (b : Float):Mat = notImplemented0(">=")
  def <= (b : Float):Mat = notImplemented0("<=")
  def == (b : Float):Mat = notImplemented0("==")
  def === (b : Float):Mat = notImplemented0("===")
  def != (b : Float):Mat = notImplemented0("!=")
  
  def max(b : Float):Mat = notImplemented0("max")
  def min(b : Float):Mat = notImplemented0("min")
  
  def \ (b : Float):Mat = notImplemented0("\\")
  def on (b : Float):Mat = notImplemented0("on")
  
  
  def + (b : Int):Mat = notImplemented0("+")
  def - (b : Int):Mat = notImplemented0("-")
  def * (b : Int):Mat = notImplemented0("*")
  def xT (b : Int):Mat = notImplemented0("xT")
  def *^ (b : Int):Mat = notImplemented0("*^")
  def Tx (b : Int):Mat = notImplemented0("Tx")
  def ^* (b : Int):Mat = notImplemented0("*^")
  def /< (b : Int):Mat = notImplemented0("/<")
  def *@ (b : Int):Mat = notImplemented0("*@")
  def ∘  (b : Int):Mat = notImplemented0("∘")
  def /  (b : Int):Mat = notImplemented0("/")
  def \\ (b : Int):Mat = notImplemented0("\\\\")
  def ^ (b : Int):Mat = notImplemented0("^") 
  def ◁ (b : Int):Mat = notImplemented0("◁")
  def ▷ (b : Int):Mat = notImplemented0("▷")
  def dot (b : Int):Mat = notImplemented0("dot")
  def dotr (b : Int):Mat = notImplemented0("dotr")
  def ∙ (b : Int):Mat = notImplemented0("dot")
  def ∙→ (b : Int):Mat = notImplemented0("dotr")
  
  def > (b : Int):Mat = notImplemented0(">")
  def < (b : Int):Mat = notImplemented0("<")
  def >= (b : Int):Mat = notImplemented0(">=")
  def <= (b : Int):Mat = notImplemented0("<=")
  def == (b : Int):Mat = notImplemented0("==")
  def === (b : Int):Mat = notImplemented0("===")
  def != (b : Int):Mat = notImplemented0("!=")
  
  def max(b : Int):Mat = notImplemented0("max")
  def min(b : Int):Mat = notImplemented0("min")
  
  def \ (b : Int):Mat = notImplemented0("\\")
  def on (b : Int):Mat = notImplemented0("on")
  
  def + (b : Long):Mat = notImplemented0("+")
  def - (b : Long):Mat = notImplemented0("-")
  def * (b : Long):Mat = notImplemented0("*")
  def xT (b : Long):Mat = notImplemented0("xT")
  def *^ (b : Long):Mat = notImplemented0("*^")
  def Tx (b : Long):Mat = notImplemented0("Tx")
  def ^* (b : Long):Mat = notImplemented0("*^")
  def /< (b : Long):Mat = notImplemented0("/<")
  def *@ (b : Long):Mat = notImplemented0("*@")
  def ∘  (b : Long):Mat = notImplemented0("∘")
  def /  (b : Long):Mat = notImplemented0("/")
  def \\ (b : Long):Mat = notImplemented0("\\\\")
  def ^ (b : Long):Mat = notImplemented0("^") 
  def ◁ (b : Long):Mat = notImplemented0("◁")
  def ▷ (b : Long):Mat = notImplemented0("▷")
  def dot (b : Long):Mat = notImplemented0("dot")
  def dotr (b : Long):Mat = notImplemented0("dotr")
  def ∙ (b : Long):Mat = notImplemented0("dot")
  def ∙→ (b : Long):Mat = notImplemented0("dotr")
  
  def > (b : Long):Mat = notImplemented0(">")
  def < (b : Long):Mat = notImplemented0("<")
  def >= (b : Long):Mat = notImplemented0(">=")
  def <= (b : Long):Mat = notImplemented0("<=")
  def == (b : Long):Mat = notImplemented0("==")
  def === (b : Long):Mat = notImplemented0("===")
  def != (b : Long):Mat = notImplemented0("!=")
  
  def max(b : Long):Mat = notImplemented0("max")
  def min(b : Long):Mat = notImplemented0("min")
  
  def \ (b : Long):Mat = notImplemented0("\\")
  def on (b : Long):Mat = notImplemented0("on")
  
  def + (b : Double):Mat = notImplemented0("+")
  def - (b : Double):Mat = notImplemented0("-")
  def * (b : Double):Mat = notImplemented0("*")
  def xT (b : Double):Mat = notImplemented0("xT")
  def *^ (b : Double):Mat = notImplemented0("*^")
  def Tx (b : Double):Mat = notImplemented0("Tx")
  def ^* (b : Double):Mat = notImplemented0("*^")
  def /< (b : Double):Mat = notImplemented0("/<")
  def *@ (b : Double):Mat = notImplemented0("*@")
  def ∘  (b : Double):Mat = notImplemented0("∘")
  def /  (b : Double):Mat = notImplemented0("/")
  def \\ (b : Double):Mat = notImplemented0("\\\\")
  def ^ (b : Double):Mat = notImplemented0("^") 
  def ◁ (b : Double):Mat = notImplemented0("◁")
  def ▷ (b : Double):Mat = notImplemented0("▷")
  def dot (b : Double):Mat = notImplemented0("dot")
  def dotr (b : Double):Mat = notImplemented0("dotr")
  def ∙ (b : Double):Mat = notImplemented0("dot")
  def ∙→ (b : Double):Mat = notImplemented0("dotr")
  
  def > (b : Double):Mat = notImplemented0(">")
  def < (b : Double):Mat = notImplemented0("<")
  def >= (b : Double):Mat = notImplemented0(">=")
  def <= (b : Double):Mat = notImplemented0("<=")
  def == (b : Double):Mat = notImplemented0("==")
  def === (b : Double):Mat = notImplemented0("===")
  def != (b : Double):Mat = notImplemented0("!=")
 
  def max(b : Double):Mat = notImplemented0("max")
  def min(b : Double):Mat = notImplemented0("min")
  
  def \ (b : Double):Mat = notImplemented0("\\")
  def on (b : Double):Mat = notImplemented0("on")
  

}


object Mat {
  import Ordered._
  
  System.setProperty("java.util.logging.SimpleFormatter.format", "%1$tH:%1$tM:%1$tS %4$s: %5$s%n");
  
  val consoleLogger = Logger.getLogger("BIDMat console logging");  
  
  def getFileLogger(logfile:String):Logger = {
  	val fileHandler = new FileHandler(logfile);
  	val formatter = new SimpleFormatter();  
  	fileHandler.setFormatter(formatter); 
    val fileLogger = Logger.getLogger("BIDMat file logging");
    fileLogger.addHandler(fileHandler);
    fileLogger.setUseParentHandlers(false);
    fileLogger;
  }

  final val SyncMethod = jcuda.runtime.JCuda.cudaStreamPerThread;
//  final val SyncMethod = jcuda.runtime.JCuda.cudaStreamLegacy;
  
  var termWidth = 80;
  var youHaveBeenWarned = false;
  
  def terminalWidth:Int = {
  		if (!youHaveBeenWarned) {
  			try {
  				math.max(jline.TerminalFactory.create.getWidth, termWidth);
  			}
  			catch {
  			case _:Throwable => {
  				println("Couldnt get terminal width via JLine, using %d" format termWidth);
  				youHaveBeenWarned = true;
  				termWidth;
  			}
  			}
  		} else {
  			termWidth;
  		}
  }
  
  var useCache = false						 // Use matrix caching
  
  var useGPUcache = false
  
  var recycleGrow = 1.2            // For caching, amount to grow re-allocated matrices
  
  var hasCUDA = 0                  // Number of available CUDA GPUs
  
  var hasCUDNN = false;            // Is CUDNN available? 

  var cudartVersion = 0;
  
  var useBLAS = true;
    
  var useOpenCL: Boolean = true    // Use OpenCL if available

  var hasOpenCL: Boolean = false   // true if OpenCL runtime available

  var clHandle: CLHandle = null    // Holds OpenCL context and command queue
  
  var useMKL = true;               // Use MKL libs
  
  var useMKLRand = false;          // Use MKL random libs
    
  var useSTLRand = false;
    
  var debugMem = false             // Debug GPU mem calls

  var debugCPUmem = false          // Debug CPU mem calls
  
  var debugMemThreshold = 1000;
  
  var compressType = 1             // For HDF5 I/O, 0=none, 1=zlib, 2=szip
  
  var compressionLevel = 3         // for HDF5 zlib
  
  var chunkSize = 1024*1024        // for HDF5 compression
  
  var szipBlock = 32               // HDF5 szip block size
  
  var numThreads = Runtime.getRuntime().availableProcessors();
  
  var numOMPthreads = numThreads;
  
  var nflops = 0L
  
  var oneBased = 0                 // Whether matrix indices are 0: zero-based (like C) or 1: one-based (like Matlab)
  
  var ioneBased = 1                // Whether sparse matrix *internal* indices are zero 0: or one-based 1:
  
  var useGPUsort = true;
  
  var hostAllocSize = 4L*2047*1024*1024; // 2 billion floats
  
  var _inline = false;
  
  def inline = _inline;
  
  def setInline {
    _inline = true;
//    System.setProperty("java.awt.headless", "true"); // Ptplot wont work in headless mode
  }
  
  def unsetInline {
    _inline = false;
//    System.setProperty("java.awt.headless", "false");
  }
  
  final val MSEED:Int = 1452462553;

  final val myrand = new java.util.Random(MSEED);
  
  val opcodes = HashMap.empty[String, Int];
  
  val _opcode = 1;
  
  final val OS_WINDOWS = 0;
  final val OS_LINUX = 1;
  final val OS_OSX = 2;
  final val OS_ANDROID = 3;
  final val OS_UNKNOWN = 4;
  
  def getOS:Int = {
    val osname = System.getProperty("os.name");
    if (osname.startsWith("Windows")) OS_WINDOWS
    else if (osname.startsWith("Linux")) {
      if (System.getProperty("java.vendor").contains("ndroid")) OS_ANDROID
      else OS_LINUX      
    } else if (osname.startsWith("Mac")) OS_OSX
    else OS_UNKNOWN   
  }
  
  val ostype = getOS;
  
  var useStdio = (! (ostype == OS_WINDOWS));  // HDF5 directive
  
  private val _cache2 = HashMap.empty[Tuple2[Long,Int], Mat];             // Matrix caches
  
  private val _cache3 = HashMap.empty[Tuple3[Long,Long,Int], Mat];
  
  private val _cache4 = HashMap.empty[Tuple4[Long,Long,Long,Int], Mat];
  
  private val _cache5 = HashMap.empty[Tuple5[Long,Long,Long,Long,Int], Mat];
  
  def cache2(key:Tuple2[Long,Int]):Mat = {
    _cache2.synchronized {
    	if (_cache2.contains(key)) {
    		_cache2(key)
    	} else {
    		null
    	}
    }
  }
  
  def cache3(key:Tuple3[Long,Long,Int]):Mat = {
    _cache3.synchronized {
    	if (_cache3.contains(key)) {
    		_cache3(key)
    	} else {
    		null
    	}
    }
  }
  
    
  def cache4(key:Tuple4[Long,Long,Long,Int]):Mat = {
    _cache4.synchronized {
    	if (_cache4.contains(key)) {
    		_cache4(key)
    	} else {
    		null
    	}
    }
  }
  
  def cache5(key:Tuple5[Long,Long,Long,Long,Int]):Mat = {
    _cache5.synchronized {
    	if (_cache5.contains(key)) {
    		_cache5(key)
    	} else {
    		null
    	}
    }
  }
  
  def cache2put(key:Tuple2[Long,Int], m:Mat):Unit = {
    _cache2.synchronized {
    	_cache2(key) = m
    }
  }
  
  def cache3put(key:Tuple3[Long,Long,Int], m:Mat):Unit = {
  	_cache3.synchronized {
  		_cache3(key) = m
  	}
  }
  
  def cache4put(key:Tuple4[Long,Long,Long,Int], m:Mat):Unit = {
  	_cache4.synchronized {
  		_cache4(key) = m
  	}
  }
  
  def cache5put(key:Tuple5[Long,Long,Long,Long,Int], m:Mat):Unit = {
  	_cache5.synchronized {
  		_cache5(key) = m
  	}
  }
  
  def clearCaches = {
    _cache2.clear
    _cache3.clear
    _cache4.clear
    _cache5.clear
  }
  
  def trimCache2(ithread:Int) = {
    _cache2.synchronized {
      val keys = _cache2.keySet
      keys.foreach((key:Tuple2[Long,Int]) => {
      	val toremove:Boolean = _cache2.get(key).get match {
      	case aa:GMat => (aa.myGPU == ithread)
      	case aa:GDMat => (aa.myGPU == ithread)
      	case aa:GIMat => (aa.myGPU == ithread)
      	case aa:GLMat => (aa.myGPU == ithread)
      	case aa:GSMat => (aa.myGPU == ithread)
      	case aa:GSDMat => (aa.myGPU == ithread)
      	case _ => false
      }
      if (toremove) _cache2.remove(key)
      })  
    }
  }
  
   def trimCache3(ithread:Int) = {
    _cache3.synchronized {
      val keys = _cache3.keySet
      keys.foreach((key:Tuple3[Long,Long,Int]) => {
      	val toremove:Boolean = _cache3.get(key).get match {
      	case aa:GMat => (aa.myGPU == ithread)
      	case aa:GDMat => (aa.myGPU == ithread)
      	case aa:GIMat => (aa.myGPU == ithread)
      	case aa:GLMat => (aa.myGPU == ithread)
      	case aa:GSMat => (aa.myGPU == ithread)
      	case aa:GSDMat => (aa.myGPU == ithread)
      	case _ => false
      }
      if (toremove) _cache3.remove(key)
      })  
    }
  }
   
  def trimCache4(ithread:Int) = {
    _cache4.synchronized {
      val keys = _cache4.keySet
      keys.foreach((key:Tuple4[Long,Long,Long,Int]) => {
      	val toremove:Boolean = _cache4.get(key).get match {
      	case aa:GMat => (aa.myGPU == ithread)
      	case aa:GDMat => (aa.myGPU == ithread)
      	case aa:GIMat => (aa.myGPU == ithread)
      	case aa:GLMat => (aa.myGPU == ithread)
      	case aa:GSMat => (aa.myGPU == ithread)
      	case aa:GSDMat => (aa.myGPU == ithread)
      	case _ => false
      }
      if (toremove) _cache4.remove(key)
      })  
    }
  }
  
  def trimCache5(ithread:Int) = {
    _cache5.synchronized {
      val keys = _cache5.keySet
      keys.foreach((key:Tuple5[Long,Long,Long,Long,Int]) => {
      	val toremove:Boolean = _cache5.get(key).get match {
      	case aa:GMat => (aa.myGPU == ithread)
      	case aa:GDMat => (aa.myGPU == ithread)
      	case aa:GIMat => (aa.myGPU == ithread)
      	case aa:GLMat => (aa.myGPU == ithread)
      	case aa:GSMat => (aa.myGPU == ithread)
      	case aa:GSDMat => (aa.myGPU == ithread)
      	case _ => false
      }
      if (toremove) _cache5.remove(key)
      })  
    }
  }
  
  def trimCaches(ithread:Int) = {
    trimCache2(ithread);
    trimCache3(ithread);
    trimCache4(ithread);
    trimCache5(ithread);
  }
  
  def getJARdir:String = {
    val path = Mat.getClass.getProtectionDomain().getCodeSource().getLocation().getPath()
    val jstr = java.net.URLDecoder.decode(path, "UTF-8")
    path.replace("BIDMat.jar","")
  }
  
  def checkMKL:Unit = checkMKL(false)
  
  def checkMKL(verb:Boolean):Unit = {
    if (verb) useMKL = true;
    if (useMKL) {
    	try {
    		edu.berkeley.bid.LibUtils.loadLibrary("bidmatcpu", true);
    	} catch {
    	case th:Throwable => {
    		println("Cant find native CPU libraries")
    		if (verb) print(th.getMessage)
    		useMKL = false
    	}
    	}
    }
    useBLAS = useMKL;
    if (useMKL) useMKLRand = (UTILS.hasMKL() == 1);
    useSTLRand = useMKL & !useMKLRand;
    try {
      edu.berkeley.bid.LibUtils.loadLibrary("jhdf5")
//      System.loadLibrary("jhdf5")
    } catch {
    case _:Throwable => {
      println("Cant find native HDF5 library")
    }
    }
  }
  
  def checkCUDA:Unit = checkCUDA(false);

  def checkCUDA(verbose:Boolean):Unit = {
  		if (hasCUDA == 0) {
  			val os = System.getProperty("os.name");
  			try {
  				if (os.equals("Linux") || os.equals("Mac OS X")) {
  					System.loadLibrary("cudart");
					if (os.equals("Mac OS X")) {
					    System.loadLibrary("cublas");
					    System.loadLibrary("curand");
					    System.loadLibrary("cusparse");
					}
  				} else {
				    val libnames = List("cudart64_80", "cudart64_75", "cudart64_70", "cudart64_65",
				    	"cudart64_55", "cudart64_50_35", "cudart64_42_9").iterator;
  					var found = false;
  					while (!found && libnames.hasNext) {
  						found = true;
  						try{
  							System.loadLibrary(libnames.next);
  						} catch {
  						case _:Throwable => found = false;
  						}
  					}
  					if (!found) throw new RuntimeException("Couldnt find a cudart lib");
  				}
  			} catch {
  			case x:Throwable =>  {
  				println("Couldn't load CUDA runtime");
  				if (verbose) {
  					val msg = x.getMessage;
  					if (msg != null) println(msg);
  				}
  				hasCUDA = -1;    		
  			}
  			}
  			if (hasCUDA >= 0) {
  				try {
//  					jcuda.LibUtils.loadLibrary("JCudaRuntime");
					cudartVersion = jcuda.runtime.JCuda.CUDART_VERSION;  // Forces load of Jcuda class
  				} catch {
  				case y:Throwable =>  {
  					println("Couldn't load JCuda");
  					if (verbose) {
  						val msg = y.getMessage;
  						if (msg != null) println(msg);
  					}
  					hasCUDA = -1;   
  				}
  				}
  			}
  		}
  		if (hasCUDA >= 0) {
  			try {
  				var cudanum = new Array[Int](1);
  				val err = jcuda.runtime.JCuda.cudaGetDeviceCount(cudanum);
				if (err != 0) {
				       println("CUDA error %s" format jcuda.runtime.JCuda.cudaGetErrorString(err));
				}
  				hasCUDA = cudanum(0);
  				printf("%d CUDA device%s found", hasCUDA, if (hasCUDA == 1) "" else "s");
  				if (hasCUDA > 0) {
  					jcuda.runtime.JCuda.cudaRuntimeGetVersion(cudanum);
  					println(", CUDA version %d.%d" format (cudanum(0)/1000, (cudanum(0)%100) / 10));
  				} else {
  					println("");
  				}
  			} catch {
  			case e:NoClassDefFoundError => println("Couldn't load the JCUDA driver");
  			case e:Exception => println("Exception while initializing JCUDA driver");
  			case z:Throwable => println("Something went wrong while loading JCUDA driver" + z.getMessage);
  			}
  			if (hasCUDA > 0) {
  				try {
  					edu.berkeley.bid.LibUtils.loadLibrary("bidmatcuda");
  				} catch {
  				case z:Throwable => println("Something went wrong while loading BIDMat CUDA library" + z.getMessage);
  				}
  			}
  			SciFunctions.initCUDArngs;
  			GFunctions.initCublas;
            GFunctions.initCUDNN();
  		}
  }

  def checkOpenCL():Unit = checkOpenCL(false)

  def checkOpenCL(verbose: Boolean):Unit = {
    if (!useOpenCL) return
    // Attempt to load the runtime library
    if (!hasOpenCL) {
      // TODO: More robust method of loading JOCL and OpenCL runtime library
      try {
        // Enable OpenCL exceptions and also attempt to load JOCL (by running
        // the static block inside org.jocl.CL)
        // Note: JOCL also loads the OpenCL runtime
        org.jocl.CL.setExceptionsEnabled(true)
        hasOpenCL = true
      } catch {
        case err:Throwable => {
          println("Couldn't load OpenCL runtime")
          if (verbose) println(err.getMessage)
          hasOpenCL = false
        }
      }
    }

    // if we loaded the runtime successfully, query available devices
    if (hasOpenCL) {
      try {
        clHandle = CLHandle()
      } catch {
        case err:Throwable => {
          println("Couldn't acquire OpenCL context and command queue")
          if (verbose) println(err.getMessage)
          hasOpenCL = false
        }
      }
    }
  }

  def copyToIntArray[@specialized(Double, Float, Long, Byte, Short) T](data:Array[T], i0:Int, idata:Array[Int], d0:Int, n:Int)
  (implicit numeric : Numeric[T]) = {
    var i = 0 
    while (i < n) {
      idata(i+d0) = numeric.toInt(data(i+i0));
      i += 1
    }
  }
  
  def copyToDoubleArray[@specialized(Int, Float, Long, Byte, Short) T](data:Array[T], i0:Int, ddata:Array[Double], d0:Int, n:Int)
  (implicit numeric : Numeric[T]) = {
    var i = 0 
    while (i < n) {
      ddata(i+d0) = numeric.toDouble(data(i+i0));
      i += 1
    }
  }
  
  def copyToFloatArray[@specialized(Int, Double, Long, Byte, Short) T](data:Array[T], i0:Int, fdata:Array[Float], d0:Int, n:Int)
  (implicit numeric : Numeric[T]) = {
    var i = 0 
    while (i < n) {
      fdata(i+d0) = numeric.toFloat(data(i+i0));
      i += 1
    }
  }
  
  def copyToLongArray[@specialized(Int, Double, Float, Byte, Short) T](data:Array[T], i0:Int, fdata:Array[Long], d0:Int, n:Int)
  (implicit numeric : Numeric[T]) = {
    var i = 0 
    while (i < n) {
      fdata(i+d0) = numeric.toLong(data(i+i0));
      i += 1
    }
  }
  
 def copyToByteArray[@specialized(Double, Float, Long, Byte, Short) T](data:Array[T], i0:Int, idata:Array[Byte], d0:Int, n:Int)
  (implicit numeric : Numeric[T]) = {
    var i = 0 
    while (i < n) {
      idata(i+d0) = numeric.toInt(data(i+i0)).toByte;
      i += 1
    }
  }
  
  def copyListToFloatArray[T](a:List[T], b:Array[Float])(implicit numeric : Numeric[T]) = {
    var i = 0; 
    var todo = a.iterator
    val alen = a.length
    while (i < alen) {
      val h = todo.next
      b(i) = numeric.toFloat(h)
      i += 1
    }
  }
  
  def ibinsearch(v:Int, x:Array[Int], istartp:Int, iendp:Int):Int = {
    var istart = istartp
    var iend = iendp
    while (iend - istart > 1) {
      var mid:Int = (istart + iend)/2
      if (v < x(mid)) iend = mid else istart = mid
    }
    if (iend > istart && v == x(istart)) istart else -1
  }

  def binsearch[T : Ordering](v:T, x:Array[T], istartp:Int, iendp:Int):Int = {
    var istart = istartp
    var iend = iendp
    while (iend - istart > 1) {
      var mid:Int = (istart + iend)/2
      if (v < x(mid)) iend = mid else istart = mid
    }
    if (v == x(istart)) istart else -1
  }
  
  def lexsort[T :Ordering](a:List[Array[T]]):Array[Int] = {
    val n = a(0).length
    val ind = new Array[Int](n)
    var i = 0; while(i < n) {ind(i) = i; i += 1}
    def comp(i:Int, j:Int):Int = {
      val alen = a.length;
      val ip = ind(i)
      val jp = ind(j)
      var c0 = 0
      var k = 0;
      while (k < alen && c0 == 0) {
        c0 = a(k)(ip) compare a(k)(jp)
        k += 1
      }
      if (c0 != 0) {
        c0
      } else {
        ip compare jp
      }
    }
    def swap(i:Int, j:Int):Unit = {
      val tmp = ind(i)
      ind(i) = ind(j)
      ind(j) = tmp
    }
    BIDMat.Sorting.quickSort(comp, swap, 0, n)
    ind
  }

  def ilexsort(a:List[Array[Int]]):Array[Int] = {
    val n = a(0).length
    val ind = new Array[Int](n)
    var i = 0; while(i < n) {ind(i) = i; i += 1}
    def comp(i:Int, j:Int):Int = {
      var k = 0;
      val alen = a.length;
      var c0 = 0
      val ip = ind(i)
      val jp = ind(j)
      while (k < alen && c0 == 0) {
        c0 = a(k)(ip) compare a(k)(jp)
        k += 1
      }
      if (c0 != 0) {
        c0
      } else {
        ip compare jp
      }
    }
    def swap(i:Int, j:Int):Unit = {
      val tmp = ind(i)
      ind(i) = ind(j)
      ind(j) = tmp
    }
    BIDMat.Sorting.quickSort(comp, swap, 0, n)
    ind
  }
  
  def ilexsort2(a:Array[Int], b:Array[Int]):Array[Int] = {
    val n = a.length
    val ind = new Array[Int](n)
    var i = 0; while(i < n) {ind(i) = i; i += 1}
    def comp(i:Int, j:Int):Int = {
      val c0 = a(i) compare a(j)
      if (c0 != 0) {
        c0
      } else {
        val c1 = b(i) compare b(j)     
        if (c1 != 0) {
        	c1
        } else {
        	ind(i) compare ind(j)
        }         
      }
    }
    def swap(i:Int, j:Int):Unit = {
      val tmpa = a(i)
      a(i) = a(j)
      a(j) = tmpa
      val tmpb = b(i)
      b(i) = b(j)
      b(j) = tmpb
      val tmpi = ind(i)
      ind(i) = ind(j)
      ind(j) = tmpi
    }
    BIDMat.Sorting.quickSort(comp, swap, 0, n)
    ind
  }
  
  def ilexsort3[T](a:Array[Int], b:Array[Int], c:Array[T]):Unit = {
    val n = a.length
    def comp(i:Int, j:Int):Int = {
      val c0 = a(i) compare a(j)
      if (c0 != 0) {
        c0
      } else {
      	b(i) compare b(j)     
      }
    }
    def swap(i:Int, j:Int):Unit = {
      val tmpa = a(i)
      a(i) = a(j)
      a(j) = tmpa
      val tmpb = b(i)
      b(i) = b(j)
      b(j) = tmpb
      val tmpc = c(i)
      c(i) = c(j)
      c(j) = tmpc
    }
    BIDMat.Sorting.quickSort(comp, swap, 0, n)
  }
  
  def ilexsort(args:Array[Int]*):Array[Int] = {
    ilexsort(args.toList)
  }

  def lexsort[T : Ordering](args:Array[T]*):Array[Int] = {
    lexsort(args.toList)
  }

	def stackTraceString(e:Throwable):String = {
	  val sw = new java.io.StringWriter;
	  e.printStackTrace(new java.io.PrintWriter(sw));
	  sw.toString;
	}
  
}
