package BIDMat
import scala.collection.mutable.HashMap
import java.lang.ref._
import jcuda.NativePointerObject
import edu.berkeley.bid.UTILS

@SerialVersionUID(100L)
class Mat(nr:Int, nc:Int) extends ND(Array(nr, nc)) with Serializable {
  override val nrows = nr
  override val ncols = nc
  
  def llength = 1L*nr*nc
  
  private var _GUID = Mat.myrand.nextLong
  
  def notImplemented0(s:String):Mat = { 
    throw new RuntimeException("operator "+s+" not implemented for "+this.mytype)
  }
  
   def notImplemented0ND(s:String):ND = { 
    throw new RuntimeException("operator "+s+" not implemented for "+this.mytype)
  }
  
  def notImplemented1(s:String,that:ND):ND = { 
    throw new RuntimeException("operator "+s+" not implemented for "+this.mytype+" and "+that.mytype)
  }
  
  def notImplemented1(s:String,that:Mat):Mat = { 
    throw new RuntimeException("operator "+s+" not implemented for "+this.mytype+" and "+that.mytype)
  }
  
  def notImplemented2(s:String,that:Float):Mat = { 
    throw new RuntimeException("operator "+s+" not implemented for "+this.mytype+" and Float")
  }
  
  def notImplemented2(s:String,that:Double):Mat = { 
    throw new RuntimeException("operator "+s+" not implemented for "+this.mytype+" and Double")
  }
  
  def notImplementedf(s:String):Float = { 
    throw new RuntimeException("operator "+s+" not implemented for "+this.mytype)
  }
  
  def notImplementedd(s:String):Double = { 
    throw new RuntimeException("operator "+s+" not implemented for "+this.mytype)
  }
  
  def notImplementedi(s:String):Int = { 
    throw new RuntimeException("operator "+s+" not implemented for "+this.mytype)
  }
  
  def notImplementedl(s:String):Long = { 
    throw new RuntimeException("operator "+s+" not implemented for "+this.mytype)
  }
  
  def applyf(indx:Int):Float  = throw new RuntimeException("1D access not supported for "+this.mytype);
  
  def apply(i1:Mat, i2:Mat, i3:Mat):ND = throw new RuntimeException("3D access not supported for "+this.mytype);
  def apply(i1:Mat, i2:Mat, i3:Mat, i4:Mat):ND = throw new RuntimeException("4D access not supported for "+this.mytype);
  
  def update(i1:Mat, i2:Mat, vv:ND):ND = throw new RuntimeException("2D updates not supported for "+this.mytype);
  def update(i1:Mat, i2:Mat, i3:Mat, vv:ND):ND = throw new RuntimeException("3D updates not supported for "+this.mytype);
  def update(i1:Mat, i2:Mat, i3:Mat, i4:Mat, vv:ND):ND = throw new RuntimeException("4D access not supported for "+this.mytype);
  
  def t = notImplemented0("t")  
  def dv:Double = throw new RuntimeException("operator dv not implemented for "+this.mytype)
  
  override def mytype = "Mat"
  def copyTo(a:Mat) = notImplemented0("copyTo");
  def copy = notImplemented0("copy");
  def newcopy = notImplemented0("newcopy");
  def set(v:Float) = notImplemented0("set");
  def set(v:Double) = notImplemented0("set");
  def zeros(nr:Int, nc:Int) = notImplemented0("zeros");
  def zeros(dims0:IMat) = notImplemented0ND("zeros");
  def zeros = notImplemented0("zeros");
  def ones(nr:Int, nc:Int) = notImplemented0("ones");
  def ones(dims0:IMat) = ones(dims0(0), dims0(1));
  def izeros(nr:Int, nc:Int) = notImplemented0("izeros");
  def iones(nr:Int, nc:Int) = notImplemented0("iones");
  def clear = notImplemented0("clear");
  def clearUpper(i:Int) = notImplemented0("clearUpper");
  def clearLower(i:Int) = notImplemented0("clearLower"); 
  def clearUpper = notImplemented0("clearUpper");
  def clearLower = notImplemented0("clearLower");
  def free = notImplemented0("free");
  def view(nr:Int, nc:Int):Mat = notImplemented0("view");
  def view(nr:Int, nc:Int, setGUID:Boolean):Mat = notImplemented0("view");
    
  def nnz:Int = {notImplemented0("nnz"); 0}
  def zeros(nr:Int, nc:Int, nnz:Int):Mat = zeros(nr, nc);
  def recycle(nr:Int, nc:Int, nnz:Int):Mat = notImplemented0("recycle");
  def contents:Mat = notImplemented0("contents");
  def colslice(a:Int, b:Int, out:Mat):Mat = notImplemented0("colslice");
  def colslice(a:Int, b:Int, out:Mat, c:Int):Mat = notImplemented0("colslice");
  def colslice(a:Int, b:Int, out:Mat, c:Int, pb:Boolean):Mat = colslice(a, b, out, c);
  def rowslice(a:Int, b:Int, out:Mat):Mat = notImplemented0("rowslice");
  def rowslice(a:Int, b:Int, out:Mat, c:Int):Mat = notImplemented0("rowslice");
  def colslice(a:Int, b:Int):Mat = notImplemented0("colslice");
  def rowslice(a:Int, b:Int):Mat = notImplemented0("rowslice");
  
  def colslice(a:Int, b:Int, out:ND):ND = notImplemented0("colslice");
  def colslice(a:Int, b:Int, out:ND, c:Int):ND = notImplemented0("colslice");
  
  def apply(a:IMat):Mat = notImplemented0("linear array access");
  def apply(a:IMat, b:IMat):Mat = notImplemented0("block array access");
  def apply(a:IMat, b:Int):Mat = notImplemented0("block array access");	
  def apply(a:Int, b:IMat):Mat = notImplemented0("block array access");
  
  def apply(a:GIMat):Mat = notImplemented0("linear array access");
  def apply(a:GIMat, b:GIMat):Mat = notImplemented0("block array access");
  def apply(a:GIMat, b:Int):Mat = notImplemented0("block array access");	
  def apply(a:Int, b:GIMat):Mat = notImplemented0("block array access");
  
  def apply(a:IMat, b:GIMat):Mat = notImplemented0("block array access");
  def apply(a:GIMat, b:IMat):Mat = notImplemented0("block array access");
  
  def apply(a:Mat):Mat = notImplemented0("linear array access");
  def apply(a:Mat, b:Int):Mat = notImplemented0("block array access");
  def apply(a:Int, b:Mat):Mat = notImplemented0("block array access");
  def apply(a:Mat, b:Mat):Mat = notImplemented0("block array access");
  
  def update(a:IMat, b:Mat) = notImplemented0("linear update");
  def update(a:IMat, b:IMat, m:Mat) = notImplemented0("block update");
  def update(a:IMat, b:Int, m:Mat) = notImplemented0("block update");	
  def update(a:Int, b:IMat, m:Mat) = notImplemented0("block update");
  
  def update(a:IMat, b:Int) = notImplemented0("linear update");
  def update(a:IMat, b:IMat, c:Int) = notImplemented0("block update");
  def update(a:IMat, b:Int, c:Int) = notImplemented0("block update");   
  def update(a:Int, b:IMat, c:Int) = notImplemented0("block update"); 
    
  def update(a:IMat, b:Float) = notImplemented0("linear update");
  def update(a:IMat, b:IMat, c:Float) = notImplemented0("block update");
  def update(a:IMat, b:Int, c:Float) = notImplemented0("block update");   
  def update(a:Int, b:IMat, c:Float) = notImplemented0("block update");
  
  def update(a:IMat, b:Double) = notImplemented0("linear update");
  def update(a:IMat, b:IMat, c:Double) = notImplemented0("block update");
  def update(a:IMat, b:Int, c:Double) = notImplemented0("block update");   
  def update(a:Int, b:IMat, c:Double) = notImplemented0("block update");
  
  def update(a:IMat, b:Long) = notImplemented0("linear update");
  def update(a:IMat, b:IMat, c:Long) = notImplemented0("block update");
  def update(a:IMat, b:Int, c:Long) = notImplemented0("block update");   
  def update(a:Int, b:IMat, c:Long) = notImplemented0("block update");

  def update(a:GIMat, b:Mat) = notImplemented0("linear update");
  def update(a:GIMat, b:GIMat, m:Mat) = notImplemented0("block update");
  def update(a:GIMat, b:Int, m:Mat) = notImplemented0("block update");	
  def update(a:Int, b:GIMat, m:Mat) = notImplemented0("block update");
  def update(a:GIMat, b:IMat, m:Mat) = notImplemented0("block update");
  def update(a:IMat, b:GIMat, m:Mat) = notImplemented0("block update");
  
  def update(a:GIMat, b:Int) = notImplemented0("linear update");
  def update(a:GIMat, b:GIMat, c:Int) = notImplemented0("block update");
  def update(a:GIMat, b:Int, c:Int) = notImplemented0("block update");  
  def update(a:Int, b:GIMat, c:Int) = notImplemented0("block update");
  def update(a:GIMat, b:IMat, c:Int) = notImplemented0("block update");
  def update(a:IMat, b:GIMat, c:Int) = notImplemented0("block update");
  
  def update(a:GIMat, b:Float) = notImplemented0("linear update");
  def update(a:GIMat, b:GIMat, c:Float) = notImplemented0("block update");
  def update(a:GIMat, b:Int, c:Float) = notImplemented0("block update");  
  def update(a:Int, b:GIMat, c:Float) = notImplemented0("block update");
  def update(a:GIMat, b:IMat, c:Float) = notImplemented0("block update");
  def update(a:IMat, b:GIMat, c:Float) = notImplemented0("block update");
  
  def update(a:GIMat, b:Double) = notImplemented0("linear update");
  def update(a:GIMat, b:GIMat, c:Double) = notImplemented0("block update");
  def update(a:GIMat, b:Int, c:Double) = notImplemented0("block update");  
  def update(a:Int, b:GIMat, c:Double) = notImplemented0("block update");
  def update(a:GIMat, b:IMat, c:Double) = notImplemented0("block update");
  def update(a:IMat, b:GIMat, c:Double) = notImplemented0("block update");
  
  def update(a:Mat, v:Mat):Mat = notImplemented0("linear update");
  def update(a:Mat, b:Mat, v:Mat):Mat = notImplemented0("block update");
  def update(a:Mat, b:Int, v:Mat):Mat = notImplemented0("block update");
  def update(a:Int, b:Mat, v:Mat):Mat = notImplemented0("block update");
  
  def update(a:Mat, v:Int):Mat = notImplemented0("linear update");
  def update(a:Mat, b:Mat, v:Int):Mat = notImplemented0("block update");
  def update(a:Mat, b:Int, v:Int):Mat = notImplemented0("block update");
  def update(a:Int, b:Mat, v:Int):Mat = notImplemented0("block update");
  
  def update(a:Mat, v:Float):Mat = notImplemented0("linear update");
  def update(a:Mat, b:Mat, v:Float):Mat = notImplemented0("block update");
  def update(a:Mat, b:Int, v:Float):Mat = notImplemented0("block update");
  def update(a:Int, b:Mat, v:Float):Mat = notImplemented0("block update");
  
  def update(a:Mat, v:Double):Mat = notImplemented0("linear update");
  def update(a:Mat, b:Mat, v:Double):Mat = notImplemented0("block update");
  def update(a:Mat, b:Int, v:Double):Mat = notImplemented0("block update");
  def update(a:Int, b:Mat, v:Double):Mat = notImplemented0("block update");
  
  def update(a:Mat, v:Long):Mat = notImplemented0("linear update");
  def update(a:Mat, b:Mat, v:Long):Mat = notImplemented0("block update");
  def update(a:Mat, b:Int, v:Long):Mat = notImplemented0("block update");
  def update(a:Int, b:Mat, v:Long):Mat = notImplemented0("block update");
  
  def update(a:Int, v:Float):Mat = notImplemented0("linear update");
  def update(a:Int, v:Double):Mat = notImplemented0("linear update");
  def update(a:Int, v:Int):Mat = notImplemented0("linear update");
  def update(a:Int, v:Long):Mat = notImplemented0("linear update");
  
  def update(a:Int, b:Int, v:Float):Mat = notImplemented0("update");
  def update(a:Int, b:Int, v:Double):Mat = notImplemented0("update");
  def update(a:Int, b:Int, v:Int):Mat = notImplemented0("update");
  def update(a:Int, b:Int, v:Long):Mat = notImplemented0("update");  
  def tileMult(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:Mat, broff:Int, bcoff:Int, c:Mat, croff:Int, ccoff:Int):Mat =
    notImplemented0("tileMult");
  def tileMultT(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:Mat, broff:Int, bcoff:Int, c:Mat, croff:Int, ccoff:Int):Mat =
    notImplemented0("tileMultT");
  
  def blockGemm(transa:Int, transb:Int, nr:Int, nc:Int, reps:Int, aoff:Int, lda:Int, astep:Int, 
      b:Mat, boff:Int, ldb:Int, bstep:Int, c:Mat, coff:Int, ldc:Int, cstep:Int):Mat = notImplemented0("blockGemm");
  
  def copyTo(a:ND):ND = notImplemented1("<--", a);
  
  val asMat:Mat = this;

  def madd(a:Mat, b:Mat, at:Boolean, bt:Boolean):Mat = notImplemented1("update", a);
  def madd(a:Mat, b:Mat):Mat = notImplemented1("update", a);
  
  def unary_-():Mat = notImplemented1("-", this)
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
  
  def <-- (b : Mat):Mat = b.copyTo(this)
  def \ (b : Mat):Mat = notImplemented1("\\", b)
  def on (b : Mat):Mat = notImplemented1("on", b)
  
  def *  (b : Float):Mat = notImplemented2("*", b)
  def +  (b : Float):Mat = notImplemented2("+", b)
  def -  (b : Float):Mat = notImplemented2("-", b)
  def *@ (b : Float):Mat = notImplemented2("*@", b)
  def ∘  (b : Float):Mat = notImplemented2("∘", b)
  def / (b : Float):Mat = notImplemented2("/", b)
  def ^  (b : Float):Mat = notImplemented2("^", b)
  
  def >  (b : Float):Mat = notImplemented2(">", b)
  def <  (b : Float):Mat = notImplemented2("<", b)
  def >= (b : Float):Mat = notImplemented2(">=", b)
  def <= (b : Float):Mat = notImplemented2("<=", b)
  def == (b : Float):Mat = notImplemented2("==", b)
  def === (b : Float):Mat = notImplemented2("===", b)
  def != (b : Float):Mat = notImplemented2("!=", b)
  
  def *  (b : Int):Mat = notImplemented2("*", b)
  def +  (b : Int):Mat = notImplemented2("+", b)
  def -  (b : Int):Mat = notImplemented2("-", b)
  def *@ (b : Int):Mat = notImplemented2("*@", b)
  def ∘  (b : Int):Mat = notImplemented2("∘", b)
  def / (b : Int):Mat = notImplemented2("/", b)
  def ^  (b : Int):Mat = notImplemented2("^", b)
  
  def >  (b : Int):Mat = notImplemented2(">", b)
  def <  (b : Int):Mat = notImplemented2("<", b)
  def >= (b : Int):Mat = notImplemented2(">=", b)
  def <= (b : Int):Mat = notImplemented2("<=", b)
  def == (b : Int):Mat = notImplemented2("==", b)
  def === (b : Int):Mat = notImplemented2("===", b)
  def != (b : Int):Mat = notImplemented2("!=", b)
  
  def + (b : Long):Mat = notImplemented2("+", b)
  def - (b : Long):Mat = notImplemented2("-", b)
  def * (b : Long):Mat = notImplemented2("*", b)
  def *@ (b : Long):Mat = notImplemented2("*@", b)
  def ∘  (b : Long):Mat = notImplemented2("∘", b)
  def /  (b : Long):Mat = notImplemented2("/", b)
  def ^ (b : Long):Mat   = notImplemented2("^", b)
  
  def > (b : Long):Mat = notImplemented2(">", b)
  def < (b : Long):Mat = notImplemented2("<", b)
  def >= (b : Long):Mat = notImplemented2(">=", b)
  def <= (b : Long):Mat = notImplemented2("<=", b)
  def == (b : Long):Mat = notImplemented2("==", b)
  def === (b : Long):Mat = notImplemented2("===", b)
  def != (b : Long):Mat = notImplemented2("!=", b)
  
  def *  (b : Double):Mat = notImplemented2("*", b)
  def +  (b : Double):Mat = notImplemented2("+", b)
  def -  (b : Double):Mat = notImplemented2("-", b)
  def *@ (b : Double):Mat = notImplemented2("*@", b)
  def ∘  (b : Double):Mat = notImplemented2("∘", b)
  def / (b : Double):Mat = notImplemented2("/", b)
  def ^  (b : Double):Mat = notImplemented2("^", b)
  
  def >  (b : Double):Mat = notImplemented2(">", b)
  def <  (b : Double):Mat = notImplemented2("<", b)
  def >= (b : Double):Mat = notImplemented2(">=", b)
  def <= (b : Double):Mat = notImplemented2("<=", b)
  def == (b : Double):Mat = notImplemented2("==", b)
  def === (b : Double):Mat = notImplemented2("===", b)
  def != (b : Double):Mat = notImplemented2("!=", b)
  
  def +  (b : ND):ND = notImplemented1("+", b)
  def -  (b : ND):ND = notImplemented1("-", b)
  def *  (b : ND):ND = notImplemented1("*", b)
  def *^ (b : ND):ND = notImplemented1("*^", b)
  def xT (b : ND):ND = notImplemented1("*", b)
  def Tx (b : ND):ND = notImplemented1("*", b)
  def ^* (b : ND):ND = notImplemented1("^*", b)
  def ** (b : ND):ND = notImplemented1("**", b)
  def ⊗  (b : ND):ND = notImplemented1("⊗", b)       // unicode 8855, 0x2297
  def /< (b : ND):ND = notImplemented1("/<", b)
  def ∘  (b : ND):ND = notImplemented1("∘", b)        // unicode 8728, 0x2218 
  def *@ (b : ND):ND = notImplemented1("*@", b)
  def /  (b : ND):ND = notImplemented1("/", b)
  def \\ (b : ND):ND = notImplemented1("\\\\", b)
  def ^  (b : ND):ND = notImplemented1("^", b) 
  def ◁  (b : ND):ND = notImplemented1("◁", b)        // unicode 9665, 0x25C1 
  def ▷  (b : ND):ND = notImplemented1("▷", b)        // unicode 9666, 0x25C2
  def dot (b : ND):ND = notImplemented1("dot", b)
  def dotr (b : ND):ND = notImplemented1("dotr", b) 
  def ∙ (b : ND):ND = notImplemented1("dot", b)       // unicode 8729, 0x2219 
  def ∙→ (b : ND):ND = notImplemented1("dotr", b)     // unicode (8729, 8594) (0x2219, 0x2192)
    
  def >  (b : ND):ND = notImplemented1(">", b)
  def <  (b : ND):ND = notImplemented1("<", b)
  def >= (b : ND):ND = notImplemented1(">=", b)
  def <= (b : ND):ND = notImplemented1("<=", b)
  def == (b : ND):ND = notImplemented1("==", b)
  def === (b : ND):ND = notImplemented1("===", b)
  def != (b : ND):ND = notImplemented1("!=", b)
  
//  def <-- (b : ND):ND = b.copyTo(this)
  def \ (b : ND):ND = notImplemented1("\\", b)
  def on (b : ND):ND = notImplemented1("on", b)
  
  
  def ddot (b : Mat):Double = {notImplemented1("ddot", b); 0}
  def ∙∙ (b : Mat):Double = {notImplemented1("ddot", b); 0}
  
  
  def ^* (b : DSPair):Mat = notImplemented0("^*")
  def Tx (b : DSPair):Mat = notImplemented0("Tx")
  def @@ (b : Mat):DSPair = (this, b) match {
    case (aa:FMat, bb:SMat) => new FDSPair(aa, bb) 
    case (aa:GMat, bb:GSMat) => new GDSPair(aa, bb)
  }
  
  def ~ (b : Mat):Pair = b match {
    case bb:FMat => new FPair(this, bb)
    case bb:DMat => new DPair(this, bb)
    case bb:IMat => new IPair(this, bb)
    case bb:SMat => new SPair(this, bb)
    case bb:SDMat => new SDPair(this, bb)
    case bb:CMat => new CPair(this, bb)
    case bb:GMat => new GPair(this, bb)
    case bb:GIMat => new GIPair(this, bb)
    case bb:GDMat => new GDPair(this, bb)
    case bb:GLMat => new GLPair(this, bb)
  }
    
  def ~ (b:ND):Pair = this ~ b.asInstanceOf[Mat]
  
}

abstract class DSPair extends Serializable {}

abstract class Pair extends NDPair {
  
  def notImplemented0(s:String):Mat = { 
    throw new RuntimeException("operator "+s+" not implemented for "+this)
  }
  
  def notImplemented1(s:String,that:Mat):Mat = { 
    throw new RuntimeException("operator "+s+" not implemented for "+this+" and "+that.mytype)
  }
  
  def notImplemented1(s:String,that:ND):Mat = { 
    throw new RuntimeException("operator "+s+" not implemented for "+this+" and "+that.mytype)
  }
  
  def t = notImplemented0("t")
  
  def + (b : Mat):Mat = notImplemented1("+", b)
  def - (b : Mat):Mat = notImplemented1("-", b)
  def * (b : Mat):Mat = notImplemented1("*", b)
  def xT (b : Mat):Mat = notImplemented1("xT", b)
  def *^ (b : Mat):Mat = notImplemented1("*^", b)
  def Tx (b : Mat):Mat = notImplemented1("Tx", b)
  def ^* (b : Mat):Mat = notImplemented1("*^", b)
  def /< (b : Mat):Mat = notImplemented1("/<", b)
  def *@ (b : Mat):Mat = notImplemented1("*@", b)
  def ∘  (b : Mat):Mat = notImplemented1("∘", b)
  def /  (b : Mat):Mat = notImplemented1("/", b)
  def \\ (b : Mat):Mat = notImplemented1("\\\\", b)
  def ^ (b : Mat):Mat = notImplemented1("^", b) 
  def ◁ (b : Mat):Mat = notImplemented1("◁", b)
  def ▷ (b : Mat):Mat = notImplemented1("▷", b)
  def dot (b : Mat):Mat = notImplemented1("dot", b)
  def dotr (b : Mat):Mat = notImplemented1("dotr", b)
  def ∙ (b : Mat):Mat = notImplemented1("dot", b)
  def ∙→ (b : Mat):Mat = notImplemented1("dotr", b)
  def ** (b : Mat):Mat = notImplemented1("**", b)
  def ⊗  (b : Mat):Mat = notImplemented1("⊗", b) 
  
  def > (b : Mat):Mat = notImplemented1(">", b)
  def < (b : Mat):Mat = notImplemented1("<", b)
  def >= (b : Mat):Mat = notImplemented1(">=", b)
  def <= (b : Mat):Mat = notImplemented1("<=", b)
  def == (b : Mat):Mat = notImplemented1("==", b)
  def === (b : Mat):Mat = notImplemented1("===", b)
  def != (b : Mat):Mat = notImplemented1("!=", b)
  
  def \ (b : Mat):Mat = notImplemented1("\\", b)
  def on (b : Mat):Mat = notImplemented1("on", b)
  
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
  
  def \ (b : Double):Mat = notImplemented0("\\")
  def on (b : Double):Mat = notImplemented0("on")
  
  
  def + (b : ND):Mat = notImplemented1("+", b)
  def - (b : ND):Mat = notImplemented1("-", b)
  def * (b : ND):Mat = notImplemented1("*", b)
  def xT (b : ND):Mat = notImplemented1("xT", b)
  def *^ (b : ND):Mat = notImplemented1("*^", b)
  def Tx (b : ND):Mat = notImplemented1("Tx", b)
  def ^* (b : ND):Mat = notImplemented1("*^", b)
  def /< (b : ND):Mat = notImplemented1("/<", b)
  def *@ (b : ND):Mat = notImplemented1("*@", b)
  def ∘  (b : ND):Mat = notImplemented1("∘", b)
  def /  (b : ND):Mat = notImplemented1("/", b)
  def \\ (b : ND):Mat = notImplemented1("\\\\", b)
  def ^ (b : ND):Mat = notImplemented1("^", b) 
  def ◁ (b : ND):Mat = notImplemented1("◁", b)
  def ▷ (b : ND):Mat = notImplemented1("▷", b)
  def dot (b : ND):Mat = notImplemented1("dot", b)
  def dotr (b : ND):Mat = notImplemented1("dotr", b)
  def ∙ (b : ND):Mat = notImplemented1("dot", b)
  def ∙→ (b : ND):Mat = notImplemented1("dotr", b)
  def ** (b : ND):Mat = notImplemented1("**", b)
  def ⊗  (b : ND):Mat = notImplemented1("⊗", b) 
  
  def > (b : ND):Mat = notImplemented1(">", b)
  def < (b : ND):Mat = notImplemented1("<", b)
  def >= (b : ND):Mat = notImplemented1(">=", b)
  def <= (b : ND):Mat = notImplemented1("<=", b)
  def == (b : ND):Mat = notImplemented1("==", b)
  def === (b : ND):Mat = notImplemented1("===", b)
  def != (b : ND):Mat = notImplemented1("!=", b)
  
  def \ (b : ND):Mat = notImplemented1("\\", b)
  def on (b : ND):Mat = notImplemented1("on", b)
  

}


object Mat {
  import Ordered._
  
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
  
  var recycleGrow = 1.2            // For caching, amount to grow re-allocated matrices
  
  var hasCUDA = 0                  // Number of available CUDA GPUs
  
  var useBLAS = true;
    
  var useOpenCL: Boolean = true    // Use OpenCL if available

  var hasOpenCL: Boolean = false   // true if OpenCL runtime available

  var clShutdownHookSet: Boolean = false // true if we've set a shutdown hook to clean up OpenCL resources

  var clResourcesFreed: Boolean = false // true if we've freed the OpenCL context & command queue

  var numOpenCLGPUs = 0            // number of available OpenCL GPUs

  var clContext: org.jocl.cl_context = null // OpenCL context

  var clQueue: org.jocl.cl_command_queue = null // Current OpenCL device command queue
  
  var useMKL = true;               // Use MKL libs
  
  var useMKLRand = false;          // Use MKL random libs
    
  var useSTLRand = false;
    
  var debugMem = false             // Debug GPU mem calls
  
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
  
  val ostype = getOS  
  
  var useStdio = (! (ostype == OS_WINDOWS))  // HDF5 directive
  
  private val _cache2 = HashMap.empty[Tuple2[Long,Int], Mat]              // Matrix caches
  
  private val _cache3 = HashMap.empty[Tuple3[Long,Long,Int], Mat]
  
  private val _cache4 = HashMap.empty[Tuple4[Long,Long,Long,Int], Mat]
  
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
  
  def clearCaches = {
    _cache2.clear
    _cache3.clear
    _cache4.clear
    ND.clearCaches
  }
  
  def trimCache2(ithread:Int) = {
    _cache2.synchronized {
      val keys = _cache2.keySet
      keys.foreach((key:Tuple2[Long,Int]) => {
      	val toremove:Boolean = _cache2.get(key).get match {
      	case aa:GMat => (aa.myGPU == ithread) 
      	case aa:GSMat => (aa.myGPU == ithread)
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
      	case aa:GSMat => (aa.myGPU == ithread)
      	case _ => false
      }
      if (toremove) _cache3.remove(key)
      })  
    }
  }
   
  def trimCache4(ithread:Int) = {
    _cache3.synchronized {
      val keys = _cache4.keySet
      keys.foreach((key:Tuple4[Long,Long,Long,Int]) => {
      	val toremove:Boolean = _cache4.get(key).get match {
      	case aa:GMat => (aa.myGPU == ithread) 
      	case aa:GSMat => (aa.myGPU == ithread)
      	case _ => false
      }
      if (toremove) _cache4.remove(key)
      })  
    }
  }
  
  def trimCaches(ithread:Int) = {
    trimCache2(ithread)
    trimCache3(ithread)
    trimCache4(ithread)
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
    	  if (ostype == OS_WINDOWS) edu.berkeley.bid.LibUtils.loadLib("libiomp5md")
    		edu.berkeley.bid.LibUtils.loadLibrary("bidmatcpu")
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
//      jcuda.LibUtils.loadLibrary("jhdf5")
      System.loadLibrary("jhdf5")
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
  				} else {
  					val libnames = List("cudart64_70", "cudart64_65", "cudart64_55", "cudart64_50_35", "cudart64_42_9").iterator;
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
  				println("Couldnt load CUDA runtime");
  				if (verbose) {
  					val msg = x.getMessage;
  					if (msg != null) println(msg);
  				}
  				hasCUDA = -1;    		
  			}
  			}
  			if (hasCUDA >= 0) {
  				try {
  					jcuda.LibUtils.loadLibrary("JCudaRuntime");
  				} catch {
  				case y:Throwable =>  {
  					println("Couldnt load JCuda");
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
  				jcuda.runtime.JCuda.cudaGetDeviceCount(cudanum);
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
  					jcuda.LibUtils.loadLibrary("bidmatcuda");
  				} catch {
  				case z:Throwable => println("Something went wrong while loading BIDMat CUDA library" + z.getMessage);
  				}
  			}
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
          println("Failed to load JOCL native shared library")
          if (verbose) {
            val msg = err.getMessage
            if (msg != null) println(msg)
          }
          hasOpenCL = false
        }
      }
    }
    // if we loaded the runtime successfully, query available devices
    if (hasOpenCL) {
      val platforms = getCLPlatforms()
      val platform = platforms(0)
      //val gpus = getCLDevices(platform, org.jocl.CL.CL_DEVICE_TYPE_GPU)
      val gpus = getCLDevices(platform, org.jocl.CL.CL_DEVICE_TYPE_CPU)
      numOpenCLGPUs = gpus.length
      clContext = createCLContext(platform, gpus)
      clQueue = createCLQueue(clContext, gpus(0))
      clResourcesFreed = false

      // Make sure to clean up before shutdown
      if (!clShutdownHookSet) {
        clShutdownHookSet = true
        Runtime.getRuntime().addShutdownHook(new Thread() {
          override def run() = {
            freeOpenCL()
          }
        })
      }
    }
  }

  def freeOpenCL():Unit = {
    synchronized {
      if (!clResourcesFreed) {
        CLKernelCache.free()
        org.jocl.CL.clReleaseCommandQueue(clQueue)
        org.jocl.CL.clReleaseContext(clContext)
        clResourcesFreed = true
        hasOpenCL = false
      }
    }
  }

  def getCLPlatforms():Array[org.jocl.cl_platform_id] = {
    val num_platforms_ptr = Array(0)
    org.jocl.CL.clGetPlatformIDs(0, null, num_platforms_ptr)
    val num_platforms = num_platforms_ptr(0)
    val platforms = Array.ofDim[org.jocl.cl_platform_id](num_platforms)
    org.jocl.CL.clGetPlatformIDs(platforms.length, platforms, null)
    platforms
  }

  def getCLDevices(platform: org.jocl.cl_platform_id, device_type: Long):Array[org.jocl.cl_device_id] = {
    val num_devices_ptr = Array(0)
    org.jocl.CL.clGetDeviceIDs(platform, device_type.toInt, 0, null, num_devices_ptr)
    val num_devices = num_devices_ptr(0)
    val devices = Array.ofDim[org.jocl.cl_device_id](num_devices)
    org.jocl.CL.clGetDeviceIDs(platform, device_type.toInt, num_devices, devices, null)
    devices
  }

  def createCLContext(platform: org.jocl.cl_platform_id, devices: Array[org.jocl.cl_device_id]):org.jocl.cl_context = {
    val properties = new org.jocl.cl_context_properties()
    properties.addProperty(org.jocl.CL.CL_CONTEXT_PLATFORM, platform)
    org.jocl.CL.clCreateContext(properties, devices.length, devices, null, null, null)
  }

  def createCLQueue(context: org.jocl.cl_context, device: org.jocl.cl_device_id):org.jocl.cl_command_queue = {
    val properties = new org.jocl.cl_queue_properties()
    properties.addProperty(org.jocl.CL.CL_QUEUE_PROPERTIES, org.jocl.CL.CL_QUEUE_PROFILING_ENABLE)
    org.jocl.CL.clCreateCommandQueueWithProperties(context, device, properties, null)
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
  
}
