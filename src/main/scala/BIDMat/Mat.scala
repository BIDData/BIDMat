package BIDMat

class Mat(nr:Int, nc:Int) {
  val nrows = nr
  val ncols = nc

  def length = nr*nc
  
  val GUID = Mat.myrand.nextLong
  
  def notImplemented0(s:String):Mat = { 
    throw new RuntimeException("operator "+s+" not implemented for "+this.mytype)
  }
  
  def notImplemented1(s:String,that:Mat):Mat = { 
    throw new RuntimeException("operator "+s+" not implemented for "+this.mytype+" and "+that.mytype)
  }
  
  def t = notImplemented0("t")  
  def dv:Double = throw new RuntimeException("operator dv not implemented for "+this.mytype)
  
  def mytype = "Mat"
  def copyTo(a:Mat) = notImplemented0("copy");
  def copy = notImplemented0("copy");
  def set(v:Float) = notImplemented0("set")
  def zeros(nr:Int, nc:Int) = notImplemented0("zeros");
  def ones(nr:Int, nc:Int) = notImplemented0("ones");
  def clearUpper(i:Int) = notImplemented0("clearUpper");
  def clearLower(i:Int) = notImplemented0("clearLower"); 
  def clearUpper = notImplemented0("clearUpper");
  def clearLower = notImplemented0("clearLower");
    
  def nnz:Int = {notImplemented0("nnz"); 0}
  def clear = notImplemented0("clear");
  def zeros(nr:Int, nc:Int, nnz:Int):Mat = zeros(nr, nc)
  def recycle(nr:Int, nc:Int, nnz:Int):Mat = notImplemented0("recycle");
  def contents:Mat = notImplemented0("recycle");
  
  def apply(a:IMat):Mat = notImplemented0("linear array access");
  def apply(a:IMat, b:IMat):Mat = notImplemented0("block array access");
  def apply(a:IMat, b:Int):Mat = notImplemented0("block array access");	
  def apply(a:Int, b:IMat):Mat = notImplemented0("block array access");
  
  def update(a:IMat, b:Mat) = notImplemented0("linear update");
  def update(a:IMat, b:IMat, m:Mat) = notImplemented0("block update");
  def update(a:IMat, b:Int, m:Mat) = notImplemented0("block update");	
  def update(a:Int, b:IMat, m:Mat) = notImplemented0("block update");
  
  def + (b : Mat):Mat = notImplemented1("+", b)
  def - (b : Mat):Mat = notImplemented1("-", b)
  def * (b : Mat):Mat = notImplemented1("*", b)
  def *^ (b : Mat):Mat = notImplemented1("*^", b)
  def xT (b : Mat):Mat = notImplemented1("*", b)
  def Tx (b : Mat):Mat = notImplemented1("*", b)
  def / (b : Mat):Mat = notImplemented1("/", b)
  def *@ (b : Mat):Mat = notImplemented1("*@", b)
  def /@ (b : Mat):Mat = notImplemented1("/@", b)
  def \\ (b : Mat):Mat = notImplemented1("\\\\", b)
  def ^ (b : Mat):Mat = notImplemented1("^", b) 
  
  def > (b : Mat):Mat = notImplemented1(">", b)
  def < (b : Mat):Mat = notImplemented1("<", b)
  def >= (b : Mat):Mat = notImplemented1(">=", b)
  def <= (b : Mat):Mat = notImplemented1("<=", b)
  def == (b : Mat):Mat = notImplemented1("==", b)
  def === (b : Mat):Mat = notImplemented1("===", b)
  def != (b : Mat):Mat = notImplemented1("!=", b)
  
  def <-- (b : Mat):Mat = b.copyTo(this)
  
  def + (b : Int):Mat = notImplemented0("+")
  def - (b : Int):Mat = notImplemented0("-")
  def * (b : Int):Mat = notImplemented0("*")
  def / (b : Int):Mat = notImplemented0("/")
  def *@ (b : Int):Mat = notImplemented0("*@")
  def /@ (b : Int):Mat = notImplemented0("/@")
  def \\ (b : Int):Mat = notImplemented0("\\\\")
  def ^ (b : Int):Mat = notImplemented0("^") 
  
  def > (b : Int):Mat = notImplemented0(">")
  def < (b : Int):Mat = notImplemented0("<")
  def >= (b : Int):Mat = notImplemented0(">=")
  def <= (b : Int):Mat = notImplemented0("<=")
  def == (b : Int):Mat = notImplemented0("==")
  def === (b : Int):Mat = notImplemented0("===")
  def != (b : Int):Mat = notImplemented0("!=")
  
  def + (b : Float):Mat = notImplemented0("+")
  def - (b : Float):Mat = notImplemented0("-")
  def * (b : Float):Mat = notImplemented0("*")
  def / (b : Float):Mat = notImplemented0("/")
  def *@ (b : Float):Mat = notImplemented0("*@")
  def /@ (b : Float):Mat = notImplemented0("/@")
  def \\ (b : Float):Mat = notImplemented0("\\\\")
  def ^ (b : Float):Mat = notImplemented0("^") 
  
  def > (b : Float):Mat = notImplemented0(">")
  def < (b : Float):Mat = notImplemented0("<")
  def >= (b : Float):Mat = notImplemented0(">=")
  def <= (b : Float):Mat = notImplemented0("<=")
  def == (b : Float):Mat = notImplemented0("==")
  def === (b : Float):Mat = notImplemented0("===")
  def != (b : Float):Mat = notImplemented0("!=")
  
  def + (b : Double):Mat = notImplemented0("+")
  def - (b : Double):Mat = notImplemented0("-")
  def * (b : Double):Mat = notImplemented0("*")
  def / (b : Double):Mat = notImplemented0("/")
  def *@ (b : Double):Mat = notImplemented0("*@")
  def /@ (b : Double):Mat = notImplemented0("/@")
  def \\ (b : Double):Mat = notImplemented0("\\\\")
  def ^ (b : Double):Mat = notImplemented0("^") 
  
  def > (b : Double):Mat = notImplemented0(">")
  def < (b : Double):Mat = notImplemented0("<")
  def >= (b : Double):Mat = notImplemented0(">=")
  def <= (b : Double):Mat = notImplemented0("<=")
  def == (b : Double):Mat = notImplemented0("==")
  def === (b : Double):Mat = notImplemented0("===")
  def != (b : Double):Mat = notImplemented0("!=")  
  
  def \ (b : Mat):Mat = notImplemented1("\\", b)
  def on (b : Mat):Mat = notImplemented1("on", b)
  def ~ (b : Mat):Pair = b match {
    case bb:FMat => new FPair(this, bb)
    case bb:DMat => new DPair(this, bb)
    case bb:IMat => new IPair(this, bb)
    case bb:SMat => new SPair(this, bb)
//    case bb:SDMat => new SDPair(this, bb)
    case bb:CMat => new CPair(this, bb)
    case bb:GMat => new GPair(this, bb)
  }
  
  def ddot (b : Mat):Double = {notImplemented1("ddot", b); 0}
  
  def dot (b : Mat):Mat = {notImplemented1("dot", b); this}

}

abstract class Pair {
  
  def notImplemented0(s:String):Mat = { 
    throw new RuntimeException("operator "+s+" not implemented for "+this)
  }
  def notImplemented1(s:String,that:Mat):Mat = { 
    throw new RuntimeException("operator "+s+" not implemented for "+this+" and "+that)
  }
  
  def t = notImplemented0("t")
  
  def + (b : Mat):Mat = notImplemented1("+", b)
  def - (b : Mat):Mat = notImplemented1("-", b)
  def * (b : Mat):Mat = notImplemented1("*", b)
  def xT (b : Mat):Mat = notImplemented1("xT", b)
  def *^ (b : Mat):Mat = notImplemented1("*^", b)
  def Tx (b : Mat):Mat = notImplemented1("Tx", b)
  def / (b : Mat):Mat = notImplemented1("/", b)
  def *@ (b : Mat):Mat = notImplemented1("*@", b)
  def /@ (b : Mat):Mat = notImplemented1("/@", b)
  def \\ (b : Mat):Mat = notImplemented1("\\\\", b)
  def ^ (b : Mat):Mat = notImplemented1("^", b) 
  
  def > (b : Mat):Mat = notImplemented1(">", b)
  def < (b : Mat):Mat = notImplemented1("<", b)
  def >= (b : Mat):Mat = notImplemented1(">=", b)
  def <= (b : Mat):Mat = notImplemented1("<=", b)
  def == (b : Mat):Mat = notImplemented1("==", b)
  def === (b : Mat):Mat = notImplemented1("===", b)
  def != (b : Mat):Mat = notImplemented1("!=", b)
  
  def \ (b : Mat):Mat = notImplemented1("\\", b)
  def on (b : Mat):Mat = notImplemented1("on", b)
  
  def + (b : Int):Mat = notImplemented0("+")
  def - (b : Int):Mat = notImplemented0("-")
  def * (b : Int):Mat = notImplemented0("*")
  def / (b : Int):Mat = notImplemented0("/")
  def *@ (b : Int):Mat = notImplemented0("*@")
  def /@ (b : Int):Mat = notImplemented0("/@")
  def \\ (b : Int):Mat = notImplemented0("\\\\")
  def ^ (b : Int):Mat = notImplemented0("^") 
  
  def > (b : Int):Mat = notImplemented0(">")
  def < (b : Int):Mat = notImplemented0("<")
  def >= (b : Int):Mat = notImplemented0(">=")
  def <= (b : Int):Mat = notImplemented0("<=")
  def == (b : Int):Mat = notImplemented0("==")
  def === (b : Int):Mat = notImplemented0("===")
  def != (b : Int):Mat = notImplemented0("!=")
  
  def + (b : Float):Mat = notImplemented0("+")
  def - (b : Float):Mat = notImplemented0("-")
  def * (b : Float):Mat = notImplemented0("*")
  def / (b : Float):Mat = notImplemented0("/")
  def *@ (b : Float):Mat = notImplemented0("*@")
  def /@ (b : Float):Mat = notImplemented0("/@")
  def \\ (b : Float):Mat = notImplemented0("\\\\")
  def ^ (b : Float):Mat = notImplemented0("^") 
  
  def > (b : Float):Mat = notImplemented0(">")
  def < (b : Float):Mat = notImplemented0("<")
  def >= (b : Float):Mat = notImplemented0(">=")
  def <= (b : Float):Mat = notImplemented0("<=")
  def == (b : Float):Mat = notImplemented0("==")
  def === (b : Float):Mat = notImplemented0("===")
  def != (b : Float):Mat = notImplemented0("!=")
  
  def + (b : Double):Mat = notImplemented0("+")
  def - (b : Double):Mat = notImplemented0("-")
  def * (b : Double):Mat = notImplemented0("*")
  def / (b : Double):Mat = notImplemented0("/")
  def *@ (b : Double):Mat = notImplemented0("*@")
  def /@ (b : Double):Mat = notImplemented0("/@")
  def \\ (b : Double):Mat = notImplemented0("\\\\")
  def ^ (b : Double):Mat = notImplemented0("^") 
  
  def > (b : Double):Mat = notImplemented0(">")
  def < (b : Double):Mat = notImplemented0("<")
  def >= (b : Double):Mat = notImplemented0(">=")
  def <= (b : Double):Mat = notImplemented0("<=")
  def == (b : Double):Mat = notImplemented0("==")
  def === (b : Double):Mat = notImplemented0("===")
  def != (b : Double):Mat = notImplemented0("!=")  
}

object Mat {
  import Ordered._
  import scala.tools.jline.TerminalFactory
  
  var compressType = 1            // 0=none, 1=zlib, 2=szip
  
  var compressionLevel = 3        // for zlib
  
  var chunkSize = 1024*1024       // for either method
  
  var szipBlock = 32              // szip block size
  
  var numThreads = 1
  
  var noMKL:Boolean = false
  
  var nflops = 0L
  
  var oneBased = 0
  
  var ioneBased = 1
  
  final val MSEED:Int = 1452462553 

  final val myrand = new java.util.Random(MSEED)
  
  var hasCUDA = 0
  
  var useStdio = (! System.getProperty("os.name").startsWith("Windows"))
  
  def checkMKL:Unit = {
    if (!noMKL) {
      try {
      	System.loadLibrary("bidmatmkl")
      }	catch {
    			case _ => {
    				println("Cant find native CPU library")
    				noMKL = true
    			}
      }
    }
  }
  
  def checkCUDA:Unit = {
    if (hasCUDA == 0) {
    	try {
    		val os = System.getProperty("os.name")
    		if (os.equals("Linux")) {
    			System.loadLibrary("cudart")
    			System.loadLibrary("JCudaRuntime-linux-x86_64")
		} else if (os.equals("Mac OS X")) {
    			System.loadLibrary("cudart")
    			System.loadLibrary("JCudaRuntime-apple-x86_64")
    		} else {
    			try {
    				System.loadLibrary("cudart64_50_35")
    				System.loadLibrary("JCudaRuntime-windows-x86_64")
    			} catch {
    			case _ => try {
    				System.loadLibrary("cudart64_42_9")
    				System.loadLibrary("JCudaRuntime-windows-x86_64")
    			} 
    			}
    		}
    	} catch {
    	case _ =>  {
    		println("Cant find CUDA SDK or JCUDA")
    		hasCUDA = -1    		
    	}
    	}
    }
    if (hasCUDA >= 0) {
    	try {
    		var cudanum = new Array[Int](1)
    		jcuda.runtime.JCuda.cudaGetDeviceCount(cudanum)
    		hasCUDA = cudanum(0)
    		printf("%d CUDA device%s found", hasCUDA, if (hasCUDA == 1) "" else "s")
    		if (hasCUDA > 0) {
    			jcuda.runtime.JCuda.cudaRuntimeGetVersion(cudanum)
    			println(", CUDA version %d.%d" format (cudanum(0)/1000, (cudanum(0)%100) / 10))
    		} else {
    			println("")
    		}
    	} catch {
    	case e:NoClassDefFoundError => println("Couldn't load the JCUDA driver")
    	case e:Exception => println("Exception while initializing JCUDA driver")
    	case _ => println("Something went wrong while loading JCUDA driver")
    	}
    }
  }
  
  var terminal = TerminalFactory.create
  
  def terminalWidth = math.max(terminal.getWidth,80)

  def copyToIntArray[@specialized(Double, Float) T](data:Array[T], i0:Int, idata:Array[Int], d0:Int, n:Int)
  (implicit numeric : Numeric[T]) = {
    var i = 0 
    while (i < n) {
      idata(i+d0) = numeric.toInt(data(i+i0));
      i += 1
    }
  }
  
  def copyToDoubleArray[@specialized(Int, Float) T](data:Array[T], i0:Int, ddata:Array[Double], d0:Int, n:Int)
  (implicit numeric : Numeric[T]) = {
    var i = 0 
    while (i < n) {
      ddata(i+d0) = numeric.toDouble(data(i+i0));
      i += 1
    }
  }
  
  def copyToFloatArray[@specialized(Int, Double) T](data:Array[T], i0:Int, fdata:Array[Float], d0:Int, n:Int)
  (implicit numeric : Numeric[T]) = {
    var i = 0 
    while (i < n) {
      fdata(i+d0) = numeric.toFloat(data(i+i0));
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
    if (v == x(istart)) istart else -1
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
  
  def ilexsort(args:Array[Int]*):Array[Int] = {
    ilexsort(args.toList)
  }

  def lexsort[T : Ordering](args:Array[T]*):Array[Int] = {
    lexsort(args.toList)
  }
  
}
