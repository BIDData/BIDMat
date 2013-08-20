package BIDMat
import scala.collection.mutable.HashMap
import java.lang.ref._
import jcuda.NativePointerObject

class Mat(nr:Int, nc:Int) {
  val nrows = nr
  val ncols = nc

  def length = nr*nc
  
  private var _GUID = Mat.myrand.nextLong
  
  def setGUID(v:Long):Unit = {_GUID = v}
  
  def GUID:Long = _GUID
  
  def notImplemented0(s:String):Mat = { 
    throw new RuntimeException("operator "+s+" not implemented for "+this.mytype)
  }
  
  def notImplemented1(s:String,that:Mat):Mat = { 
    throw new RuntimeException("operator "+s+" not implemented for "+this.mytype+" and "+that.mytype)
  }
  
  def t = notImplemented0("t")  
  def dv:Double = throw new RuntimeException("operator dv not implemented for "+this.mytype)
  
  def mytype = "Mat"
  def copyTo(a:Mat) = notImplemented0("copyTo");
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
  def contents:Mat = notImplemented0("contents");
  def colslice(a:Int, b:Int, out:Mat):Mat = notImplemented0("colslice");
  def colslice(a:Int, b:Int, out:Mat, c:Int):Mat = notImplemented0("colslice");
  def rowslice(a:Int, b:Int, out:Mat):Mat = notImplemented0("rowslice");
  def rowslice(a:Int, b:Int, out:Mat, c:Int):Mat = notImplemented0("rowslice");
  
  def apply(a:IMat):Mat = notImplemented0("linear array access");
  def apply(a:IMat, b:IMat):Mat = notImplemented0("block array access");
  def apply(a:IMat, b:Int):Mat = notImplemented0("block array access");	
  def apply(a:Int, b:IMat):Mat = notImplemented0("block array access");
  
  def update(a:IMat, b:Mat) = notImplemented0("linear update");
  def update(a:IMat, b:IMat, m:Mat) = notImplemented0("block update");
  def update(a:IMat, b:Int, m:Mat) = notImplemented0("block update");	
  def update(a:Int, b:IMat, m:Mat) = notImplemented0("block update");
  
  def unary_-():Mat = notImplemented1("-", this)
  def + (b : Mat):Mat = notImplemented1("+", b)
  def - (b : Mat):Mat = notImplemented1("-", b)
  def * (b : Mat):Mat = notImplemented1("*", b)
  def *^ (b : Mat):Mat = notImplemented1("*^", b)
  def xT (b : Mat):Mat = notImplemented1("*", b)
  def Tx (b : Mat):Mat = notImplemented1("*", b)
  def ^* (b : Mat):Mat = notImplemented1("*^", b)
  def /< (b : Mat):Mat = notImplemented1("/<", b)
  def *@ (b : Mat):Mat = notImplemented1("*@", b)
  def ∘  (b : Mat):Mat = notImplemented1("∘", b)
  def / (b : Mat):Mat = notImplemented1("/", b)
  def \\ (b : Mat):Mat = notImplemented1("\\\\", b)
  def ^ (b : Mat):Mat = notImplemented1("^", b) 
  def ◁ (b : Mat):Mat = notImplemented1("◁", b)
  def ▷ (b : Mat):Mat = notImplemented1("▷", b)
  def dot (b : Mat):Mat = notImplemented1("dot", b)
  def dotr (b : Mat):Mat = notImplemented1("dotr", b) 
  def ∙ (b : Mat):Mat = notImplemented1("dot", b)
  def ∙∙ (b : Mat):Mat = notImplemented1("dotr", b)
  
  def > (b : Mat):Mat = notImplemented1(">", b)
  def < (b : Mat):Mat = notImplemented1("<", b)
  def >= (b : Mat):Mat = notImplemented1(">=", b)
  def <= (b : Mat):Mat = notImplemented1("<=", b)
  def == (b : Mat):Mat = notImplemented1("==", b)
  def === (b : Mat):Mat = notImplemented1("===", b)
  def != (b : Mat):Mat = notImplemented1("!=", b)
  
  def <-- (b : Mat):Mat = b.copyTo(this)
  
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
  

  
  def ^* (b : DSPair):Mat = notImplemented0("^*")
  def Tx (b : DSPair):Mat = notImplemented0("Tx")
  def @@ (b : Mat):DSPair = (this, b) match {
    case (aa:FMat, bb:SMat) => new FDSPair(aa, bb) 
    case (aa:GMat, bb:GSMat) => new GDSPair(aa, bb)
  }
}

abstract class DSPair {}

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
  def ∙∙ (b : Mat):Mat = notImplemented1("dotr", b)
  
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
  def ∙∙ (b : Float):Mat = notImplemented0("dotr")
  
  def > (b : Float):Mat = notImplemented0(">")
  def < (b : Float):Mat = notImplemented0("<")
  def >= (b : Float):Mat = notImplemented0(">=")
  def <= (b : Float):Mat = notImplemented0("<=")
  def == (b : Float):Mat = notImplemented0("==")
  def === (b : Float):Mat = notImplemented0("===")
  def != (b : Float):Mat = notImplemented0("!=")
  
  def \ (b : Float):Mat = notImplemented0("\\")
  def on (b : Float):Mat = notImplemented0("on")
  

}


object Mat {
  import Ordered._
  import scala.tools.jline.TerminalFactory
  
  var useCache = false						// Use expression caching
  
  var compressType = 1            // 0=none, 1=zlib, 2=szip
  
  var compressionLevel = 3        // for zlib
  
  var chunkSize = 1024*1024       // for either method
  
  var szipBlock = 32              // szip block size
  
  var numThreads = Runtime.getRuntime().availableProcessors();
  
  var numOMPthreads = numThreads;
  
  var noMKL:Boolean = false
  
  var nflops = 0L
  
  var oneBased = 0
  
  var ioneBased = 1
  
  var recycleGrow = 1.2
  
  var debugMem = false
  
  var useGPUsort = true
  
  final val MSEED:Int = 1452462553 

  final val myrand = new java.util.Random(MSEED)
  
  private val _cache2 = HashMap.empty[Tuple2[Long,Int], Mat]
  
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
  }
  
  def trimCache2(ithread:Int) = {
    _cache2.synchronized {
      val keys = _cache2.keySet
      keys.foreach((key:Tuple2[Long,Int]) => {
      	val toremove:Boolean = _cache2.get(key) match {
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
      	val toremove:Boolean = _cache3.get(key) match {
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
      	val toremove:Boolean = _cache4.get(key) match {
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
  
  val opcodes = HashMap.empty[String, Int]
  
  val _opcode = 1
  
  var hasCUDA = 0
  
  var useStdio = (! System.getProperty("os.name").startsWith("Windows"))
  
  def getJARdir:String = {
    val path = Mat.getClass.getProtectionDomain().getCodeSource().getLocation().getPath()
    val jstr = java.net.URLDecoder.decode(path, "UTF-8")
    path.replace("BIDMat.jar","")
  }
  
  def checkMKL:Unit = {
    if (!noMKL) {
      try {      	
      	System.loadLibrary("bidmatmkl")
      }	catch {
      case _ =>
      try {
      	val os = System.getProperty("os.name")
      	if (os.equals("Linux")) {
      		System.load(getJARdir+"lib/linux64/libbidmatmkl.so")
      		System.load(getJARdir+"lib/linux64/libjhdf5.so")
      	} else if (os.equals("Mac OS X")) {
      		System.load(getJARdir+"lib/osx64/libbidmatmkl.so")
      		System.load(getJARdir+"lib/osx64/libjhdf5.so")      		      		
      	} else {
      		System.load(getJARdir.replace("/c:","C:")+"lib/win64/bidmatmkl.dll")
      		System.load(getJARdir.replace("/c:","C:")+"lib/win64/jhdf5.dll")
      	}
      } catch {
      case _ => {
      	println("Cant find native CPU libraries")
      	noMKL = true
      }
      }
      }
    }
  }
  
  def checkCUDA:Unit = {
    if (hasCUDA == 0) {
    	val os = System.getProperty("os.name")
    	try {
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
    	case _ => try {
    		val JCuda_libs = List("JCudaRuntime", "JCublas2", "JCublas", "JCudaDriver", "JCufft", "JCurand", "JCusparse2", "JCusparse")
    		if (os.equals("Linux")) {
    		  JCuda_libs.foreach(i => System.load(getJARdir+"lib/linux64/JCUDA5.0/lib"+i+"-linux-x86_64.so"))
    		  System.load(getJARdir+"lib/linux64/JCUDA5.0/libbidmatcuda.so")
    		} else if (os.contains("Mac")) {
    		  JCuda_libs.foreach(i => System.load(getJARdir+"lib/osx64/JCUDA5.0/lib"+i+"-apple-x86_64.dylib"))
    		  System.load(getJARdir+"lib/osx64/JCUDA5.0/libbidmatcuda.dylib")
    		} else {
    		  JCuda_libs.foreach(i => System.load(getJARdir+"lib/win64/JCUDA5.0/"+i+"-windows-x86_64.dll"))
    		  System.load(getJARdir+"lib/win64/JCUDA5.0/bidmatcuda.dll")
    		}
    	} catch {
    	case _ =>  {
    		println("Cant find CUDA SDK or JCUDA")
    		hasCUDA = -1    		
    	}
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
