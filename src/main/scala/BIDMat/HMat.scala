package BIDMat

import java.io._
import java.util.zip._
import scala.util.matching.Regex
import Regex._
import scala.collection.mutable._
import scala.actors._
import scala.actors.Actor._
import MatFunctions._
import MatHDF5._
import edu.berkeley.bid.UTILS._

case class HMat(nr:Int, nc:Int, fileList:List[String], varname:String, blkinds:Array[Int], catdim:Int) extends Mat(nr, nc) {

  var fnameCache:String = null

  var fmatCache:Mat = null
  
  override def mytype = "HMat"
  
// Implement slicing from a hard disk matrix
  override def apply(a:IMat, b:IMat):Mat = { 
    var ilast:Int = 0
    def findindx(ind:Int):Int = {
    	while (ilast >= 0 && ind < blkinds(ilast)) ilast -= 1
    	while (ilast < blkinds.length && ind >= blkinds(ilast)) ilast += 1
    	if (ilast >= 0 && blkinds(ilast) <= ind && blkinds(ilast+1) > ind) {
    		ilast
    	} else {
    		-1
      }
    }
    
    val locs = IMat(1,b.length)
    var i = 0
    var iblk = 0
    var out:Mat = null
    while (i <= b.length) {
    	if (i < b.length) locs(i) = findindx(b(i))
    	if (i == b.length || locs(i) != locs(iblk)) {
    		if (fnameCache == null || fileList(locs(iblk)) != fnameCache) {
    			fmatCache = MatHDF5.hload(fileList(locs(iblk)), varname).asInstanceOf[Mat] 
    			fnameCache = fileList(locs(iblk))
    		}
      	val newmat = fmatCache(a, b(MatFunctions.irow(iblk->i)))
      	if (out.asInstanceOf[AnyRef] != null) {
    	  	out = out \ newmat
        } else {
          out = newmat
        }
        iblk = i + 1
      }
      i += 1
    }
    out
  }
}

object HMat {
  
  def readSomeInts(din:InputStream, a:Array[Int], buf:Array[Byte], n:Int) {
    var nread = 0
    while (nread < 4*n) {
      val readnow = din.read(buf, 0, math.min(buf.length, 4*n-nread))
      memcpybi(readnow, buf, 0, a, nread)
      nread += readnow
    }
  }
  
  def readSomeFloats(din:InputStream, a:Array[Float], buf:Array[Byte], n:Int) {
    var nread = 0
    while (nread < 4*n) {
      val readnow = din.read(buf, 0, math.min(buf.length, 4*n-nread))
      memcpybf(readnow, buf, 0, a, nread)
      nread += readnow
    }
  }
  
  def readSomeDoubles(din:InputStream, a:Array[Double], buf:Array[Byte], n:Int) {
    var nread = 0
    while (nread < 8*n) {
      val readnow = din.read(buf, 0, math.min(buf.length, 8*n-nread))
      memcpybd(readnow, buf, 0, a, nread)
      nread += readnow
    }
  }
  
  def writeSomeInts(dout:OutputStream, a:Array[Int], buf:Array[Byte], n:Int) {
    var nwritten = 0
    while (nwritten < 4*n) {
      val todo = math.min(4*n-nwritten, buf.length)
    	memcpyib(todo, a, nwritten, buf, 0)
      dout.write(buf, 0, todo)
      nwritten += todo
    }
  }
  
  def writeSomeFloats(dout:OutputStream, a:Array[Float], buf:Array[Byte], n:Int) {
    var nwritten = 0
    while (nwritten < 4*n) {
      val todo = math.min(4*n-nwritten, buf.length)
    	memcpyfb(todo, a, nwritten, buf, 0)
      dout.write(buf, 0, todo)
      nwritten += todo
    }
  }
  
  def writeSomeDoubles(dout:OutputStream, a:Array[Double], buf:Array[Byte], n:Int) {
    var nwritten = 0
    while (nwritten < 8*n) {
      val todo = math.min(8*n-nwritten, buf.length)
    	memcpydb(todo, a, nwritten, buf, 0)
      dout.write(buf, 0, todo)
      nwritten += todo
    }
  }
  
  def getInputStream(fname:String, compressed:Boolean):InputStream = {
    val fin = new FileInputStream(fname)
    if (compressed) {
      new GZIPInputStream(fin, 1024*1024)
    } else {
      new BufferedInputStream(fin, 1024*1024)
    }
  }
  
  def getOutputStream(fname:String, compressed:Boolean):OutputStream = {
    import edu.berkeley.bid.UTILS._
  	_getOutputStream(fname, compressed, Mat.compressionLevel)
  }
  
  def loadFMat(fname:String, compressed:Boolean=true):FMat = {
    val gin = getInputStream(fname, compressed)
    val buff = new Array[Byte](1024*1024)
    val hints = new Array[Int](4)
    readSomeInts(gin, hints, buff, 4)
    val ftype = hints(0)
    val nrows = hints(1)
    val ncols = hints(2)
    val out = FMat(nrows, ncols)
    readSomeFloats(gin, out.data, buff, ncols*nrows)
    gin.close
    out
  }
   
  def loadIMat(fname:String, compressed:Boolean=true):IMat = {
    val gin = getInputStream(fname, compressed)
    val buff = new Array[Byte](1024*1024)
    val hints = new Array[Int](4)
    readSomeInts(gin, hints, buff, 4)
    val ftype = hints(0)
    val nrows = hints(1)
    val ncols = hints(2)
    val out = IMat(nrows, ncols)
    readSomeInts(gin, out.data, buff, ncols*nrows)
    gin.close
    out
  }
   
  def loadDMat(fname:String, compressed:Boolean=true):DMat = {
    val gin = getInputStream(fname, compressed)
    val buff = new Array[Byte](1024*1024)
    val hints = new Array[Int](4)
    readSomeInts(gin, hints, buff, 4)
    val ftype = hints(0)
    val nrows = hints(1)
    val ncols = hints(2)
    val out = DMat(nrows, ncols)
    readSomeDoubles(gin, out.data, buff, ncols*nrows)
    gin.close
    out
  }
  
  def saveFMat(fname:String, m:FMat, compressed:Boolean=true):Unit = {
    val gout = getOutputStream(fname, compressed)
    val hints = new Array[Int](4)
    val tbuf = new Array[Byte](16)
    hints(0) = 130 // 1=dense, 3=float
    hints(1) = m.nrows
    hints(2) = m.ncols
    hints(3) = 0
    writeSomeInts(gout, hints, tbuf, 4)
    val buff = new Array[Byte](math.min(1024*1024, 4*m.ncols*m.nrows))
    writeSomeFloats(gout, m.data, buff, m.nrows*m.ncols)
    gout.close
  }
  
  def saveIMat(fname:String, m:IMat, compressed:Boolean=true):Unit = {
  	val gout = getOutputStream(fname, compressed)
    val hints = new Array[Int](4)
    val tbuf = new Array[Byte](16)
    hints(0) = 110 // 1=dense, 1=int
    hints(1) = m.nrows
    hints(2) = m.ncols
    hints(3) = 0
    writeSomeInts(gout, hints, tbuf, 4)
    val buff = new Array[Byte](math.min(1024*1024, 4*m.ncols*m.nrows))
    writeSomeInts(gout, m.data, buff, m.nrows*m.ncols)
    gout.close
  }
  
  def saveDMat(fname:String, m:DMat, compressed:Boolean=true):Unit = {
    val gout = getOutputStream(fname, compressed)
    val hints = new Array[Int](4)
    val tbuf = new Array[Byte](16)
    hints(0) = 140 // 1=dense, 4=double
    hints(1) = m.nrows
    hints(2) = m.ncols
    hints(3) = 0
    writeSomeInts(gout, hints, tbuf, 4)
    val buff = new Array[Byte](math.min(1024*1024, 4*m.ncols*m.nrows))
    writeSomeDoubles(gout, m.data, buff, m.nrows*m.ncols)
    gout.close
  }
  
  def loadSMat(fname:String, compressed:Boolean=true):SMat = {
    val gin = getInputStream(fname, compressed)
    val buff = new Array[Byte](1024*1024)
    val hints = new Array[Int](4)
    readSomeInts(gin, hints, buff, 4)
    val ftype = hints(0)
    val nrows = hints(1)
    val ncols = hints(2)
    val nnz = hints(3)
    val out = SMat(nrows, ncols, nnz)
    readSomeInts(gin, out.jc, buff, ncols+1)
    readSomeInts(gin, out.ir, buff, nnz)
    readSomeFloats(gin, out.data, buff, nnz)
    MatHDF5.addOne(out.jc)
    MatHDF5.addOne(out.ir)
    gin.close
    out
  }
  
  def saveSMat(fname:String, m:SMat, compressed:Boolean=true):Unit = {
    val gout = getOutputStream(fname, compressed)
    val hints = new Array[Int](4)
    val tbuf = new Array[Byte](16)
    hints(0) = 231 // 2=sparse, 3=float, 1=int
    hints(1) = m.nrows
    hints(2) = m.ncols
    hints(3) = m.nnz
    writeSomeInts(gout, hints, tbuf, 4)
    val buff = new Array[Byte](math.min(1024*1024, 4*math.max(m.ncols+1, m.nnz)))
    try {
    	MatHDF5.subOne(m.jc)
    	MatHDF5.subOne(m.ir)
    	writeSomeInts(gout, m.jc, buff, m.ncols+1)
    	writeSomeInts(gout, m.ir, buff, m.nnz)
    	writeSomeFloats(gout, m.data, buff, m.nnz)
    } catch {
      case e:Exception => {
      	MatHDF5.addOne(m.jc)
      	MatHDF5.addOne(m.ir)
      	throw new RuntimeException("Exception in saveSMat "+e)
      }
      case _ => {
      	MatHDF5.addOne(m.jc)
      	MatHDF5.addOne(m.ir)
      	throw new RuntimeException("Problem in saveSMat")
      }
    }
    MatHDF5.addOne(m.jc)
    MatHDF5.addOne(m.ir)
    gout.close
  } 
  
  def testLoad(fname:String, varname:String, n:Int) = {
    val a = new Array[SMat](n)
    var ndone = izeros(n,1)
    for (i <- 0 until n) {
      actor {
        a(i) = loadSMat(("/disk%02d/" format i)+fname)
        ndone(i) = 1
      }
    }
    while (SciFunctions.sum(ndone).v < n) {Thread.sleep(10)}
    a
  }
  
  def apply(dirname:String, filepat:String, varname:String, catd:Int) { 
    var files:ListBuffer[String] = new ListBuffer[String]
    val dir:File = new File(dirname)
    val slen = dir.getName.length + 1
    
    def searchDir(dir:File) {
      for (f <- dir.listFiles) {
        if (f.isDirectory) {
	      searchDir(f)
        } else { 
          if (f.getName.substring(slen).matches(filepat)) { 
            files.append(f.getName)
          }
        }
      }
    }

    searchDir(dir)
    val blkinds = new Array[Int](files.length+1)
    var i = 0
    var nrows = -1
    var ncols = -1
    files.foreach((fn:String) => { 
      val (nr, nc) = MatHDF5.readMatDims(fn, varname)
      if (catd == 2) { 
        if (nrows >= 0) { 
          if (nr != nrows) { 
            throw new RuntimeException("incorrect number of rows in file "+fn)
          }
        } else { 
          nrows = nr.asInstanceOf[Int]
        }
        blkinds(i+1) = blkinds(i) + nc.asInstanceOf[Int]
        i += 1
      } else if (catd == 1) { 
        if (ncols >= 0) { 
          if (nc != 1) { 
            throw new RuntimeException("incorrect number of cols in file "+fn)
          }
        } else { 
          ncols = 1
        }
        blkinds(i+1) = blkinds(i) + nr.asInstanceOf[Int]
        i += 1
      }
    })
    if (catd == 2) { 
      HMat(nrows, blkinds(files.length), files.toList, varname, blkinds, 2)
    } else {
      if (catd == 1) { 
      	HMat(blkinds(files.length), ncols, files.toList, varname, blkinds, 1)
      } else {
      	throw new RuntimeException("cat dimension must be 1 or 2")
      }
    }
  }
}






