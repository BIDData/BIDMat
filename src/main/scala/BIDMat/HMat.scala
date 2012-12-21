package BIDMat

import java.io._
import java.util.zip._
import java.nio._
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
  
  val byteOrder = ByteOrder.LITTLE_ENDIAN
  
  def readSomeInts(din:InputStream, a:Array[Int], buf:ByteBuffer, n:Int) {
    var nread = 0
    val ibuff = buf.asIntBuffer
    val bbuff = buf.array
    while (nread < n) {
      val todo = if (n - nread > ibuff.capacity) ibuff.capacity else (n - nread)
      var readnow = din.read(bbuff, 0, todo*4)
      while (readnow % 4 != 0) {
        readnow += din.read(bbuff, readnow, 4 - (readnow % 4))
        println("blah")
      }
//      println("%d %d %d %d %d" format (nread, readnow, todo, ibuff.capacity, bbuff.length))
      ibuff.get(a, nread, readnow/4)
      ibuff.position(0)
      nread += readnow/4
    }
  }
  
  def readSomeFloats(din:InputStream, a:Array[Float], buf:ByteBuffer, n:Int) {
    var nread = 0
    val fbuff = buf.asFloatBuffer
    val bbuff = buf.array
    while (nread < n) {
      val todo = if (n - nread > fbuff.capacity) fbuff.capacity else (n - nread)
      var readnow = din.read(bbuff, 0, todo*4)
      while (readnow % 4 != 0) {
        readnow += din.read(bbuff, readnow, 4 - (readnow % 4))
      }
      fbuff.get(a, nread, readnow/4)
      fbuff.position(0)
      nread += readnow/4  
    }
  }
  
  def readSomeDoubles(din:InputStream, a:Array[Double], buf:ByteBuffer, n:Int) {
    var nread = 0
    val dbuff = buf.asDoubleBuffer
    val bbuff = buf.array
    while (nread < n) {
      val todo = if (n - nread > dbuff.capacity) dbuff.capacity else (n - nread)
      var readnow = din.read(bbuff, 0, todo*8)
      while (readnow % 8 != 0) {
        readnow += din.read(bbuff, readnow, 8 - (readnow % 8))
      }
      dbuff.get(a, nread, readnow/8)
      dbuff.position(0)
      nread += readnow/8 
    }
  }
  
  def writeSomeInts(dout:OutputStream, a:Array[Int], buf:ByteBuffer, n:Int) {
    var nwritten = 0
    val ibuff = buf.asIntBuffer
    val bbuff = buf.array
    while (nwritten < n) {
    	val todo = if (n - nwritten > ibuff.capacity) ibuff.capacity else (n - nwritten)
      ibuff.put(a, nwritten, todo)
      dout.write(bbuff, 0, todo*4)
      nwritten += todo
    }
  }
  
  def writeSomeFloats(dout:OutputStream, a:Array[Float], buf:ByteBuffer, n:Int) {
    var nwritten = 0
    val fbuff = buf.asFloatBuffer
    val bbuff = buf.array
    while (nwritten < n) {
    	val todo = if (n - nwritten > fbuff.capacity) fbuff.capacity else (n - nwritten)
      fbuff.put(a, nwritten, todo)
      dout.write(bbuff, 0, todo*4)
      nwritten += todo
    }
  }
  
  def writeSomeDoubles(dout:OutputStream, a:Array[Double], buf:ByteBuffer, n:Int) {
    var nwritten = 0
    val dbuff = buf.asDoubleBuffer
    val bbuff = buf.array
    while (nwritten < n) {
    	val todo = if (n - nwritten > dbuff.capacity) dbuff.capacity else (n - nwritten)
    	dbuff.put(a, nwritten, todo)
      dout.write(bbuff, 0, todo*8)
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
    val buff = ByteBuffer.allocate(1024*1024).order(byteOrder)
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
    val buff = ByteBuffer.allocate(1024*1024).order(byteOrder)
    val hints = new Array[Int](4)
    readSomeInts(gin, hints, buff, 4)
    val ftype = hints(0)
    val nrows = hints(1)
    val ncols = hints(2)
    println("%d %d %d\n" format (ftype, nrows, ncols))
    val out = IMat(nrows, ncols)
    readSomeInts(gin, out.data, buff, ncols*nrows)
    gin.close
    out
  }
   
  def loadDMat(fname:String, compressed:Boolean=true):DMat = {
    val gin = getInputStream(fname, compressed)
    val bytebuff = ByteBuffer.allocate(1024*1024).order(byteOrder)
    val hints = new Array[Int](4)
    readSomeInts(gin, hints, bytebuff, 4)
    val ftype = hints(0)
    val nrows = hints(1)
    val ncols = hints(2)
    val out = DMat(nrows, ncols)
    readSomeDoubles(gin, out.data, bytebuff, ncols*nrows)
    gin.close
    out
  }
  
  def saveFMat(fname:String, m:FMat, compressed:Boolean=true):Unit = {
    val gout = getOutputStream(fname, compressed)
    val hints = new Array[Int](4)
    val tbuf = ByteBuffer.allocate(16).order(byteOrder)
    hints(0) = 130 // 1=dense, 3=float
    hints(1) = m.nrows
    hints(2) = m.ncols
    hints(3) = 0
    writeSomeInts(gout, hints, tbuf, 4)
    val buff = ByteBuffer.allocate(math.min(1024*1024, 4*m.ncols*m.nrows)).order(byteOrder)
    writeSomeFloats(gout, m.data, buff, m.nrows*m.ncols)
    gout.close
  }
  
  def saveIMat(fname:String, m:IMat, compressed:Boolean=true):Unit = {
  	val gout = getOutputStream(fname, compressed)
    val hints = new Array[Int](4)
    val tbuf = ByteBuffer.allocate(16).order(byteOrder)
    hints(0) = 110 // 1=dense, 1=int
    hints(1) = m.nrows
    hints(2) = m.ncols
    hints(3) = 0
    writeSomeInts(gout, hints, tbuf, 4)
    val buff = ByteBuffer.allocate(math.min(1024*1024, 4*m.ncols*m.nrows)).order(byteOrder)
    writeSomeInts(gout, m.data, buff, m.nrows*m.ncols)
    gout.close
  }
  
  def saveDMat(fname:String, m:DMat, compressed:Boolean=true):Unit = {
    val gout = getOutputStream(fname, compressed)
    val hints = new Array[Int](4)
    val tbuf = ByteBuffer.allocate(16).order(byteOrder)
    hints(0) = 140 // 1=dense, 4=double
    hints(1) = m.nrows
    hints(2) = m.ncols
    hints(3) = 0
    writeSomeInts(gout, hints, tbuf, 4)
    val buff = ByteBuffer.allocate(math.min(1024*1024, 4*m.ncols*m.nrows)).order(byteOrder)
    writeSomeDoubles(gout, m.data, buff, m.nrows*m.ncols)
    gout.close
  }
  
  def loadSMat(fname:String, compressed:Boolean=true):SMat = {
    val gin = getInputStream(fname, compressed)
    val buff = ByteBuffer.allocate(1024*1024).order(byteOrder)
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
    val tbuf = ByteBuffer.allocate(16).order(byteOrder)
    hints(0) = 231 // 2=sparse, 3=float, 1=int
    hints(1) = m.nrows
    hints(2) = m.ncols
    hints(3) = m.nnz
    writeSomeInts(gout, hints, tbuf, 4)
    val buff = ByteBuffer.allocate(math.min(1024*1024, 4*math.max(m.ncols+1, m.nnz))).order(byteOrder)
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






