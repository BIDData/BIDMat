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
import net.jpountz.lz4._


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
  
  def readSomeBytes(din:InputStream, a:Array[Byte], n:Int) {
    var nread = 0
    while (nread < n) {
      nread += din.read(a, nread, n - nread)
    }
  }
  
  def readSomeInts(din:InputStream, a:Array[Int], buf:ByteBuffer, n:Int) {
    var nread = 0
    val ibuff = buf.asIntBuffer
    val bbuff = buf.array
    var readnow = 0
    while (nread < n) {
      val todo = if (n - nread > ibuff.capacity) ibuff.capacity else (n - nread)
      readnow += din.read(bbuff, readnow, todo*4 - readnow)
//      println("%d %d %d %d %d" format (nread, readnow, todo, ibuff.capacity, bbuff.length))
      ibuff.get(a, nread, readnow/4)
      ibuff.position(0)
      nread += readnow/4
      if (readnow % 4 != 0) {
        System.arraycopy(bbuff, 4*(readnow/4), bbuff, 0, readnow % 4)
      }
      readnow = readnow % 4
    }
  }
  
  def readSomeFloats(din:InputStream, a:Array[Float], buf:ByteBuffer, n:Int) {
    var nread = 0
    val fbuff = buf.asFloatBuffer
    val bbuff = buf.array
    var readnow = 0
    while (nread < n) {
      val todo = if (n - nread > fbuff.capacity) fbuff.capacity else (n - nread)
      readnow += din.read(bbuff, readnow, todo*4 - readnow)
      fbuff.get(a, nread, readnow/4)
      fbuff.position(0)
      nread += readnow/4  
      if (readnow % 4 != 0) {
        System.arraycopy(bbuff, 4*(readnow/4), bbuff, 0, readnow % 4)
      }
      readnow = readnow % 4
    }
  }
  
  def readSomeDoubles(din:InputStream, a:Array[Double], buf:ByteBuffer, n:Int) {
    var nread = 0
    val dbuff = buf.asDoubleBuffer
    val bbuff = buf.array
    var readnow = 0
    while (nread < n) {
      val todo = if (n - nread > dbuff.capacity) dbuff.capacity else (n - nread)
      readnow += din.read(bbuff, readnow, todo*8 - readnow)
      dbuff.get(a, nread, readnow/8)
      dbuff.position(0)
      nread += readnow/8 
      if (readnow % 8 != 0) {
        System.arraycopy(bbuff, 8*(readnow/8), bbuff, 0, readnow % 8)
      }
      readnow = readnow % 8
    }
  }
  
  def writeSomeInts(dout:OutputStream, a:Array[Int], buf:ByteBuffer, n:Int) {
    var nwritten = 0
    val ibuff = buf.asIntBuffer
    val bbuff = buf.array
    while (nwritten < n) {
    	val todo = if (n - nwritten > ibuff.capacity) ibuff.capacity else (n - nwritten)
      ibuff.put(a, nwritten, todo)
      ibuff.position(0)
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
      fbuff.position(0)
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
    	dbuff.position(0)
      dout.write(bbuff, 0, todo*8)
      nwritten += todo
    }
  }
  
  def getInputStream(fname:String, compressed:Int):InputStream = {
    val fin = new FileInputStream(fname)
    if (compressed == 2 || (compressed == 0 && fname.endsWith(".gz"))) {
      new BufferedInputStream (new GZIPInputStream(fin, 1024*1024))      
    } else  if (compressed == 3 || (compressed == 0 && fname.endsWith(".lz4"))) {
      new LZ4BlockInputStream(new BufferedInputStream(fin, 1024*1024))   
    } else {
    	new BufferedInputStream(fin, 1024*1024)
    } 
  }
  
  def getOutputStream(fname:String, compressed:Int):OutputStream = {
    if (compressed == 2 || (compressed == 0 && fname.endsWith(".gz"))) {
    	import edu.berkeley.bid.UTILS._
    	_getOutputStream(fname, Mat.compressionLevel)
    } else  if (compressed == 3 || (compressed == 0 && fname.endsWith(".lz4"))) {
      val fout = new FileOutputStream(fname)
      if (Mat.compressionLevel >= 6) {
        val hc = LZ4Factory.fastestInstance.highCompressor
        new BufferedOutputStream(new LZ4BlockOutputStream(fout, 1 << 16, hc), 1024*1024)
      } else {
      	new BufferedOutputStream(new LZ4BlockOutputStream(fout), 1024*1024)   
      }
    } else {
    	val fout = new FileOutputStream(fname)
    	new BufferedOutputStream(fout, 1024*1024)
    }
  }
  
  def loadMat(fname:String, omat:Mat, compressed:Int):Mat = {
    val gin = getInputStream(fname, compressed)
    val buff = ByteBuffer.allocate(1024).order(byteOrder)
    val hints = new Array[Int](4)
    readSomeInts(gin, hints, buff, 4)
    val ftype = hints(0)
    gin.close
    ftype match {
      case 130 => loadFMat(fname, compressed, omat)
      case 110 => loadIMat(fname, compressed, omat)
      case 140 => loadDMat(fname, compressed, omat)
      case 231 => loadSMat(fname, compressed)
      case 331 => loadSMat(fname, compressed)
      case 201 => loadBMat(fname, compressed)
      case 301 => loadBMat(fname, compressed)
    }
  }
  
  def loadMat(fname:String):Mat = loadMat(fname, null, 0)
  
  def loadMat(fname:String, omat:Mat):Mat = loadMat(fname, omat, 0)
  
  def loadFMat(fname:String, compressed:Int, omat:Mat):FMat = {
    val gin = getInputStream(fname, compressed)
    val buff = ByteBuffer.allocate(1024*1024).order(byteOrder)
    val hints = new Array[Int](4)
    readSomeInts(gin, hints, buff, 4)
    val ftype = hints(0)
    val nrows = hints(1)
    val ncols = hints(2)
    if (ftype != 130) {
      throw new RuntimeException("loadFMat expected type field 130 but was %d" format ftype)
    }
    val out = FMat.newOrCheckFMat(nrows, ncols, omat)
    readSomeFloats(gin, out.data, buff, ncols*nrows)
    gin.close
    out
  }
  
  def loadFMat(fname:String):FMat = loadFMat(fname, 0, null)
  
  def loadFMat(fname:String, compressed:Int):FMat = loadFMat(fname, compressed, null)
   
  def loadIMat(fname:String, compressed:Int, omat:Mat):IMat = {
    val gin = getInputStream(fname, compressed)
    val buff = ByteBuffer.allocate(1024*1024).order(byteOrder)
    val hints = new Array[Int](4)
    readSomeInts(gin, hints, buff, 4)
    val ftype = hints(0)
    val nrows = hints(1)
    val ncols = hints(2)
    if (ftype != 110) {
      throw new RuntimeException("loadIMat expected type field 110 but was %d" format ftype)
    }
//    println("%d %d %d\n" format (ftype, nrows, ncols))
    val out = IMat.newOrCheckIMat(nrows, ncols, omat)
    readSomeInts(gin, out.data, buff, ncols*nrows)
    gin.close
    out
  }
  
    
  def loadIMat(fname:String):IMat = loadIMat(fname, 0, null)
  
  def loadIMat(fname:String, compressed:Int):IMat = loadIMat(fname, compressed, null)
   
  def loadDMat(fname:String, compressed:Int, omat:Mat):DMat = {
    val gin = getInputStream(fname, compressed)
    val bytebuff = ByteBuffer.allocate(1024*1024).order(byteOrder)
    val hints = new Array[Int](4)
    readSomeInts(gin, hints, bytebuff, 4)
    val ftype = hints(0)
    val nrows = hints(1)
    val ncols = hints(2)
    if (ftype != 140) {
      throw new RuntimeException("loadDMat expected type field 140 but was %d" format ftype)
    }
    val out = DMat.newOrCheckDMat(nrows, ncols, omat)
    readSomeDoubles(gin, out.data, bytebuff, ncols*nrows)
    gin.close
    out
  } 
    
  def loadDMat(fname:String):DMat = loadDMat(fname, 0, null)
  
  def loadDMat(fname:String, compressed:Int):DMat = loadDMat(fname, compressed, null)
  
  def saveFMat(fname:String, m:FMat, compressed:Int=0):Unit = {
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
  
  def saveIMat(fname:String, m:IMat, compressed:Int=0):Unit = {
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
  
  def saveDMat(fname:String, m:DMat, compressed:Int=0):Unit = {
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
  
  def loadSMat(fname:String, compressed:Int=0):SMat = {
    val gin = getInputStream(fname, compressed)
    val buff = ByteBuffer.allocate(1024*1024).order(byteOrder)
    val hints = new Array[Int](4)
    readSomeInts(gin, hints, buff, 4)
    val ftype = hints(0)
    val nrows = hints(1)
    val ncols = hints(2)
    val nnz = hints(3)
    if (ftype != 231 && ftype != 331) {
      throw new RuntimeException("loadSMat expected type field 231 or 331 but was %d" format ftype)
    }
    val norows:Boolean = (ftype/100 == 3)
    val out = if (norows) {
    	SMat.SnoRows(nrows, ncols, nnz)
    } else {
      SMat(nrows, ncols, nnz)
    }
    readSomeInts(gin, out.jc, buff, ncols+1)
    if (!norows) readSomeInts(gin, out.ir, buff, nnz)
    readSomeFloats(gin, out.data, buff, nnz)
    if (Mat.ioneBased == 1) {
      MatHDF5.addOne(out.jc)
      if (!norows) MatHDF5.addOne(out.ir)
    }
    gin.close
    out
  }
  
  def loadBMat(fname:String, compressed:Int=0):BMat = {
    val gin = getInputStream(fname, compressed)
    val buff = ByteBuffer.allocate(1024*1024).order(byteOrder)
    val hints = new Array[Int](4)
    readSomeInts(gin, hints, buff, 4)
    val ftype = hints(0)
    val nrows = hints(1)
    val ncols = hints(2)
    val nnz = hints(3)
    if (ftype != 201 && ftype != 301) {
      throw new RuntimeException("loadBMat expected type field 201 or 301 but was %d" format ftype)
    }
    val norows:Boolean = (ftype/100 == 3)
    val out = if (norows) {
    	BMat.SnoRows(nrows, ncols, nnz)
    } else {
    	BMat(nrows, ncols, nnz)
    }
    readSomeInts(gin, out.jc, buff, ncols+1)
    if (!norows) readSomeInts(gin, out.ir, buff, nnz)
    readSomeBytes(gin, out.data, nnz)
    if (Mat.ioneBased == 1) {
      MatHDF5.addOne(out.jc)
      if (!norows) MatHDF5.addOne(out.ir)
    }
    gin.close
    out
  }
  
  def saveSMat(fname:String, m:SMat, compressed:Int=0):Unit = {
    val gout = getOutputStream(fname, compressed)
    val hints = new Array[Int](4)
    val tbuf = ByteBuffer.allocate(16).order(byteOrder)
    if (m.ir != null) {
    	hints(0) = 231 // 2=sparse, 3=float, 1=int
    } else {
      hints(0) = 331 // 3=sparse:norows, 3=float, 1=int
    }
    hints(1) = m.nrows
    hints(2) = m.ncols
    hints(3) = m.nnz
    writeSomeInts(gout, hints, tbuf, 4)
    val buff = ByteBuffer.allocate(math.min(1024*1024, 4*math.max(m.ncols+1, m.nnz))).order(byteOrder)
    try {
    	MatHDF5.subOne(m.jc)
    	writeSomeInts(gout, m.jc, buff, m.ncols+1)
    	if (m.ir != null) {
    	  MatHDF5.subOne(m.ir)
    	  writeSomeInts(gout, m.ir, buff, m.nnz)
    	}
    	writeSomeFloats(gout, m.data, buff, m.nnz)
    } catch {
      case e:Exception => {
      	MatHDF5.addOne(m.jc)
      	if (m.ir != null) MatHDF5.addOne(m.ir)
      	throw new RuntimeException("Exception in saveSMat "+e)
      }
      case _ => {
      	MatHDF5.addOne(m.jc)
      	if (m.ir != null) MatHDF5.addOne(m.ir)
      	throw new RuntimeException("Problem in saveSMat")
      }
    }
    MatHDF5.addOne(m.jc)
    if (m.ir != null) MatHDF5.addOne(m.ir)
    gout.close
  } 
  
   def saveBMat(fname:String, m:BMat, compressed:Int=0):Unit = {
    val gout = getOutputStream(fname, compressed)
    val hints = new Array[Int](4)
    val tbuf = ByteBuffer.allocate(16).order(byteOrder)
    if (m.ir != null) {
    	hints(0) = 201 // 2=sparse, 0=byte, 1=int
    } else {
      hints(0) = 301 // 3=sparse:norows, 0=byte, 1=int
    }
    hints(1) = m.nrows
    hints(2) = m.ncols
    hints(3) = m.nnz
    writeSomeInts(gout, hints, tbuf, 4)
    val buff = ByteBuffer.allocate(math.min(1024*1024, 4*math.max(m.ncols+1, m.nnz))).order(byteOrder)
    try {
    	MatHDF5.subOne(m.jc)
    	writeSomeInts(gout, m.jc, buff, m.ncols+1)
    	if (m.ir != null) {
    	  MatHDF5.subOne(m.ir)
    	  writeSomeInts(gout, m.ir, buff, m.nnz)
    	}
    	gout.write(m.data, 0, m.nnz)
    } catch {
      case e:Exception => {
      	MatHDF5.addOne(m.jc)
      	if (m.ir != null) MatHDF5.addOne(m.ir)
      	throw new RuntimeException("Exception in saveSMat "+e)
      }
      case _ => {
      	MatHDF5.addOne(m.jc)
      	if (m.ir != null) MatHDF5.addOne(m.ir)
      	throw new RuntimeException("Problem in saveSMat")
      }
    }
    MatHDF5.addOne(m.jc)
    if (m.ir != null) MatHDF5.addOne(m.ir)
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






