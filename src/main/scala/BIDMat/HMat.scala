package BIDMat

import java.io._
import java.util.zip._
import java.nio._
import scala.util.matching.Regex
import Regex._
import scala.collection.mutable._
import scala.concurrent.Future
import scala.concurrent.ExecutionContext.Implicits.global
import MatFunctions._
import MatHDF5._
import edu.berkeley.bid.UTILS._
import net.jpountz.lz4._


// Matrix format key
// Type = WXYZ00ABC (decimal digits)
// WXYZ = version number (currently zero)
// A = matrix type: 1 (dense), 2 (sparse), 3 (sparse, norows), 4 (1-tensor), 5 (2-tensor), 6 (3-tensor),...
// B = data type: 0 (byte), 1 (int), 2 (long), 3 (float), 4 (double), 5 (complex float), 6 (complex double)
// C = index type (sparse matrices only): 1 (int), 2 (long)

case class HMat(nr:Int, nc:Int, fileList:List[String], varname:String, blkinds:Array[Int], catdim:Int) extends Mat(nr, nc) {

  var fnameCache:String = null

  var fmatCache:Mat = null
  
  override def mytype = "HMat"
  
  override def apply(a:IMat, b:IMat) = applyx(a,b)
  
  override def apply(a:Mat, b:Mat) = applyx(a.asInstanceOf[IMat], b.asInstanceOf[IMat])
// Implement slicing from a hard disk matrix
  def applyx(a:IMat, b:IMat):Mat = { 
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
  
  var DEFAULT_BUFSIZE = 64*1024
  
  def getInputStream(fname:String, compressed:Int):InputStream = {
    val fin = new FileInputStream(fname)
    if (compressed == 2 || (compressed == 0 && fname.endsWith(".gz"))) {
      new BufferedInputStream (new GZIPInputStream(fin, DEFAULT_BUFSIZE))      
    } else  if (compressed == 3 || (compressed == 0 && fname.endsWith(".lz4"))) {
      new LZ4BlockInputStream(new BufferedInputStream(fin, DEFAULT_BUFSIZE))   
    } else {
      new BufferedInputStream(fin, DEFAULT_BUFSIZE)
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
        new BufferedOutputStream(new LZ4BlockOutputStream(fout, 1 << 16, hc), DEFAULT_BUFSIZE)
      } else {
        new BufferedOutputStream(new LZ4BlockOutputStream(fout), DEFAULT_BUFSIZE)   
      }
    } else {
      val fout = new FileOutputStream(fname)
      new BufferedOutputStream(fout, DEFAULT_BUFSIZE)
    }
  }
  
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
  
  def loadMat(fname:String, omat:Mat, compressed:Int):Mat = {
    val gin = getInputStream(fname, compressed)
    val buff = ByteBuffer.allocate(1024).order(byteOrder)
    val hints = new Array[Int](4)
    readSomeInts(gin, hints, buff, 4)
    val ftype = hints(0)
    gin.close
    ftype match {
      case 130 => loadFMat(fname, omat, compressed)
      case 110 => loadIMat(fname, omat, compressed)
      case 140 => loadDMat(fname, omat, compressed)
      case 231 => loadSMat(fname, compressed)
      case 331 => loadSMat(fname, compressed)
      case 241 => loadSDMat(fname, compressed)
      case 341 => loadSDMat(fname, compressed)
      case 201 => loadSBMat(fname, compressed)
      case 301 => loadSBMat(fname, compressed)
    }
  }
  
  def loadMat(fname:String):Mat = loadMat(fname, null, 0)
  
  def loadMat(fname:String, omat:Mat):Mat = loadMat(fname, omat, 0)
  
  def saveMat(fname:String, m:Mat, compressed:Int=0):Unit = {
    m match {
      case a:FMat => saveFMat(fname, a, compressed)
      case a:DMat => saveDMat(fname, a, compressed)
      case a:IMat => saveIMat(fname, a, compressed)
      case a:SBMat => saveSBMat(fname, a, compressed)
      case a:SDMat => saveSDMat(fname, a, compressed)
      case a:SMat => saveSMat(fname, a, compressed)
    }
  }
  
  def loadFMatTxt(fname:String, omat:Mat, compressed:Int):FMat = {
    val fin = new BufferedReader(new InputStreamReader(getInputStream(fname, compressed)))
    var nrows = 0
    var firstline = fin.readLine()
    val parts = firstline.split("[\t ,:]+")
    while (firstline != null && firstline.length > 0) {
      firstline = fin.readLine()
      nrows += 1  
    }
    fin.close
    val din = new BufferedReader(new InputStreamReader(getInputStream(fname, compressed)))
    val ncols = parts.length
    val out = FMat.newOrCheckFMat(nrows, ncols, omat)
    var irow = 0
    while (irow < nrows) {
      val parts = din.readLine().split("[\t ,:]+")
      var icol = 0
      while (icol < ncols) {
        out.data(irow + icol*out.nrows) = parts(icol).toFloat
        icol += 1
      }     
      irow += 1
    } 
    din.close
    out    
  }
  
  def loadFMat(fname:String, omat:Mat, compressed:Int):FMat = {
    if (fname.endsWith(".txt") || fname.endsWith(".txt.gz") || fname.endsWith(".txt.lz4")) {
      loadFMatTxt(fname, omat, compressed)
    } else {
      val gin = getInputStream(fname, compressed)
      val buff = ByteBuffer.allocate(DEFAULT_BUFSIZE).order(byteOrder)
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
  }
  
  def loadFMat(fname:String):FMat = loadFMat(fname, null, 0)
  
  def loadFMat(fname:String, omat:Mat):FMat = loadFMat(fname, omat, 0)
   
  def loadIMatTxt(fname:String, omat:Mat, compressed:Int):IMat = {
    val fin = new BufferedReader(new InputStreamReader(getInputStream(fname, compressed)))
    var nrows = 0
    var firstline = fin.readLine()
    val parts = firstline.split("[\t ,:]+")
    while (firstline != null && firstline.length > 0) {
      firstline = fin.readLine()
      nrows += 1  
    }
    fin.close
    val din = new BufferedReader(new InputStreamReader(getInputStream(fname, compressed)))
    val ncols = parts.length
    val out = IMat.newOrCheckIMat(nrows, ncols, omat)
    var irow = 0
    while (irow < nrows) {
      val parts = din.readLine().split("[\t ,:]+")
      var icol = 0
      while (icol < ncols) {
        out.data(irow + icol*out.nrows) = parts(icol).toInt
        icol += 1
      }     
      irow += 1
    } 
    din.close
    out    
  }
  
  def loadIMat(fname:String, omat:Mat, compressed:Int):IMat = {
    if (fname.endsWith(".txt") || fname.endsWith(".txt.gz") || fname.endsWith(".txt.lz4")) {
      loadIMatTxt(fname, omat, compressed)
    } else {
      val gin = getInputStream(fname, compressed)
      val buff = ByteBuffer.allocate(DEFAULT_BUFSIZE).order(byteOrder)
      val hints = new Array[Int](4)
      readSomeInts(gin, hints, buff, 4)
      val ftype = hints(0)
      val nrows = hints(1)
      val ncols = hints(2)
      if (ftype != 110) {
        throw new RuntimeException("loadIMat expected type field 110 but was %d" format ftype)
      }
      val out = IMat.newOrCheckIMat(nrows, ncols, omat)
      readSomeInts(gin, out.data, buff, ncols*nrows)
      gin.close
      out
    }
  }  
    
  def loadIMat(fname:String):IMat = loadIMat(fname, null, 0)
  
  def loadIMat(fname:String, omat:Mat):IMat = loadIMat(fname, omat, 0)
   
  def loadDMatTxt(fname:String, omat:Mat, compressed:Int):DMat = {
    val fin = new BufferedReader(new InputStreamReader(getInputStream(fname, compressed)))
    var nrows = 0
    var firstline = fin.readLine()
    val parts = firstline.split("[\t ,:]+")
    while (firstline != null && firstline.length > 0) {
      firstline = fin.readLine()
      nrows += 1  
    }
    fin.close
    val din = new BufferedReader(new InputStreamReader(getInputStream(fname, compressed)))
    val ncols = parts.length
    val out = DMat.newOrCheckDMat(nrows, ncols, omat)
    var irow = 0
    while (irow < nrows) {
      val parts = din.readLine().split("[\t ,:]+")
      var icol = 0
      while (icol < ncols) {
        out.data(irow + icol*out.nrows) = parts(icol).toDouble
        icol += 1
      }     
      irow += 1
    } 
    din.close
    out    
  }
  
  def loadDMat(fname:String, omat:Mat, compressed:Int):DMat = {
    if (fname.endsWith(".txt") || fname.endsWith(".txt.gz") || fname.endsWith(".txt.lz4")) {
      loadDMatTxt(fname, omat, compressed)
    } else {
      val gin = getInputStream(fname, compressed)
      val bytebuff = ByteBuffer.allocate(DEFAULT_BUFSIZE).order(byteOrder)
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
  } 
  
  def loadIDX(fname:String):FND = loadIDX(fname, 0);
  
  def loadIDX(fname:String, compressed:Int, byteOrder:ByteOrder=ByteOrder.BIG_ENDIAN):FND = {
    val gin = getInputStream(fname, compressed);
    val bytebuff = ByteBuffer.allocate(DEFAULT_BUFSIZE).order(byteOrder);
    val magicnum = new Array[Int](1);
    readSomeInts(gin, magicnum, bytebuff, 1);
    val mnum = magicnum(0);
    val mtype = (mnum >>> 8) & 0xff;
    val ndims = mnum & 0xff;
    val idims = new Array[Int](ndims);
    readSomeInts(gin, idims, bytebuff, ndims);
    val dims = idims.reverse
    val length = dims.reduce(_*_);
    val result = mtype match {
      case 0xD => {
        val out = FND(dims);
        readSomeFloats(gin, out.data, bytebuff, length);
        out;
      }      
      case 0x8 => {
        val out = FND(dims);
        val btmp = new Array[Byte](length);
        readSomeBytes(gin, btmp, length);
        var i = 0;
        while (i < length) {
          out.data(i) = (btmp(i).toInt & 0xff).toFloat;
          i += 1;
        }
        out;
      }
      case 0x9 => {
        val out = FND(dims);
        val btmp = new Array[Byte](length);
        readSomeBytes(gin, btmp, length);
        Mat.copyToFloatArray(btmp, 0, out.data, 0, length);
        out;
      }
    }
    gin.close;
    result;
  } 
    
  def loadDMat(fname:String):DMat = loadDMat(fname, null, 0)
  
  def loadDMat(fname:String, omat:Mat):DMat = loadDMat(fname, omat, 0)
  
  def saveFMat(fname:String, m:FMat, compressed:Int=0):Unit = {
    if (fname.endsWith(".txt") || fname.endsWith(".txt.gz") || fname.endsWith(".txt.lz4")) {
      saveFMatTxt(fname, m, compressed)
    } else {
      val gout = getOutputStream(fname, compressed)
      val hints = new Array[Int](4)
      val tbuf = ByteBuffer.allocate(16).order(byteOrder)
      hints(0) = 130 // 1=dense, 3=float
      hints(1) = m.nrows
      hints(2) = m.ncols
      hints(3) = 0
      writeSomeInts(gout, hints, tbuf, 4)
      val bsize = 4*m.ncols*m.nrows
      val buff = ByteBuffer.allocate(if (bsize > 0 && bsize < DEFAULT_BUFSIZE) bsize else DEFAULT_BUFSIZE).order(byteOrder)
      writeSomeFloats(gout, m.data, buff, m.nrows*m.ncols)
      gout.close
    }
  }

  def saveFND(fname:String, m:FND, compressed:Int=0):Unit = {
    val gout = getOutputStream(fname, compressed);
    val dims = m.dims;
    val ndims = dims.length;
    val hints = new Array[Int](1);
    val tbuf = ByteBuffer.allocate(dims.length*4).order(byteOrder);
    hints(0) = 30 + 100 * (ndims + 3);
    writeSomeInts(gout, hints, tbuf, 1);
    writeSomeInts(gout, dims.data, tbuf, ndims);
    hints(0) = 0;
    writeSomeInts(gout, hints, tbuf, 1);
    val bsize = 4*m.length;
    val buff = ByteBuffer.allocate(if (bsize > 0 && bsize < DEFAULT_BUFSIZE) bsize else DEFAULT_BUFSIZE).order(byteOrder);
    writeSomeFloats(gout, m.data, buff, m.length);
    gout.close;
  }

  def loadFND(fname:String, omat:ND, compressed:Int):FND = {
    val gin = getInputStream(fname, compressed);
    val buff = ByteBuffer.allocate(DEFAULT_BUFSIZE).order(byteOrder);
    val hints = new Array[Int](1);
    readSomeInts(gin, hints, buff, 1);
    if (hints(0) % 100 != 30) throw new RuntimeException("loadFND: bad type code " + hints(0));
    val ndims = hints(0) / 100 - 3;
    val dims = new Array[Int](ndims);
    readSomeInts(gin, dims, buff, ndims);
    readSomeInts(gin, hints, buff, 1);
    val out = FND.newOrCheckFND(dims, omat);
    readSomeFloats(gin, out.data, buff, out.length);
    gin.close;
    out;
  }
  
  def saveFMatTxt(fname:String, m:FMat, compressed:Int=0, delim:String="\t"):Unit = {
    val gout = getOutputStream(fname, compressed)
    val fout = new BufferedWriter(new OutputStreamWriter(gout))
    var i = 0
    while (i < m.nrows) {
      if (m.ncols > 0) {
        fout.write(m(i,0).toString)
      }
      var j = 1
      while (j < m.ncols) {
        fout.write(delim + m(i,j).toString)
        j += 1
      }
      fout.write("\n")
      i += 1
    }
    fout.close
  }
  
  def saveIMat(fname:String, m:IMat, compressed:Int=0):Unit = {
    if (fname.endsWith(".txt") || fname.endsWith(".txt.gz") || fname.endsWith(".txt.lz4")) {
      saveIMatTxt(fname, m, compressed)
    } else {
      val gout = getOutputStream(fname, compressed)
      val hints = new Array[Int](4)
      val tbuf = ByteBuffer.allocate(16).order(byteOrder)
      hints(0) = 110 // 1=dense, 1=int
      hints(1) = m.nrows
      hints(2) = m.ncols
      hints(3) = 0
      writeSomeInts(gout, hints, tbuf, 4)
      val bsize = 4*m.ncols*m.nrows
      val buff = ByteBuffer.allocate(if (bsize > 0 && bsize < DEFAULT_BUFSIZE) bsize else DEFAULT_BUFSIZE).order(byteOrder)
      writeSomeInts(gout, m.data, buff, m.nrows*m.ncols)
      gout.close
    }
  }
  
  def saveIMatTxt(fname:String, m:IMat, compressed:Int=0, delim:String="\t"):Unit = {
    val gout = getOutputStream(fname, compressed)
    val fout = new BufferedWriter(new OutputStreamWriter(gout))
    var i = 0
    while (i < m.nrows) {
      if (m.ncols > 0) {
        fout.write(m(i,0).toString)
      }
      var j = 1
      while (j < m.ncols) {
        fout.write(delim + m(i,j).toString)
        j += 1
      }
      fout.write("\n")
      i += 1
    }
    fout.close
  }
  
  def saveDMat(fname:String, m:DMat, compressed:Int=0):Unit = {
    if (fname.endsWith(".txt") || fname.endsWith(".txt.gz") || fname.endsWith(".txt.lz4")) {
      saveDMatTxt(fname, m, compressed)
    } else {
      val gout = getOutputStream(fname, compressed)
      val hints = new Array[Int](4)
      val tbuf = ByteBuffer.allocate(16).order(byteOrder)
      hints(0) = 140 // 1=dense, 4=double
      hints(1) = m.nrows
      hints(2) = m.ncols
      hints(3) = 0
      writeSomeInts(gout, hints, tbuf, 4)
      val bsize = 8*m.ncols*m.nrows
      val buff = ByteBuffer.allocate(if (bsize > 0 && bsize < DEFAULT_BUFSIZE) bsize else DEFAULT_BUFSIZE).order(byteOrder)
      writeSomeDoubles(gout, m.data, buff, m.nrows*m.ncols)
      gout.close
    }
  }
  
  def saveDMatTxt(fname:String, m:DMat, compressed:Int=0, delim:String="\t"):Unit = {
    val gout = getOutputStream(fname, compressed)
    val fout = new BufferedWriter(new OutputStreamWriter(gout))
    var i = 0
    while (i < m.nrows) {
      if (m.ncols > 0) {
        fout.write(m(i,0).toString)
      }
      var j = 1
      while (j < m.ncols) {
        fout.write(delim + m(i,j).toString)
        j += 1
      }
      fout.write("\n")
      i += 1
    }
    fout.close
  }
  
  def loadSMat(fname:String, compressed:Int=0):SMat = {
    val gin = getInputStream(fname, compressed)
    val buff = ByteBuffer.allocate(DEFAULT_BUFSIZE).order(byteOrder)
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
  
  def loadSDMat(fname:String, compressed:Int=0):SDMat = {
    val gin = getInputStream(fname, compressed)
    val buff = ByteBuffer.allocate(DEFAULT_BUFSIZE).order(byteOrder)
    val hints = new Array[Int](4)
    readSomeInts(gin, hints, buff, 4)
    val ftype = hints(0)
    val nrows = hints(1)
    val ncols = hints(2)
    val nnz = hints(3)
    if (ftype != 241 && ftype != 341) {
      throw new RuntimeException("loadSDMat expected type field 241 or 341 but was %d" format ftype)
    }
    val norows:Boolean = (ftype/100 == 3)
    val out = if (norows) {
        SDMat.SnoRows(nrows, ncols, nnz)
    } else {
      SDMat(nrows, ncols, nnz)
    }
    readSomeInts(gin, out.jc, buff, ncols+1)
    if (!norows) readSomeInts(gin, out.ir, buff, nnz)
    readSomeDoubles(gin, out.data, buff, nnz)
    if (Mat.ioneBased == 1) {
      MatHDF5.addOne(out.jc)
      if (!norows) MatHDF5.addOne(out.ir)
    }
    gin.close
    out
  }
  
  def loadSBMat(fname:String, compressed:Int=0):SBMat = {
    val gin = getInputStream(fname, compressed)
    val buff = ByteBuffer.allocate(DEFAULT_BUFSIZE).order(byteOrder)
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
    	SBMat.SnoRows(nrows, ncols, nnz)
    } else {
    	SBMat(nrows, ncols, nnz)
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
    val buff = ByteBuffer.allocate(4*math.min(DEFAULT_BUFSIZE/4, math.max(m.ncols+1, m.nnz))).order(byteOrder)
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
      case _:Throwable => {
      	MatHDF5.addOne(m.jc)
      	if (m.ir != null) MatHDF5.addOne(m.ir)
      	throw new RuntimeException("Problem in saveSMat")
      }
    }
    MatHDF5.addOne(m.jc)
    if (m.ir != null) MatHDF5.addOne(m.ir)
    gout.close
  } 
  
   def saveSDMat(fname:String, m:SDMat, compressed:Int=0):Unit = {
    val gout = getOutputStream(fname, compressed)
    val hints = new Array[Int](4)
    val tbuf = ByteBuffer.allocate(16).order(byteOrder)
    if (m.ir != null) {
        hints(0) = 241 // 2=sparse, 4=double, 1=int
    } else {
      hints(0) = 341 // 3=sparse:norows, 4=double, 1=int
    }
    hints(1) = m.nrows
    hints(2) = m.ncols
    hints(3) = m.nnz
    writeSomeInts(gout, hints, tbuf, 4)
    val buff = ByteBuffer.allocate(4*math.min(DEFAULT_BUFSIZE/4, math.max(m.ncols+1, m.nnz))).order(byteOrder)
    try {
        MatHDF5.subOne(m.jc)
        writeSomeInts(gout, m.jc, buff, m.ncols+1)
        if (m.ir != null) {
          MatHDF5.subOne(m.ir)
          writeSomeInts(gout, m.ir, buff, m.nnz)
        }
        writeSomeDoubles(gout, m.data, buff, m.nnz)
    } catch {
      case e:Exception => {
        MatHDF5.addOne(m.jc)
        if (m.ir != null) MatHDF5.addOne(m.ir)
        throw new RuntimeException("Exception in saveSDMat "+e)
      }
      case _:Throwable => {
        MatHDF5.addOne(m.jc)
        if (m.ir != null) MatHDF5.addOne(m.ir)
        throw new RuntimeException("Problem in saveSDMat")
      }
    }
    MatHDF5.addOne(m.jc)
    if (m.ir != null) MatHDF5.addOne(m.ir)
    gout.close
  } 
  
   def saveSBMat(fname:String, m:SBMat, compressed:Int=0):Unit = {
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
    val buff = ByteBuffer.allocate(4*math.min(DEFAULT_BUFSIZE/4, math.max(m.ncols+1, m.nnz))).order(byteOrder)
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
      case _:Throwable => {
      	MatHDF5.addOne(m.jc)
      	if (m.ir != null) MatHDF5.addOne(m.ir)
      	throw new RuntimeException("Problem in saveSMat")
      }
    }
    MatHDF5.addOne(m.jc)
    if (m.ir != null) MatHDF5.addOne(m.ir)
    gout.close
  } 
   
  /** Load a file in ASCII LibSVM format */
  /* Outputs a data matrix first, and then a matrix c of cat labels where c(i) = label of instance i. */
  /* data is an SMat, labels are integer, weights are FMat. */
  def loadLibSVM(fname:String, nrows:Int, compressed:Int = 0):(SMat, IMat, FMat) = {
    var fin = new BufferedReader(new InputStreamReader (getInputStream(fname, compressed)));
    var firstline = fin.readLine();
    var parts = if (firstline != null) firstline.split("[\t ]+") else null;
    var ncols = 0;
    var nnz = 0;
    var maxindx = 0;
    while (firstline != null && parts.length > 0) {
      nnz += parts.length - 1;
      firstline = fin.readLine();
      if (firstline != null) parts = firstline.split("[\t ]+");
      ncols += 1;
    }
    var datamat = SMat(nrows, ncols, nnz);
    var cmat = IMat(1, ncols);
    var wmat = ones(1, ncols);
    fin.close();
    fin = new BufferedReader(new InputStreamReader (getInputStream(fname, compressed)));
    firstline = fin.readLine();
    if (firstline != null) parts = firstline.split("[\t ]+");
    ncols = 0;
    nnz = 0;
    val ioneBased = Mat.ioneBased;
    datamat.jc(0) = ioneBased; 
    while (firstline != null && parts.length > 0) {
      val parts0 = parts(0).split(":");
      cmat.data(ncols) = parts0(0).toInt;
      if (parts0.length > 1) wmat.data(ncols) = parts0(1).toFloat;
      var i = 1;
      while (i < parts.length) {
        val pair = parts(i).split(":");
        val indx = pair(0).toInt;
        val vval = pair(1).toFloat;
        datamat.ir(nnz) = indx + ioneBased;
        datamat.data(nnz) = vval;
        i += 1;
        nnz += 1;
      }
      firstline = fin.readLine();
      if (firstline != null) parts = firstline.split("[\t ]+");
      ncols += 1;
      datamat.jc(ncols) = nnz + ioneBased;
    }
    fin.close();
    (datamat, cmat, wmat);    
  }
  
  def saveLibSVM(fname:String, data:SMat, labels:SMat):Unit = saveLibSVM(fname, data, labels, null, 0);
  
  def saveLibSVM(fname:String, data:SMat, labels:SMat, weights:FMat):Unit = saveLibSVM(fname, data, labels, weights, 0);
   
  def saveLibSVM(fname:String, data:SMat, labels:SMat, weights:FMat, compressed:Int):Unit = {
    val fout = new BufferedWriter (new OutputStreamWriter (getOutputStream(fname, compressed)));
    val (ilab, jlab) = find2(labels);
    val jc = data.jc;
    val ioffset = Mat.ioneBased;
    var i = 0;
    while (i < ilab.length) {
      fout.write(ilab.data(i).toString);
      val icol = jlab.data(i);
      if (weights.asInstanceOf[AnyRef] != null) {
        fout.write(":"+weights.data(icol))
      } 
      var j = jc(icol) - ioffset
      while (j < jc(icol+1) - ioffset) {
        fout.write(" "+(data.ir(j)-ioffset).toString + ":" + data.data(j).toString);
        j += 1
      }
      fout.write("\n")
      i += 1
    }
    fout.close()
  }

  def testLoad(fname:String, varname:String, n:Int) = {
    val a = new Array[SMat](n)
    var ndone = izeros(n,1)
    for (i <- 0 until n) {
      Future {
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
      val (nr, nc) = MatHDF5.readMat2Dims(fn, varname)
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






