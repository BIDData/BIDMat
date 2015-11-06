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

/**
 *
 * Matrix format key
 * Type = WXYZ00ABC (decimal digits)
 * WXYZ = version number (currently zero)
 * A = matrix type: 1 (dense), 2 (sparse), 3 (sparse, norows), 4 (1-tensor), 5 (2-tensor), 6 (3-tensor),...
 * B = data type: 0 (byte), 1 (int), 2 (long), 3 (float), 4 (double), 5 (complex float), 6 (complex double)
 * C = index type (sparse matrices only): 1 (int), 2 (long)
 *
 */

/**
 * Notes: Matrix I/O routines are now based on DataInput and DataOutput classes to provide direct support for HDFS.
 */

case class HMat(nr: Int, nc: Int, fileList: List[String], varname: String, blkinds: Array[Int], catdim: Int) extends Mat(nr, nc) {

  var fnameCache: String = null

  var fmatCache: Mat = null

  override def mytype = "HMat"

  override def apply(a: IMat, b: IMat) = applyx(a, b)

  override def apply(a: Mat, b: Mat) = applyx(a.asInstanceOf[IMat], b.asInstanceOf[IMat])
  // Implement slicing from a hard disk matrix
  def applyx(a: IMat, b: IMat): Mat = {
    var ilast: Int = 0
    def findindx(ind: Int): Int = {
      while (ilast >= 0 && ind < blkinds(ilast)) ilast -= 1
      while (ilast < blkinds.length && ind >= blkinds(ilast)) ilast += 1
      if (ilast >= 0 && blkinds(ilast) <= ind && blkinds(ilast + 1) > ind) {
        ilast
      } else {
        -1
      }
    }

    val locs = IMat(1, b.length)
    var i = 0
    var iblk = 0
    var out: Mat = null
    while (i <= b.length) {
      if (i < b.length) locs(i) = findindx(b(i))
      if (i == b.length || locs(i) != locs(iblk)) {
        if (fnameCache == null || fileList(locs(iblk)) != fnameCache) {
          fmatCache = MatHDF5.hload(fileList(locs(iblk)), varname).asInstanceOf[Mat]
          fnameCache = fileList(locs(iblk))
        }
        val newmat = fmatCache(a, b(MatFunctions.irow(iblk -> i)))
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

  var DEFAULT_BUFSIZE = 64 * 1024;

  val classLoader = HMat.getClass.getClassLoader();
  var HDFSIOclass: Class[_] = null;
  var HDFSIOinstance: AnyRef = null;
  var writeMatMethod: java.lang.reflect.Method = null;
  var writeNDMethod: java.lang.reflect.Method = null;
  var readMatMethod: java.lang.reflect.Method = null;
  var readNDMethod: java.lang.reflect.Method = null;

  var highCompressor: AnyRef = null;
  var fastCompressor: AnyRef = null;

  def initLZ4 = {
    if (highCompressor == null) {
      highCompressor = if (Mat.ostype == Mat.OS_ANDROID) {
        LZ4Factory.safeInstance.highCompressor
      } else {
        LZ4Factory.fastestInstance.highCompressor
      }
      fastCompressor = if (Mat.ostype == Mat.OS_ANDROID) {
        LZ4Factory.safeInstance.fastCompressor
      } else {
        LZ4Factory.fastestInstance.fastCompressor
      }
    }
  }

  def checkHDFSloaded = {
    if (HDFSIOclass == null) {
      HDFSIOclass = classLoader.loadClass("BIDMat.HDFSIO");
      HDFSIOinstance = HDFSIOclass.newInstance().asInstanceOf[AnyRef];
      writeMatMethod = HDFSIOclass.getMethod("writeMat", classOf[String], classOf[Mat], classOf[java.lang.Integer]);
      writeNDMethod = HDFSIOclass.getMethod("writeND", classOf[String], classOf[ND], classOf[java.lang.Integer]);
      readMatMethod = HDFSIOclass.getMethod("readMat", classOf[String], classOf[Mat]);
      readNDMethod = HDFSIOclass.getMethod("readND", classOf[String], classOf[ND]);
    }
  }

  def HDFSwriteMat(fname: String, m: Mat, compress: Int): Unit = {
    checkHDFSloaded;
    writeMatMethod.invoke(HDFSIOinstance, fname, m, new java.lang.Integer(compress));
  }

  def HDFSwriteND(fname: String, m: ND, compress: Int): Unit = {
    checkHDFSloaded;
    writeNDMethod.invoke(HDFSIOinstance, fname, m, new java.lang.Integer(compress));
  }

  def HDFSreadMat(fname: String, omat: Mat): Mat = {
    checkHDFSloaded;
    readMatMethod.invoke(HDFSIOinstance, fname, omat).asInstanceOf[Mat];
  }

  def HDFSreadFND(fname: String, ond: ND): ND = {
    checkHDFSloaded;
    readNDMethod.invoke(HDFSIOinstance, fname, ond).asInstanceOf[ND];
  }

  def getInputStream(fname: String, compressed: Int): DataInputStream = {
    val fin = new FileInputStream(fname)
    if (compressed == 2 || (compressed == 0 && fname.endsWith(".gz"))) {
      new DataInputStream(new BufferedInputStream(new GZIPInputStream(fin, DEFAULT_BUFSIZE)));
    } else if (compressed == 3 || (compressed == 0 && fname.endsWith(".lz4"))) {
      new DataInputStream(new LZ4BlockInputStream(new BufferedInputStream(fin, DEFAULT_BUFSIZE)));
    } else {
      new DataInputStream(new BufferedInputStream(fin, DEFAULT_BUFSIZE));
    }
  }

  def getOutputStream(fname: String, compressed: Int): DataOutputStream = {
    if (compressed == 2 || (compressed == 0 && fname.endsWith(".gz"))) {
      import edu.berkeley.bid.UTILS._
      new DataOutputStream(_getOutputStream(fname, Mat.compressionLevel));
    } else if (compressed == 3 || (compressed == 0 && fname.endsWith(".lz4"))) {
      val fout = new FileOutputStream(fname)
      initLZ4
      if (Mat.compressionLevel >= 6) {
        new DataOutputStream(new BufferedOutputStream(new LZ4BlockOutputStream(fout, 1 << 16, highCompressor.asInstanceOf[net.jpountz.lz4.LZ4Compressor]), DEFAULT_BUFSIZE));
      } else {
        new DataOutputStream(new BufferedOutputStream(new LZ4BlockOutputStream(fout, 1 << 16, fastCompressor.asInstanceOf[net.jpountz.lz4.LZ4Compressor]), DEFAULT_BUFSIZE));
      }
    } else {
      val fout = new FileOutputStream(fname)
      new DataOutputStream(new BufferedOutputStream(fout, DEFAULT_BUFSIZE));
    }
  }

  def readSomeBytes(din: DataInput, a: Array[Byte], n: Int) {
    din.readFully(a, 0, n);
  }

  def readSomeInts(din: DataInput, a: Array[Int], buf: ByteBuffer, n: Int) {
    val lbytes = 2;
    val mask = (1 << lbytes) - 1;
    var nread = 0;
    val ibuff = buf.asIntBuffer;
    val bbuff = buf.array;
    while (nread < n) {
      val todo = math.min(ibuff.capacity, n - nread);
      din.readFully(bbuff, 0, (todo << lbytes));
      ibuff.get(a, nread, todo);
      ibuff.position(0)
      nread += todo;
    }
  }

  def readSomeLongs(din: DataInput, a: Array[Long], buf: ByteBuffer, n: Int) {
    val lbytes = 3;
    var nread = 0;
    val lbuff = buf.asLongBuffer;
    val bbuff = buf.array;
    while (nread < n) {
      val todo = math.min(lbuff.capacity, n - nread);
      din.readFully(bbuff, 0, (todo << lbytes));
      lbuff.get(a, nread, todo);
      lbuff.position(0);
      nread += todo;
    }
  }

  def readSomeFloats(din: DataInput, a: Array[Float], buf: ByteBuffer, n: Int) {
    val lbytes = 2;
    var nread = 0;
    val fbuff = buf.asFloatBuffer;
    val bbuff = buf.array;
    while (nread < n) {
      val todo = math.min(fbuff.capacity, n - nread);
      din.readFully(bbuff, 0, (todo << lbytes));
      fbuff.get(a, nread, todo);
      fbuff.position(0);
      nread += todo;
    }
  }

  def readSomeDoubles(din: DataInput, a: Array[Double], buf: ByteBuffer, n: Int) {
    val lbytes = 3;
    var nread = 0;
    val dbuff = buf.asDoubleBuffer;
    val bbuff = buf.array;
    while (nread < n) {
      val todo = math.min(dbuff.capacity, n - nread);
      din.readFully(bbuff, 0, (todo << lbytes));
      dbuff.get(a, nread, todo);
      dbuff.position(0);
      nread += todo;
    }
  }

  def writeSomeInts(dout: DataOutput, a: Array[Int], buf: ByteBuffer, n: Int) {
    var nwritten = 0
    val ibuff = buf.asIntBuffer
    val bbuff = buf.array
    while (nwritten < n) {
      val todo = if (n - nwritten > ibuff.capacity) ibuff.capacity else (n - nwritten)
      ibuff.put(a, nwritten, todo)
      ibuff.position(0)
      dout.write(bbuff, 0, todo * 4)
      nwritten += todo
    }
  }

  def writeSomeLongs(dout: DataOutput, a: Array[Long], buf: ByteBuffer, n: Int) {
    var nwritten = 0
    val ibuff = buf.asLongBuffer
    val bbuff = buf.array
    while (nwritten < n) {
      val todo = if (n - nwritten > ibuff.capacity) ibuff.capacity else (n - nwritten)
      ibuff.put(a, nwritten, todo)
      ibuff.position(0)
      dout.write(bbuff, 0, todo * 8)
      nwritten += todo
    }
  }

  def writeSomeFloats(dout: DataOutput, a: Array[Float], buf: ByteBuffer, n: Int) {
    var nwritten = 0
    val fbuff = buf.asFloatBuffer
    val bbuff = buf.array
    while (nwritten < n) {
      val todo = if (n - nwritten > fbuff.capacity) fbuff.capacity else (n - nwritten)
      fbuff.put(a, nwritten, todo)
      fbuff.position(0)
      dout.write(bbuff, 0, todo * 4)
      nwritten += todo
    }
  }

  def writeSomeDoubles(dout: DataOutput, a: Array[Double], buf: ByteBuffer, n: Int) {
    var nwritten = 0
    val dbuff = buf.asDoubleBuffer
    val bbuff = buf.array
    while (nwritten < n) {
      val todo = if (n - nwritten > dbuff.capacity) dbuff.capacity else (n - nwritten)
      dbuff.put(a, nwritten, todo)
      dbuff.position(0)
      dout.write(bbuff, 0, todo * 8)
      nwritten += todo
    }
  }

  def loadMat(fname: String, omat: Mat, compressed: Int): Mat = {
    val gin = getInputStream(fname, compressed)
    val buff = ByteBuffer.allocate(1024).order(byteOrder)
    val hints = new Array[Int](4)
    readSomeInts(gin, hints, buff, 4)
    val ftype = hints(0)
    gin.close
    ftype match {
      case 130 => loadFMat(fname, omat, compressed)
      case 110 => loadIMat(fname, omat, compressed)
      case 120 => loadLMat(fname, omat, compressed)
      case 140 => loadDMat(fname, omat, compressed)
      case 231 => loadSMat(fname, compressed)
      case 331 => loadSMat(fname, compressed)
      case 241 => loadSDMat(fname, compressed)
      case 341 => loadSDMat(fname, compressed)
      case 201 => loadSBMat(fname, compressed)
      case 202 => loadCSMat(fname, compressed)
      case 301 => loadSBMat(fname, compressed)
      case 302 => loadCSMat(fname, compressed)
    }
  }

  def loadMat(fname: String): Mat = loadMat(fname, null, 0)

  def loadMat(fname: String, omat: Mat): Mat = loadMat(fname, omat, 0)

  def saveMat(fname: String, m: Mat, compressed: Int = 0): Unit = {
    m match {
      case a: FMat => saveFMat(fname, a, compressed)
      case a: DMat => saveDMat(fname, a, compressed)
      case a: IMat => saveIMat(fname, a, compressed)
      case a: LMat => saveLMat(fname, a, compressed)
      case a: SBMat => saveSBMat(fname, a, compressed)
      case a: SDMat => saveSDMat(fname, a, compressed)
      case a: SMat => saveSMat(fname, a, compressed)
      case a: CSMat => saveCSMat(fname, a, compressed)
    }
  }

  def loadFMatTxt(fname: String, omat: Mat, compressed: Int): FMat = {
    val fin = new BufferedReader(new InputStreamReader(getInputStream(fname, compressed).asInstanceOf[DataInputStream]))
    var nrows = 0
    var firstline = fin.readLine()
    val parts = firstline.split("[\t ,:]+")
    while (firstline != null && firstline.length > 0) {
      firstline = fin.readLine()
      nrows += 1
    }
    fin.close
    val din = new BufferedReader(new InputStreamReader(getInputStream(fname, compressed).asInstanceOf[DataInputStream]))
    val ncols = parts.length
    val out = FMat.newOrCheckFMat(nrows, ncols, omat)
    var irow = 0
    while (irow < nrows) {
      val parts = din.readLine().split("[\t ,:]+")
      var icol = 0
      while (icol < ncols) {
        out.data(irow + icol * out.nrows) = parts(icol).toFloat
        icol += 1
      }
      irow += 1
    }
    din.close
    out
  }

  def loadFMat(fname: String, omat: Mat, compressed: Int): FMat = {
    if (fname.startsWith("hdfs:")) {
      HDFSreadMat(fname, omat).asInstanceOf[FMat];
    } else if (fname.endsWith(".txt") || fname.endsWith(".txt.gz") || fname.endsWith(".txt.lz4")) {
      loadFMatTxt(fname, omat, compressed)
    } else {
      val gin = getInputStream(fname, compressed);
      val out = loadFMat(gin, omat);
      gin.close;
      out
    }
  }

  def loadFMat(gin: DataInput, omat: Mat): FMat = {
    val buff = ByteBuffer.allocate(DEFAULT_BUFSIZE).order(byteOrder);
    val hints = new Array[Int](4);
    readSomeInts(gin, hints, buff, 4);
    val ftype = hints(0);
    val nrows = hints(1);
    val ncols = hints(2);
    if (ftype != 130) {
      throw new RuntimeException("loadFMat expected type field 130 but was %d" format ftype);
    }
    val out = FMat.newOrCheckFMat(nrows, ncols, omat);
    readSomeFloats(gin, out.data, buff, ncols * nrows);
    out;
  }

  def loadFMat(fname: String): FMat = loadFMat(fname, null, 0)

  def loadFMat(fname: String, omat: Mat): FMat = loadFMat(fname, omat, 0)

  def loadFMat(gin: DataInput): FMat = loadFMat(gin, null)

  def loadIMatTxt(fname: String, omat: Mat, compressed: Int): IMat = {
    val fin = new BufferedReader(new InputStreamReader(getInputStream(fname, compressed).asInstanceOf[DataInputStream]))
    var nrows = 0
    var firstline = fin.readLine()
    val parts = firstline.split("[\t ,:]+")
    while (firstline != null && firstline.length > 0) {
      firstline = fin.readLine()
      nrows += 1
    }
    fin.close
    val din = new BufferedReader(new InputStreamReader(getInputStream(fname, compressed).asInstanceOf[DataInputStream]))
    val ncols = parts.length
    val out = IMat.newOrCheckIMat(nrows, ncols, omat)
    var irow = 0
    while (irow < nrows) {
      val parts = din.readLine().split("[\t ,:]+")
      var icol = 0
      while (icol < ncols) {
        out.data(irow + icol * out.nrows) = parts(icol).toInt
        icol += 1
      }
      irow += 1
    }
    din.close
    out
  }

  def loadIMat(fname: String, omat: Mat, compressed: Int): IMat = {
    if (fname.startsWith("hdfs:")) {
      HDFSreadMat(fname, omat).asInstanceOf[IMat];
    } else if (fname.endsWith(".txt") || fname.endsWith(".txt.gz") || fname.endsWith(".txt.lz4")) {
      loadIMatTxt(fname, omat, compressed)
    } else {
      val gin = getInputStream(fname, compressed)
      val out = loadIMat(gin, omat)
      gin.close
      out
    }
  }

  def loadIMat(gin: DataInput, omat: Mat): IMat = {
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
    readSomeInts(gin, out.data, buff, ncols * nrows)
    out
  }

  def loadIMat(fname: String): IMat = loadIMat(fname, null, 0)

  def loadIMat(fname: String, omat: Mat): IMat = loadIMat(fname, omat, 0)

  def loadLMatTxt(fname: String, omat: Mat, compressed: Int): LMat = {
    val fin = new BufferedReader(new InputStreamReader(getInputStream(fname, compressed).asInstanceOf[DataInputStream]))
    var nrows = 0
    var firstline = fin.readLine()
    val parts = firstline.split("[\t ,:]+")
    while (firstline != null && firstline.length > 0) {
      firstline = fin.readLine()
      nrows += 1
    }
    fin.close
    val din = new BufferedReader(new InputStreamReader(getInputStream(fname, compressed).asInstanceOf[DataInputStream]))
    val ncols = parts.length
    val out = LMat.newOrCheckLMat(nrows, ncols, omat)
    var irow = 0
    while (irow < nrows) {
      val parts = din.readLine().split("[\t ,:]+")
      var icol = 0
      while (icol < ncols) {
        out.data(irow + icol * out.nrows) = parts(icol).toLong
        icol += 1
      }
      irow += 1
    }
    din.close
    out
  }

  def loadLMat(fname: String, omat: Mat, compressed: Int): LMat = {
    if (fname.startsWith("hdfs:")) {
      HDFSreadMat(fname, omat).asInstanceOf[LMat];
    } else if (fname.endsWith(".txt") || fname.endsWith(".txt.gz") || fname.endsWith(".txt.lz4")) {
      loadLMatTxt(fname, omat, compressed)
    } else {
      val gin = getInputStream(fname, compressed)
      val out = loadLMat(gin, omat)
      gin.close
      out
    }
  }

  def loadLMat(gin: DataInput, omat: Mat): LMat = {
    val buff = ByteBuffer.allocate(DEFAULT_BUFSIZE).order(byteOrder)
    val hints = new Array[Int](4)
    readSomeInts(gin, hints, buff, 4)
    val ftype = hints(0)
    val nrows = hints(1)
    val ncols = hints(2)
    if (ftype != 120) {
      throw new RuntimeException("loadLMat expected type field 120 but was %d" format ftype)
    }
    val out = LMat.newOrCheckLMat(nrows, ncols, omat)
    readSomeLongs(gin, out.data, buff, ncols * nrows)
    out
  }

  def loadLMat(fname: String): LMat = loadLMat(fname, null, 0)

  def loadLMat(fname: String, omat: Mat): LMat = loadLMat(fname, omat, 0)

  def loadDMatTxt(fname: String, omat: Mat, compressed: Int): DMat = {
    val fin = new BufferedReader(new InputStreamReader(getInputStream(fname, compressed).asInstanceOf[DataInputStream]))
    var nrows = 0
    var firstline = fin.readLine()
    val parts = firstline.split("[\t ,:]+")
    while (firstline != null && firstline.length > 0) {
      firstline = fin.readLine()
      nrows += 1
    }
    fin.close
    val din = new BufferedReader(new InputStreamReader(getInputStream(fname, compressed).asInstanceOf[DataInputStream]))
    val ncols = parts.length
    val out = DMat.newOrCheckDMat(nrows, ncols, omat)
    var irow = 0
    while (irow < nrows) {
      val parts = din.readLine().split("[\t ,:]+")
      var icol = 0
      while (icol < ncols) {
        out.data(irow + icol * out.nrows) = parts(icol).toDouble
        icol += 1
      }
      irow += 1
    }
    din.close
    out
  }

  def loadDMat(fname: String, omat: Mat, compressed: Int): DMat = {
    if (fname.startsWith("hdfs:")) {
      HDFSreadMat(fname, omat).asInstanceOf[DMat];
    } else if (fname.endsWith(".txt") || fname.endsWith(".txt.gz") || fname.endsWith(".txt.lz4")) {
      loadDMatTxt(fname, omat, compressed)
    } else {
      val gin = getInputStream(fname, compressed)
      val out = loadDMat(gin, omat)
      gin.close
      out
    }
  }

  def loadDMat(gin: DataInput, omat: Mat): DMat = {
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
    readSomeDoubles(gin, out.data, bytebuff, ncols * nrows)
    out
  }

  def loadIDX(fname: String): FND = loadIDX(fname, 0);

  def loadIDX(fname: String, compressed: Int, byteOrder: ByteOrder = ByteOrder.BIG_ENDIAN): FND = {
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
    val length = dims.reduce(_ * _);
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

  def loadDMat(fname: String): DMat = loadDMat(fname, null, 0)

  def loadDMat(fname: String, omat: Mat): DMat = loadDMat(fname, omat, 0)

  def saveFMat(fname: String, m: FMat, compressed: Int = 0): Unit = {
    if (fname.startsWith("hdfs:")) {
      HDFSwriteMat(fname, m, compressed);
    } else if (fname.endsWith(".txt") || fname.endsWith(".txt.gz") || fname.endsWith(".txt.lz4")) {
      saveFMatTxt(fname, m, compressed)
    } else {
      val gout = getOutputStream(fname, compressed);
      saveFMat(gout, m);
      gout.close;
    }
  }

  def saveFMat(gout: DataOutput, m: FMat) = {
    val hints = new Array[Int](4);
    val tbuf = ByteBuffer.allocate(16).order(byteOrder);
    hints(0) = 130; // 1=dense, 3=float
    hints(1) = m.nrows;
    hints(2) = m.ncols;
    hints(3) = 0;
    writeSomeInts(gout, hints, tbuf, 4);
    val bsize = 4 * m.ncols * m.nrows;
    val buff = ByteBuffer.allocate(if (bsize > 0 && bsize < DEFAULT_BUFSIZE) bsize else DEFAULT_BUFSIZE).order(byteOrder);
    writeSomeFloats(gout, m.data, buff, m.nrows * m.ncols);
  }

  def saveFND(fname: String, m: FND, compressed: Int = 0): Unit = {
    if (fname.startsWith("hdfs:")) {
      HDFSwriteND(fname, m, compressed)
    } else {
      val gout = getOutputStream(fname, compressed);
      saveFND(gout, m);
      gout.close;
    }
  }

  def saveFND(gout: DataOutput, m: FND): Unit = {
    val dims = m.dims;
    val ndims = dims.length;
    val hints = new Array[Int](1);
    val tbuf = ByteBuffer.allocate(dims.length * 4).order(byteOrder);
    hints(0) = 30 + 100 * (ndims + 3);
    writeSomeInts(gout, hints, tbuf, 1);
    writeSomeInts(gout, dims.data, tbuf, ndims);
    hints(0) = 0;
    writeSomeInts(gout, hints, tbuf, 1);
    val bsize = 4 * m.length;
    val buff = ByteBuffer.allocate(if (bsize > 0 && bsize < DEFAULT_BUFSIZE) bsize else DEFAULT_BUFSIZE).order(byteOrder);
    writeSomeFloats(gout, m.data, buff, m.length);
  }

  def loadFND(fname: String, omat: ND, compressed: Int): FND = {
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

  def saveFMatTxt(fname: String, m: FMat, compressed: Int = 0, delim: String = "\t"): Unit = {
    val gout = getOutputStream(fname, compressed)
    val fout = new BufferedWriter(new OutputStreamWriter(gout.asInstanceOf[DataOutputStream]))
    var i = 0
    while (i < m.nrows) {
      if (m.ncols > 0) {
        fout.write(m(i, 0).toString)
      }
      var j = 1
      while (j < m.ncols) {
        fout.write(delim + m(i, j).toString)
        j += 1
      }
      fout.write("\n")
      i += 1
    }
    fout.close
  }

  def saveIMat(fname: String, m: IMat, compressed: Int = 0): Unit = {
    if (fname.startsWith("hdfs:")) {
      HDFSwriteMat(fname, m, compressed);
    } else if (fname.endsWith(".txt") || fname.endsWith(".txt.gz") || fname.endsWith(".txt.lz4")) {
      saveIMatTxt(fname, m, compressed)
    } else {
      val gout = getOutputStream(fname, compressed)
      saveIMat(gout, m);
      gout.close
    }
  }

  def saveIMat(gout: DataOutput, m: IMat): Unit = {
    val hints = new Array[Int](4)
    val tbuf = ByteBuffer.allocate(16).order(byteOrder)
    hints(0) = 110 // 1=dense, 1=int
    hints(1) = m.nrows
    hints(2) = m.ncols
    hints(3) = 0
    writeSomeInts(gout, hints, tbuf, 4)
    val bsize = 4 * m.ncols * m.nrows
    val buff = ByteBuffer.allocate(if (bsize > 0 && bsize < DEFAULT_BUFSIZE) bsize else DEFAULT_BUFSIZE).order(byteOrder)
    writeSomeInts(gout, m.data, buff, m.nrows * m.ncols)
  }

  def saveIMatTxt(fname: String, m: IMat, compressed: Int = 0, delim: String = "\t"): Unit = {
    val gout = getOutputStream(fname, compressed)
    val fout = new BufferedWriter(new OutputStreamWriter(gout.asInstanceOf[DataOutputStream]))
    var i = 0
    while (i < m.nrows) {
      if (m.ncols > 0) {
        fout.write(m(i, 0).toString)
      }
      var j = 1
      while (j < m.ncols) {
        fout.write(delim + m(i, j).toString)
        j += 1
      }
      fout.write("\n")
      i += 1
    }
    fout.close
  }

  def saveLMat(fname: String, m: LMat, compressed: Int = 0): Unit = {
    if (fname.startsWith("hdfs:")) {
      HDFSwriteMat(fname, m, compressed);
    } else if (fname.endsWith(".txt") || fname.endsWith(".txt.gz") || fname.endsWith(".txt.lz4")) {
      saveLMatTxt(fname, m, compressed)
    } else {
      val gout = getOutputStream(fname, compressed)
      saveLMat(gout, m)
      gout.close
    }
  }

  def saveLMat(gout: DataOutput, m: LMat): Unit = {
    val hints = new Array[Int](4)
    val tbuf = ByteBuffer.allocate(16).order(byteOrder)
    hints(0) = 120 // 1=dense, 2=long
    hints(1) = m.nrows
    hints(2) = m.ncols
    hints(3) = 0
    writeSomeInts(gout, hints, tbuf, 4)
    val bsize = 4 * m.ncols * m.nrows
    val buff = ByteBuffer.allocate(if (bsize > 0 && bsize < DEFAULT_BUFSIZE) bsize else DEFAULT_BUFSIZE).order(byteOrder)
    writeSomeLongs(gout, m.data, buff, m.nrows * m.ncols)
  }

  def saveLMatTxt(fname: String, m: LMat, compressed: Int = 0, delim: String = "\t"): Unit = {
    val gout = getOutputStream(fname, compressed)
    val fout = new BufferedWriter(new OutputStreamWriter(gout.asInstanceOf[DataOutputStream]))
    var i = 0
    while (i < m.nrows) {
      if (m.ncols > 0) {
        fout.write(m(i, 0).toString)
      }
      var j = 1
      while (j < m.ncols) {
        fout.write(delim + m(i, j).toString)
        j += 1
      }
      fout.write("\n")
      i += 1
    }
    fout.close
  }

  def saveDMat(fname: String, m: DMat, compressed: Int = 0): Unit = {
    if (fname.startsWith("hdfs:")) {
      HDFSwriteMat(fname, m, compressed);
    } else if (fname.endsWith(".txt") || fname.endsWith(".txt.gz") || fname.endsWith(".txt.lz4")) {
      saveDMatTxt(fname, m, compressed)
    } else {
      val gout = getOutputStream(fname, compressed)
      saveDMat(gout, m)
      gout.close
    }
  }

  def saveDMat(gout: DataOutput, m: DMat): Unit = {
    val hints = new Array[Int](4)
    val tbuf = ByteBuffer.allocate(16).order(byteOrder)
    hints(0) = 140 // 1=dense, 4=double
    hints(1) = m.nrows
    hints(2) = m.ncols
    hints(3) = 0
    writeSomeInts(gout, hints, tbuf, 4)
    val bsize = 8 * m.ncols * m.nrows
    val buff = ByteBuffer.allocate(if (bsize > 0 && bsize < DEFAULT_BUFSIZE) bsize else DEFAULT_BUFSIZE).order(byteOrder)
    writeSomeDoubles(gout, m.data, buff, m.nrows * m.ncols)
  }

  def saveDMatTxt(fname: String, m: DMat, compressed: Int = 0, delim: String = "\t"): Unit = {
    val gout = getOutputStream(fname, compressed)
    val fout = new BufferedWriter(new OutputStreamWriter(gout.asInstanceOf[DataOutputStream]))
    var i = 0
    while (i < m.nrows) {
      if (m.ncols > 0) {
        fout.write(m(i, 0).toString)
      }
      var j = 1
      while (j < m.ncols) {
        fout.write(delim + m(i, j).toString)
        j += 1
      }
      fout.write("\n")
      i += 1
    }
    fout.close
  }

  def loadSMat(fname: String, compressed: Int = 0): SMat = {
    if (fname.startsWith("hdfs:")) {
      HDFSreadMat(fname, null).asInstanceOf[SMat];
    } else {
      val gin = getInputStream(fname, compressed);
      val out = loadSMat(gin, null);
      gin.close;
      out
    }
  }

  def loadSMat(gin: DataInput, omat: Mat): SMat = {
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
    val norows: Boolean = (ftype / 100 == 3)
    val out = if (norows) {
      SMat.SnoRows(nrows, ncols, nnz)
    } else {
      SMat(nrows, ncols, nnz)
    }
    readSomeInts(gin, out.jc, buff, ncols + 1)
    if (!norows) readSomeInts(gin, out.ir, buff, nnz)
    readSomeFloats(gin, out.data, buff, nnz)
    if (Mat.ioneBased == 1) {
      MatHDF5.addOne(out.jc)
      if (!norows) MatHDF5.addOne(out.ir)
    }
    out
  }

  def loadSDMat(fname: String, compressed: Int = 0): SDMat = {
    if (fname.startsWith("hdfs:")) {
      HDFSreadMat(fname, null).asInstanceOf[SDMat];
    } else {
      val gin = getInputStream(fname, compressed)
      val out = loadSDMat(gin, null);
      gin.close
      out
    }
  }

  def loadSDMat(gin: DataInput, omat: Mat): SDMat = {
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
    val norows: Boolean = (ftype / 100 == 3)
    val out = if (norows) {
      SDMat.SnoRows(nrows, ncols, nnz)
    } else {
      SDMat(nrows, ncols, nnz)
    }
    readSomeInts(gin, out.jc, buff, ncols + 1)
    if (!norows) readSomeInts(gin, out.ir, buff, nnz)
    readSomeDoubles(gin, out.data, buff, nnz)
    if (Mat.ioneBased == 1) {
      MatHDF5.addOne(out.jc)
      if (!norows) MatHDF5.addOne(out.ir)
    }
    out
  }

  def loadSBMat(fname: String, compressed: Int = 0): SBMat = {
    if (fname.startsWith("hdfs:")) {
      HDFSreadMat(fname, null).asInstanceOf[SBMat];
    } else {
      val gin = getInputStream(fname, compressed)
      val out = loadSBMat(gin, null);
      gin.close
      out
    }
  }

  def loadSBMat(gin: DataInput, omat: Mat): SBMat = {
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
    val norows: Boolean = (ftype / 100 == 3)
    val out = if (norows) {
      SBMat.SnoRows(nrows, ncols, nnz)
    } else {
      SBMat(nrows, ncols, nnz)
    }
    readSomeInts(gin, out.jc, buff, ncols + 1)
    if (!norows) readSomeInts(gin, out.ir, buff, nnz)
    readSomeBytes(gin, out.data, nnz)
    if (Mat.ioneBased == 1) {
      MatHDF5.addOne(out.jc)
      if (!norows) MatHDF5.addOne(out.ir)
    }
    out
  }

  def loadCSMatTxt(fname: String, compressed: Int = 0): CSMat = {
    val gin = getInputStream(fname, compressed);
    val buf = new ArrayBuffer[String];
    val istream = new BufferedReader(new InputStreamReader(gin));
    var r: String = null;
    while ({ r = istream.readLine(); r != null }) {
      buf.append(r);
    }
    val cs = new CSMat(buf.length, 1, buf.toArray);
    cs
  }

  def loadCSMat(fname: String, compressed: Int = 0): CSMat = {
    if (fname.startsWith("hdfs:")) {
      HDFSreadMat(fname, null).asInstanceOf[CSMat];
    } else if (fname.endsWith(".txt") || fname.endsWith(".txt.gz") || fname.endsWith(".txt.lz4")) {
      loadCSMatTxt(fname, compressed);
    } else {
      val gin = getInputStream(fname, compressed);
      val out = loadCSMat(gin, null);
      gin.close
      out
    }
  }

  def loadCSMat(gin: DataInput, omat: Mat): CSMat = {
    val buff = ByteBuffer.allocate(DEFAULT_BUFSIZE).order(byteOrder);
    val hints = new Array[Int](4);
    readSomeInts(gin, hints, buff, 4);
    val ftype = hints(0);
    val nrows = hints(1);
    val ncols = hints(2);
    val nnz = hints(3);
    val len = nrows * ncols;
    if (ftype != 202 && ftype != 302) {
      throw new RuntimeException("loadCSMat expected type field 202 or 302 but was %d" format ftype);
    }
    val out = CSMat(nrows, ncols);
    val jc = new Array[Long](len + 1);

    readSomeLongs(gin, jc, buff, len + 1);
    var i = 0;
    var maxlen = 0;
    while (i < len) {
      val siz = (jc(i + 1) - jc(i)).toInt;
      if (siz > maxlen) maxlen = siz;
      i += 1;
    }
    val chbuf = new Array[Byte](maxlen);
    i = 0;
    while (i < len) {
      val siz = (jc(i + 1) - jc(i)).toInt;
      readSomeBytes(gin, chbuf, siz);
      out.data(i) = new String(chbuf, 0, siz, SBMat.encoding);
      i += 1;
    }
    out
  }

  def saveSMat(fname: String, m: SMat, compressed: Int = 0): Unit = {
    if (fname.startsWith("hdfs:")) {
      HDFSwriteMat(fname, m, compressed)
    } else {
      val gout = getOutputStream(fname, compressed)
      saveSMat(gout, m)
      gout.close
    }
  }

  def saveSMat(gout: DataOutput, m: SMat) = {
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
    val buff = ByteBuffer.allocate(4 * math.min(DEFAULT_BUFSIZE / 4, math.max(m.ncols + 1, m.nnz))).order(byteOrder)
    try {
      MatHDF5.subOne(m.jc)
      writeSomeInts(gout, m.jc, buff, m.ncols + 1)
      if (m.ir != null) {
        MatHDF5.subOne(m.ir)
        writeSomeInts(gout, m.ir, buff, m.nnz)
      }
      writeSomeFloats(gout, m.data, buff, m.nnz)
    } catch {
      case e: Exception => {
        MatHDF5.addOne(m.jc)
        if (m.ir != null) MatHDF5.addOne(m.ir)
        throw new RuntimeException("Exception in saveSMat " + e)
      }
      case _: Throwable => {
        MatHDF5.addOne(m.jc)
        if (m.ir != null) MatHDF5.addOne(m.ir)
        throw new RuntimeException("Problem in saveSMat")
      }
    }
    MatHDF5.addOne(m.jc)
    if (m.ir != null) MatHDF5.addOne(m.ir)
  }

  def saveSDMat(fname: String, m: SDMat, compressed: Int = 0): Unit = {
    if (fname.startsWith("hdfs:")) {
      HDFSwriteMat(fname, m, compressed)
    } else {
      val gout = getOutputStream(fname, compressed)
      saveSDMat(gout, m);
      gout.close
    }
  }

  def saveSDMat(gout: DataOutput, m: SDMat): Unit = {
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
    val buff = ByteBuffer.allocate(4 * math.min(DEFAULT_BUFSIZE / 4, math.max(m.ncols + 1, m.nnz))).order(byteOrder)
    try {
      MatHDF5.subOne(m.jc)
      writeSomeInts(gout, m.jc, buff, m.ncols + 1)
      if (m.ir != null) {
        MatHDF5.subOne(m.ir)
        writeSomeInts(gout, m.ir, buff, m.nnz)
      }
      writeSomeDoubles(gout, m.data, buff, m.nnz)
    } catch {
      case e: Exception => {
        MatHDF5.addOne(m.jc)
        if (m.ir != null) MatHDF5.addOne(m.ir)
        throw new RuntimeException("Exception in saveSDMat " + e)
      }
      case _: Throwable => {
        MatHDF5.addOne(m.jc)
        if (m.ir != null) MatHDF5.addOne(m.ir)
        throw new RuntimeException("Problem in saveSDMat")
      }
    }
    MatHDF5.addOne(m.jc)
    if (m.ir != null) MatHDF5.addOne(m.ir)
  }

  def saveSBMat(fname: String, m: SBMat, compressed: Int = 0): Unit = {
    if (fname.startsWith("hdfs:")) {
      HDFSwriteMat(fname, m, compressed)
    } else {
      val gout = getOutputStream(fname, compressed)
      saveSBMat(gout, m)
      gout.close
    }
  }

  def saveSBMat(gout: DataOutput, m: SBMat): Unit = {
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
    val buff = ByteBuffer.allocate(4 * math.min(DEFAULT_BUFSIZE / 4, math.max(m.ncols + 1, m.nnz))).order(byteOrder)
    try {
      MatHDF5.subOne(m.jc)
      writeSomeInts(gout, m.jc, buff, m.ncols + 1)
      if (m.ir != null) {
        MatHDF5.subOne(m.ir)
        writeSomeInts(gout, m.ir, buff, m.nnz)
      }
      gout.write(m.data, 0, m.nnz)
    } catch {
      case e: Exception => {
        MatHDF5.addOne(m.jc)
        if (m.ir != null) MatHDF5.addOne(m.ir)
        throw new RuntimeException("Exception in saveSMat " + e)
      }
      case _: Throwable => {
        MatHDF5.addOne(m.jc)
        if (m.ir != null) MatHDF5.addOne(m.ir)
        throw new RuntimeException("Problem in saveSMat")
      }
    }
    MatHDF5.addOne(m.jc)
    if (m.ir != null) MatHDF5.addOne(m.ir)

  }

  def saveCSMat(fname: String, m: CSMat, compressed: Int = 0): Unit = {
    if (fname.startsWith("hdfs:")) {
      HDFSwriteMat(fname, m, compressed)
    } else {
      val gout = getOutputStream(fname, compressed);
      saveCSMat(gout, m)
      gout.close
    }
  }

  def saveCSMat(gout: DataOutput, m: CSMat): Unit = {
    val jc = new Array[Long](m.length + 1);
    jc(0) = 0;
    var i = 0;
    var nnz = 0L;
    while (i < m.length) {
      nnz += m.data(i).getBytes(SBMat.encoding).length;
      jc(i + 1) = nnz;
      i += 1;
    }
    val hints = new Array[Int](4)
    val tbuf = ByteBuffer.allocate(16).order(byteOrder)
    hints(0) = 302 // 3=sparse:norows, 0=byte, 2=long
    hints(1) = m.nrows
    hints(2) = m.ncols
    hints(3) = math.min(0x7fffffffL, nnz).toInt
    writeSomeInts(gout, hints, tbuf, 4)
    val buff = ByteBuffer.allocate(4 * math.min(DEFAULT_BUFSIZE / 4, math.max(m.ncols + 1, nnz)).toInt).order(byteOrder);
    try {
      writeSomeLongs(gout, jc, buff, m.length + 1);
      i = 0;
      while (i < m.length) {
        val bytes = m.data(i).getBytes(SBMat.encoding);
        gout.write(bytes, 0, bytes.length);
        i += 1;
      }
    } catch {
      case e: Exception => {
        throw new RuntimeException("Exception in saveCSMat " + e)
      }
      case _: Throwable => {
        throw new RuntimeException("Problem in saveCSMat")
      }
    }
  }

  /** Load a file in ASCII LibSVM format */
  /* Outputs a data matrix first, and then a matrix c of cat labels where c(i) = label of instance i. */
  /* data is an SMat, labels are integer, weights are FMat. */
  def loadLibSVM(fname: String, nrows: Int, compressed: Int = 0, oneBased: Int = 0): (SMat, IMat, FMat) = {
    var fin = new BufferedReader(new InputStreamReader(getInputStream(fname, compressed).asInstanceOf[DataInputStream]));
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
    fin = new BufferedReader(new InputStreamReader(getInputStream(fname, compressed).asInstanceOf[DataInputStream]));
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
        datamat.ir(nnz) = indx - oneBased + ioneBased;
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

  def saveLibSVM(fname: String, data: SMat, labels: SMat): Unit = saveLibSVM(fname, data, labels, null, 0);

  def saveLibSVM(fname: String, data: SMat, labels: SMat, weights: FMat): Unit = saveLibSVM(fname, data, labels, weights, 0);

  def saveLibSVM(fname: String, data: SMat, labels: SMat, weights: FMat, compressed: Int, oneBased: Int = 0): Unit = {
    val fout = new BufferedWriter(new OutputStreamWriter(getOutputStream(fname, compressed).asInstanceOf[DataOutputStream]));
    val (ilab, jlab) = find2(labels);
    val jc = data.jc;
    val ioffset = Mat.ioneBased;
    var i = 0;
    while (i < ilab.length) {
      fout.write(ilab.data(i).toString);
      val icol = jlab.data(i);
      if (weights.asInstanceOf[AnyRef] != null) {
        fout.write(":" + weights.data(icol))
      }
      var j = jc(icol) - ioffset
      while (j < jc(icol + 1) - ioffset) {
        fout.write(" " + (data.ir(j) - ioffset + oneBased).toString + ":" + data.data(j).toString);
        j += 1
      }
      fout.write("\n")
      i += 1
    }
    fout.close()
  }

  def saveVW(fname: String, sdata: SMat, labels: SMat): Unit = saveVW(fname, sdata, labels, 0);

  def saveVW(fname: String, sdata: SMat, labels: SMat, compressed: Int, oneBased: Int = 0): Unit = {
    val fout = new BufferedWriter(new OutputStreamWriter(getOutputStream(fname, compressed).asInstanceOf[DataOutputStream]));
    val jc = sdata.jc;
    val ir = sdata.ir;
    val data = sdata.data;
    val jci = labels.jc;
    val iri = labels.ir;
    val ioffset = Mat.ioneBased;
    var i = 0;
    while (i < sdata.ncols) {
      var j = jci(i);
      while (j < jci(i + 1)) {
        if (j > jci(i)) fout.write(",");
        fout.write((iri(j - ioffset) - ioffset + oneBased).toString);
        j += 1;
      }
      fout.write(" |");
      j = jc(i);
      while (j < jc(i + 1)) {
        fout.write(" %d:%f" format (ir(j - ioffset) - ioffset + oneBased, data(j - ioffset)));
        j += 1;
      }
      fout.write("\n")
      i += 1
    }
    fout.close()
  }

  def testLoad(fname: String, varname: String, n: Int) = {
    val a = new Array[SMat](n)
    var ndone = izeros(n, 1)
    for (i <- 0 until n) {
      Future {
        a(i) = loadSMat(("/disk%02d/" format i) + fname)
        ndone(i) = 1
      }
    }
    while (SciFunctions.sum(ndone).v < n) { Thread.sleep(10) }
    a
  }

  def apply(dirname: String, filepat: String, varname: String, catd: Int) {
    var files: ListBuffer[String] = new ListBuffer[String]
    val dir: File = new File(dirname)
    val slen = dir.getName.length + 1

    def searchDir(dir: File) {
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
    val blkinds = new Array[Int](files.length + 1)
    var i = 0
    var nrows = -1
    var ncols = -1
    files.foreach((fn: String) => {
      val (nr, nc) = MatHDF5.readMat2Dims(fn, varname)
      if (catd == 2) {
        if (nrows >= 0) {
          if (nr != nrows) {
            throw new RuntimeException("incorrect number of rows in file " + fn)
          }
        } else {
          nrows = nr.asInstanceOf[Int]
        }
        blkinds(i + 1) = blkinds(i) + nc.asInstanceOf[Int]
        i += 1
      } else if (catd == 1) {
        if (ncols >= 0) {
          if (nc != 1) {
            throw new RuntimeException("incorrect number of cols in file " + fn)
          }
        } else {
          ncols = 1
        }
        blkinds(i + 1) = blkinds(i) + nr.asInstanceOf[Int]
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

