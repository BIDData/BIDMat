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
  
  def readSomeInts(din:DataInputStream, a:Array[Int], buf:Array[Byte], n:Int) {
    var nread = 0
    while (nread < 4*n) {
      val readnow = din.read(buf, 0, math.min(buf.length, 4*n-nread))
      memcpybi(readnow, buf, 0, a, nread)
      nread += readnow
    }
  }
  
  def readSomeFloats(din:DataInputStream, a:Array[Float], buf:Array[Byte], n:Int) {
    var nread = 0
    while (nread < 4*n) {
      val readnow = din.read(buf, 0, math.min(buf.length, 4*n-nread))
      memcpybf(readnow, buf, 0, a, nread)
      nread += readnow
    }
  }
  
  def saveSMat(fname:String, m:SMat):Unit = {
    val fout = new FileOutputStream(fname)
    val bout = new BufferedOutputStream(fout, 1024*1024)
    val gout = new GZIPOutputStream(bout)
    val dout = new DataOutputStream(gout)
    dout.writeInt(1)
    dout.writeInt(m.nrows)
    dout.writeInt(m.ncols)
    dout.writeInt(m.nnz)
    val buff = new Array[Byte](4*math.max(m.ncols+1, m.nnz))
    try {
    	MatHDF5.subOne(m.jc)
    	MatHDF5.subOne(m.ir)
    	memcpyib(4(m.ncols+1), m.jc, 0, buff, 0)
    	dout.write(buff, 0, 4*(m.ncols+1))
    	dout.flush
    	memcpyib(4*m.nnz, m.ir, 0, buff, 0)
    	dout.write(buff, 0, 4*m.nnz)
    	dout.flush
    	memcpyfb(4*m.nnz, m.data, 0, buff, 0)
    	dout.write(buff, 0, 4*m.nnz)
    	dout.flush
    } catch {
      case _ => {
      	MatHDF5.addOne(m.jc)
      	MatHDF5.addOne(m.ir)
      	throw new RuntimeException("Problem in saveSMat")
      }
    }
    MatHDF5.addOne(m.jc)
    MatHDF5.addOne(m.ir)
    dout.close
  }
  
  def loadSMat(fname:String):SMat = {
    val fin = new FileInputStream(fname)
    val bin = new BufferedInputStream(fin, 1024*1024)
    val gin = new GZIPInputStream(bin)
    val din = new DataInputStream(gin)
    val ftype = din.readInt
    val nrows = din.readInt
    val ncols = din.readInt
    val nnz = din.readInt
    val out = SMat(nrows, ncols, nnz)
    val buff = new Array[Byte](1024*1024)
    readSomeInts(din, out.jc, buff, ncols+1)
    readSomeInts(din, out.ir, buff, nnz)
    readSomeFloats(din, out.data, buff, nnz)
    MatHDF5.addOne(out.jc)
    MatHDF5.addOne(out.ir)
    din.close
    out
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
    while (SciFunctions.sum(ndone).v < n) {Thread.sleep(100)}
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






