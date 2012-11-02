package BIDMat

import java.io._
import java.util.zip._
import scala.util.matching.Regex
import Regex._
import scala.collection.mutable._
import scala.actors._
import scala.actors.Actor._
import MatFunctions._

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
  
  def saveSMat(fname:String, m:SMat):Unit = {
    val fout = new FileOutputStream(fname)
    val gout = new GZIPOutputStream(fout)
    val dout = new DataOutputStream(gout)
    dout.writeInt(1)
    dout.writeInt(m.nrows)
    dout.writeInt(m.ncols)
    dout.writeInt(m.nnz)
    var i = 0; while(i <= m.ncols) { dout.writeInt(m.jc(i)); i+=1 }
    i = 0; while(i < m.nnz) { dout.writeInt(m.ir(i)); i+=1 }
    i = 0; while(i < m.nnz) { dout.writeFloat(m.data(i)); i+=1 }
    dout.close
    gout.close
    fout.close
  }
  
  def loadSMat(fname:String):SMat = {
    val fin = new FileInputStream(fname)
    val gin = new GZIPInputStream(fin)
    val din = new DataInputStream(gin)
    val ftype = din.readInt
    val nrows = din.readInt
    val ncols = din.readInt
    val nnz = din.readInt
    val buff = new Array[Byte](4*math.max(ncols+1, nnz))
    din.read(buff, 0, 4*(ncols+1))
    din.read(buff, 0, 4*nnz)
    din.read(buff, 0, 4*nnz)
    val out = SMat(nrows, ncols, nnz)
/*    var i = 0; while(i <= ncols) { out.jc(i) = din.readInt; i+=1 }
    i = 0; while(i < nnz) { out.ir(i) = din.readInt; i+=1 }
    i = 0; while(i < nnz) { out.data(i) = din.readFloat; i+=1 } */
    din.close
    gin.close
    fin.close
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






