/* Copyright (c) 2012, Regents of the University of California                     */
/* All rights reserved.                                                            */

/* Redistribution and use in source and binary forms, with or without              */
/* modification, are permitted provided that the following conditions are met:     */
/*     * Redistributions of source code must retain the above copyright            */
/*       notice, this list of conditions and the following disclaimer.             */
/*     * Redistributions in binary form must reproduce the above copyright         */
/*       notice, this list of conditions and the following disclaimer in the       */
/*       documentation and/or other materials provided with the distribution.      */
/*     * Neither the name of the <organization> nor the                            */
/*       names of its contributors may be used to endorse or promote products      */
/*       derived from this software without specific prior written permission.     */

/* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND */
/* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED   */
/* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE          */
/* DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY              */
/* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES      */
/* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;    */
/* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND     */
/* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT      */
/* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS   */
/* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                    */

package BIDMat

import java.io._
import scala.util.matching.Regex
import Regex._
import scala.collection.mutable._

case class HMat(nr:Int, nc:Int, fileList:List[String], varname:String, blkinds:Array[Int], catdim:Int) extends Mat(nr, nc) {

  var fnameCache:String = null

  var fmatCache:FMat = null
  
// Implement slicing from a hard disk matrix
  override def apply(a:IMat, b:IMat):FMat = { 
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
    var out:FMat = null
    while (i <= b.length) {
    	if (i < b.length) locs(i) = findindx(b(i))
    	if (i == b.length || locs(i) != locs(iblk)) {
    		if (fnameCache == null || fileList(locs(iblk)) != fnameCache) {
    			fmatCache = MatHDF5.hload(fileList(locs(iblk)), varname).asInstanceOf[FMat] 
    			fnameCache = fileList(locs(iblk))
    		}
      	val newmat = fmatCache(a, b(MatFunctions.irow(iblk->i)))
      	if (out != null) {
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






