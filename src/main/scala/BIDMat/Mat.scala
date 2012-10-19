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

abstract class Mat(nr:Int, nc:Int) {
  val nrows = nr
  val ncols = nc

  def length = nr*nc

  def notImplemented0(s:String):Mat = { 
    throw new RuntimeException("operator "+s+" not implemented for "+this)
  }
  
  def notImplemented1(s:String,that:Mat):Mat = { 
    throw new RuntimeException("operator "+s+" not implemented for "+this+" and "+that)
  }
  
  def t = notImplemented0("t")  
  
  def copy(a:Mat) = notImplemented0("copy");
  def copy = notImplemented0("copy");
  def zeros(nr:Int, nc:Int) = notImplemented0("zeros");
  def ones(nr:Int, nc:Int) = notImplemented0("ones");
  def clearUpper(i:Int) = notImplemented0("clearUpper");
  def clearLower(i:Int) = notImplemented0("clearLower"); 
  def clearUpper = notImplemented0("clearUpper");
  def clearLower = notImplemented0("clearLower");
  
  def apply(a:IMat):Mat = notImplemented0("linear array access");
  def apply(a:IMat, b:IMat):Mat = notImplemented0("block array access");
  def apply(a:IMat, b:Int):Mat = notImplemented0("block array access");	
  def apply(a:Int, b:IMat):Mat = notImplemented0("block array access");
  

  def + (b : Mat):Mat = notImplemented1("+", b)
  def - (b : Mat):Mat = notImplemented1("-", b)
  def * (b : Mat):Mat = notImplemented1("*", b)
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
  def ~ (b : Mat):Pair = throw new RuntimeException("operator ~ not implemented for "+this+" and "+b)
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
}

object Mat {
  import Ordered._
  import scala.tools.jline.TerminalFactory
  
  var compress:Boolean = true
  
  var noMKL:Boolean = false
  
  var nflops = 0L
  
  var oneBased = 0
  
  var ioneBased = 1
  
  var hasCUDA = 0;
  
  val terminalWidth = TerminalFactory.create.getWidth

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
    while (i < a.length) {
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
