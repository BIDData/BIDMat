package BIDMat

import scala.reflect.ClassManifest
import scala.math.Ordering

object Sorting {
  
  def quickSort2[@specialized(Double,Float,Int,Byte) T](a:Array[T], ii:Array[Int], lo:Int, hi:Int, stride:Int)
  (implicit ord:Ordering[T]):Unit = {
    if ((hi - lo)/stride > 0) {
    	if ((hi - lo)/stride <= 10) {
    		isort(a, ii, lo, hi, stride)
    	} else {
    		val ip = partition(a, ii, lo, hi, stride)
//    		println("part %d, %f" format ((hi - lo) , 1.0f * (ip - lo) / (hi - lo)))
    		quickSort2(a, ii, lo, ip, stride)
    		quickSort2(a, ii, ip, hi, stride)
    	}
    }
  }
  
  def isort[@specialized(Double,Float,Int,Byte) T](a:Array[T], ii:Array[Int], lo:Int, hi:Int, stride:Int)
  (implicit ord:Ordering[T]):Unit = {
    var i = lo
    while (i != hi) {
      var j = i+stride
      var imin = i
      var vmin = a(i)
      while (j != hi) {
        if (ord.lteq(a(j), vmin) && (ord.lt(a(j), vmin) || ii(j) < ii(imin))) {
          vmin = a(j)
          imin = j
        }
        j += stride
      }
      a(imin) = a(i)
      a(i) = vmin
      val itmp = ii(imin)
      ii(imin) = ii(i)
      ii(i) = itmp
      i += stride
    }
  }
  
  def med3[@specialized(Double,Float,Int,Byte) T](a:Array[T], ii:Array[Int], lo:Int, hi:Int, stride:Int)
  (implicit ord:Ordering[T]):Int = {
    val nv = (hi - lo)/stride
    val i1 = lo + stride*(math.floor(nv*java.lang.Math.random()).asInstanceOf[Int])
    val i2 = lo + stride*(math.floor(nv*java.lang.Math.random()).asInstanceOf[Int])
    val i3 = lo + stride*(math.floor(nv*java.lang.Math.random()).asInstanceOf[Int])
    val v1 = a(i1)
    val v2 = a(i2)
    val v3 = a(i3)
    val ii1 = ii(i1)
    val ii2 = ii(i2)
    val ii3 = ii(i3)
    if (ord.gteq(v2,v1) && (ord.gt(v2,v1) || ii2 > ii1)) {
    	if (ord.gteq(v3, v2) && (ord.gt(v3,v2) || ii3 > ii2)) i2 else {
    		if (ord.gteq(v3, v1) && (ord.gt(v3,v1) || ii3 > ii1)) i3 else i1
    	}
    } else {
    	if (ord.gteq(v3, v1) && (ord.gt(v3,v1) || ii3 > ii1)) i1 else {
    		if (ord.gteq(v3, v2) && (ord.gt(v3,v2) || ii3 > ii2)) i3 else i2
    	}
    }
  }
  
  def med9[@specialized(Double,Float,Int,Byte) T](a:Array[T], ii:Array[Int], lo:Int, hi:Int, stride:Int)
  (implicit ord:Ordering[T]):Int = {
    val i1 = med3(a, ii, lo, hi, stride)
    val i2 = med3(a, ii, lo, hi, stride)
    val i3 = med3(a, ii, lo, hi, stride)
    val v1 = a(i1)
    val v2 = a(i2)
    val v3 = a(i3)
    val ii1 = ii(i1)
    val ii2 = ii(i2)
    val ii3 = ii(i3)
    if (ord.gteq(v2,v1) && (ord.gt(v2,v1) || ii2 > ii1)) {
    	if (ord.gteq(v3, v2) && (ord.gt(v3,v2) || ii3 > ii2)) i2 else {
    		if (ord.gteq(v3, v1) && (ord.gt(v3,v1) || ii3 > ii1)) i3 else i1
    	}
    } else {
    	if (ord.gteq(v3, v1) && (ord.gt(v3,v1) || ii3 > ii1)) i1 else {
    		if (ord.gteq(v3, v2) && (ord.gt(v3,v2) || ii3 > ii2)) i3 else i2
    	}
    }
  }
  
   def partition[@specialized(Double,Float,Int,Byte) T](a:Array[T], ii:Array[Int], lo:Int, hi:Int, stride:Int)
   (implicit ord:Ordering[T]):Int = {
      val sstride = math.signum(stride)
  		val im = if ((hi - lo)/stride > 100) {
  			med9(a, ii, lo, hi, stride)	  
  		} else {
  		  med3(a, ii, lo, hi, stride)
  		}
  		var v = a(im)
  		var iv = ii(im)
//  		println("med %d, %d, %d" format (im, lo, hi))
  		var done = false
  		var i = lo - stride
  		var j = hi 
  		while (! done) { 
  			i += stride
  			j -= stride
  			while ((hi-i)*sstride > sstride*stride && (ord.lteq(a(i), v) && (ord.lt(a(i), v) || ii(i) <= iv))) {i += stride}
  			while ((j-lo)*sstride > 0              && (ord.gteq(a(j), v) && (ord.gt(a(j), v) || ii(j) > iv))) {j -= stride}
  			if ((i - j)*sstride >= 0) {
  				done = true
  			} else {
  				val atmp = a(i)
  				a(i) = a(j)
  				a(j) = atmp
  				val itmp = ii(i)
  				ii(i) = ii(j)
  				ii(j) = itmp
  			}
  		}
  		j + stride
   }
  
  def quickSort[@specialized(Double, Float, Int, Byte) T](a:Array[T])(implicit ord:Ordering[T]) = { 
    def comp(i:Int, j:Int):Int = {
      ord.compare(a(i),a(j))
    }
    def swap(i:Int, j:Int) = {
      val tmp = a(i)
      a(i) = a(j)
      a(j) = tmp
    }
    sort1(comp, swap, 0, a.length) 
  }
  
  def quickSort(comp:(Int,Int)=>Int, swap: (Int,Int) => Unit, start:Int, len:Int) { sort1(comp, swap, start, len) }
  
  private def sort1(comp: (Int, Int) => Int, swap: (Int,Int) => Unit, off: Int, len: Int) {
    
    def vecswap(_a: Int, _b: Int, n: Int) {
      var a = _a
      var b = _b
      var i = 0
      while (i < n) {
        swap(a, b)
        i += 1
        a += 1
        b += 1
      }
    }
    def med3(a: Int, b: Int, c: Int) = {
      if (comp(a,b) < 0) {
        if (comp(b,c) < 0) b else if (comp(a,c) < 0) c else a
      } else {
        if (comp(b,c) > 0) b else if (comp(a,c) > 0) c else a
      }
    }
    def sort2(off: Int, len: Int) {
      // Insertion sort on smallest arrays
      if (len < 7) {
        var i = off
        while (i < len + off) {
          var j = i
          while (j > off && comp(j-1,j) > 0) {
            swap(j, j-1)
            j -= 1
          }
          i += 1
        }
      } else {
        // Choose a partition element, v
        var m = off + (len >> 1) // Small arrays, middle element
        if (len > 30) {
          var l = off
          var n = off + len - 1
          if (len > 300) { // Big arrays, pseudomedian of 9
            val s = len / 8
            l = med3(l, l+s, l+2*s)
            m = med3(m-s, m, m+s)
            n = med3(n-2*s, n-s, n)
          }
          m = med3(l, m, n) // Mid-size, med of 3
        }

        // Establish Invariant: v* (<v)* (>v)* v*
        var a = off
        var b = a
        var c = off + len - 1
        var d = c
        var done = false
        while (!done) {
          var pp = -1
          while (b <= c && pp <= 0) {
            pp = comp(b, m) 
            if (pp == 0) {
              swap(a, b)
              m = a
              a += 1
            }
            if (pp <= 0) b += 1
          }
          pp = 1
          while (c >= b &&  pp >= 0) {
            pp = comp(c, m)
            if (pp == 0) {
              swap(c, d)
              m = d
              d -= 1
            }
            if (pp >= 0) c -= 1
          }
          if (b > c) {
            done = true
          } else {
            swap(b, c)
            c -= 1
            b += 1
          }
        }

        // Swap partition elements back to middle
        val n = off + len
        var s = math.min(a-off, b-a)
        vecswap(off, b-s, s)
        s = math.min(d-c, n-d-1)
        vecswap(b, n-s, s)

        // Recursively sort non-partition-elements
        s = b - a
        if (s > 1)
          sort2(off, s)
        s = d - c
        if (s > 1)
          sort2(n-s, s)
      }
    }
    sort2(off, len)
  }
  
  def main(args:Array[String]) = {
    import BIDMat.SciFunctions._
    import BIDMat.MatFunctions._
    val n = args(0).toInt
    val a = SciFunctions.rand(n, 1)
    val ii = MatFunctions.icol(0->n)
    quickSort2(a.data, ii.data, 0, n, 1)
    println("check %d" format find(a(1->n,0) < a(0->(n-1),0)).length)
  }
}
