/*                     __                                               *\
**     ________ ___   / /  ___     Scala API                            **
**    / __/ __// _ | / /  / _ |    (c) 2006-2009, Ross Judson           **
**  __\ \/ /__/ __ |/ /__/ __ |    http://scala-lang.org/               **
** /____/\___/_/ |_/____/_/ | |                                         **
**                          |/                                          **
\*                                                                      */

package BIDMat

import scala.reflect.ClassManifest
import scala.math.Ordering

object Sorting {
  
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
}
