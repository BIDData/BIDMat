package BIDMat

import scala.reflect.ClassManifest
import scala.math.Ordering
import scala.concurrent.future
import scala.concurrent.ExecutionContext.Implicits.global

object Sorting {

	def quickSort2[T](ga:Array[T], ii:Array[Int], lo:Int, hi:Int, stride:Int):Unit = {
			ga match {
			case a:Array[Float] => quickSort2(a, ii, lo, hi, stride, Mat.numThreads/2)
			case a:Array[Double] => quickSort2(a, ii, lo, hi, stride, Mat.numThreads/2)
			case a:Array[Int] => quickSort2(a, ii, lo, hi, stride, Mat.numThreads/2)
			case a:Array[Long] => quickSort2(a, ii, lo, hi, stride, Mat.numThreads/2)
			}
	}

	def quickSort2(a:Array[Float], ii:Array[Int], lo:Int, hi:Int, stride:Int, nthreads:Int):Unit = {
			if ((hi - lo)/stride > 0) {
				if ((hi - lo)/stride <= 16) {
					isort(a, ii, lo, hi, stride)
				} else {
					val ip = partition(a, ii, lo, hi, stride)
					if (nthreads > 1 && (hi-lo)/stride > 400) {
						var done0 = false
						var done1 = false
						future { quickSort2(a, ii, lo, ip, stride, nthreads/2); done0 = true }
						future { quickSort2(a, ii, ip, hi, stride, nthreads/2); done1 = true }
						while (!done0 || !done1) {Thread.`yield`}
					} else {
						quickSort2(a, ii, lo, ip, stride, nthreads/2)
						quickSort2(a, ii, ip, hi, stride, nthreads/2)
					}
				}
			}
	}

	def isort(a:Array[Float], ii:Array[Int], lo:Int, hi:Int, stride:Int):Unit = {
			var i = lo
			while (i != hi) {
				var j = i+stride
				var imin = i
				var vmin = a(i)
				while (j != hi) {
					if (a(j) <= vmin && ((a(j) < vmin) || ii(j) < ii(imin))) {
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

	def med3(a:Array[Float], ii:Array[Int], lo:Int, hi:Int, stride:Int):Int = {
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
			if ((v2 >= v1) && ((v2 > v1) || ii2 > ii1)) {
				if ((v3 >= v2) && ((v3 > v2) || ii3 > ii2)) i2 else {
					if ((v3 >= v1) && ((v3 > v1) || ii3 > ii1)) i3 else i1
				}
			} else {
				if ((v3 >= v1) && ((v3 > v1) || ii3 > ii1)) i1 else {
					if ((v3 >= v2) && ((v3 > v2) || ii3 > ii2)) i3 else i2
				}
			}
	}

	def med9(a:Array[Float], ii:Array[Int], lo:Int, hi:Int, stride:Int):Int = {
			val i1 = med3(a, ii, lo, hi, stride)
			val i2 = med3(a, ii, lo, hi, stride)
			val i3 = med3(a, ii, lo, hi, stride)
			val v1 = a(i1)
			val v2 = a(i2)
			val v3 = a(i3)
			val ii1 = ii(i1)
			val ii2 = ii(i2)
			val ii3 = ii(i3)
			if ((v2 >= v1) && ((v2 > v1) || ii2 > ii1)) {
				if ((v3 >= v2) && ((v3 > v2) || ii3 > ii2)) i2 else {
					if ((v3 >= v1) && ((v3 > v1) || ii3 > ii1)) i3 else i1
				}
			} else {
				if ((v3 >= v1) && ((v3 > v1) || ii3 > ii1)) i1 else {
					if ((v3 >= v2) && ((v3 > v2) || ii3 > ii2)) i3 else i2
				}
			}
	}

	def partition(a:Array[Float], ii:Array[Int], lo:Int, hi:Int, stride:Int):Int = {
			val sstride = math.signum(stride)
			val nvals = (hi - lo)/stride
			val im = if (nvals > 600) {
				med9(a, ii, lo, hi, stride)	  
			} else if (nvals > 100) {
				med3(a, ii, lo, hi, stride)
			} else {
				lo + stride*(math.floor(nvals*java.lang.Math.random()).asInstanceOf[Int])
			}
			var v = a(im)
			var iv = ii(im)
			var done = false
			var i = lo - stride
			var j = hi 
			while (! done) { 
				i += stride
				j -= stride
				while ((hi-i)*sstride > sstride*stride && ((a(i) <= v) && ((a(i) < v) || ii(i) <= iv))) {i += stride}
				while (                                   ((a(j) >= v) && ((a(j) > v) || ii(j) > iv)))  {j -= stride}
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

	def quickSort2(a:Array[Long], ii:Array[Int], lo:Int, hi:Int, stride:Int, nthreads:Int):Unit = {
			if ((hi - lo)/stride > 0) {
				if ((hi - lo)/stride <= 16) {
					isort(a, ii, lo, hi, stride)
				} else {
					val ip = partition(a, ii, lo, hi, stride)
					if (nthreads > 1 && (hi-lo)/stride > 400) {
						var done0 = false
						var done1 = false
						future { quickSort2(a, ii, lo, ip, stride, nthreads/2); done0 = true }
						future { quickSort2(a, ii, ip, hi, stride, nthreads/2); done1 = true }
						while (!done0 || !done1) {Thread.`yield`}
					} else {
						quickSort2(a, ii, lo, ip, stride, nthreads/2)
						quickSort2(a, ii, ip, hi, stride, nthreads/2)
					}
				}
			}
	}

	def isort(a:Array[Long], ii:Array[Int], lo:Int, hi:Int, stride:Int):Unit = {
			var i = lo
			while (i != hi) {
				var j = i+stride
				var imin = i
				var vmin = a(i)
				while (j != hi) {
					if (a(j) <= vmin && ((a(j) < vmin) || ii(j) < ii(imin))) {
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

	def med3(a:Array[Long], ii:Array[Int], lo:Int, hi:Int, stride:Int):Int = {
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
			if ((v2 >= v1) && ((v2 > v1) || ii2 > ii1)) {
				if ((v3 >= v2) && ((v3 > v2) || ii3 > ii2)) i2 else {
					if ((v3 >= v1) && ((v3 > v1) || ii3 > ii1)) i3 else i1
				}
			} else {
				if ((v3 >= v1) && ((v3 > v1) || ii3 > ii1)) i1 else {
					if ((v3 >= v2) && ((v3 > v2) || ii3 > ii2)) i3 else i2
				}
			}
	}

	def med9(a:Array[Long], ii:Array[Int], lo:Int, hi:Int, stride:Int):Int = {
			val i1 = med3(a, ii, lo, hi, stride)
			val i2 = med3(a, ii, lo, hi, stride)
			val i3 = med3(a, ii, lo, hi, stride)
			val v1 = a(i1)
			val v2 = a(i2)
			val v3 = a(i3)
			val ii1 = ii(i1)
			val ii2 = ii(i2)
			val ii3 = ii(i3)
			if ((v2 >= v1) && ((v2 > v1) || ii2 > ii1)) {
				if ((v3 >= v2) && ((v3 > v2) || ii3 > ii2)) i2 else {
					if ((v3 >= v1) && ((v3 > v1) || ii3 > ii1)) i3 else i1
				}
			} else {
				if ((v3 >= v1) && ((v3 > v1) || ii3 > ii1)) i1 else {
					if ((v3 >= v2) && ((v3 > v2) || ii3 > ii2)) i3 else i2
				}
			}
	}

	def partition(a:Array[Long], ii:Array[Int], lo:Int, hi:Int, stride:Int):Int = {
			val sstride = math.signum(stride)
			val nvals = (hi - lo)/stride
			val im = if (nvals > 600) {
				med9(a, ii, lo, hi, stride)	  
			} else if (nvals > 100) {
				med3(a, ii, lo, hi, stride)
			} else {
				lo + stride*(math.floor(nvals*java.lang.Math.random()).asInstanceOf[Int])
			}
			var v = a(im)
			var iv = ii(im)
			var done = false
			var i = lo - stride
			var j = hi 
			while (! done) { 
				i += stride
				j -= stride
				while ((hi-i)*sstride > sstride*stride && ((a(i) <= v) && ((a(i) < v) || ii(i) <= iv))) {i += stride}
				while (                                   ((a(j) >= v) && ((a(j) > v) || ii(j) > iv)))  {j -= stride}
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

    def quickSort2(a:Array[Double], ii:Array[Int], lo:Int, hi:Int, stride:Int, nthreads:Int):Unit = {
            if ((hi - lo)/stride > 0) {
                if ((hi - lo)/stride <= 16) {
                    isort(a, ii, lo, hi, stride)
                } else {
                    val ip = partition(a, ii, lo, hi, stride)
                    if (nthreads > 1 && (hi-lo)/stride > 400) {
                        var done0 = false
                        var done1 = false
                        future { quickSort2(a, ii, lo, ip, stride, nthreads/2); done0 = true }
                        future { quickSort2(a, ii, ip, hi, stride, nthreads/2); done1 = true }
                        while (!done0 || !done1) {Thread.`yield`}
                    } else {
                        quickSort2(a, ii, lo, ip, stride, nthreads/2)
                        quickSort2(a, ii, ip, hi, stride, nthreads/2)
                    }
                }
            }
    }

    def isort(a:Array[Double], ii:Array[Int], lo:Int, hi:Int, stride:Int):Unit = {
            var i = lo
            while (i != hi) {
                var j = i+stride
                var imin = i
                var vmin = a(i)
                while (j != hi) {
                    if (a(j) <= vmin && ((a(j) < vmin) || ii(j) < ii(imin))) {
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

    def med3(a:Array[Double], ii:Array[Int], lo:Int, hi:Int, stride:Int):Int = {
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
            if ((v2 >= v1) && ((v2 > v1) || ii2 > ii1)) {
                if ((v3 >= v2) && ((v3 > v2) || ii3 > ii2)) i2 else {
                    if ((v3 >= v1) && ((v3 > v1) || ii3 > ii1)) i3 else i1
                }
            } else {
                if ((v3 >= v1) && ((v3 > v1) || ii3 > ii1)) i1 else {
                    if ((v3 >= v2) && ((v3 > v2) || ii3 > ii2)) i3 else i2
                }
            }
    }

    def med9(a:Array[Double], ii:Array[Int], lo:Int, hi:Int, stride:Int):Int = {
            val i1 = med3(a, ii, lo, hi, stride)
            val i2 = med3(a, ii, lo, hi, stride)
            val i3 = med3(a, ii, lo, hi, stride)
            val v1 = a(i1)
            val v2 = a(i2)
            val v3 = a(i3)
            val ii1 = ii(i1)
            val ii2 = ii(i2)
            val ii3 = ii(i3)
            if ((v2 >= v1) && ((v2 > v1) || ii2 > ii1)) {
                if ((v3 >= v2) && ((v3 > v2) || ii3 > ii2)) i2 else {
                    if ((v3 >= v1) && ((v3 > v1) || ii3 > ii1)) i3 else i1
                }
            } else {
                if ((v3 >= v1) && ((v3 > v1) || ii3 > ii1)) i1 else {
                    if ((v3 >= v2) && ((v3 > v2) || ii3 > ii2)) i3 else i2
                }
            }
    }

    def partition(a:Array[Double], ii:Array[Int], lo:Int, hi:Int, stride:Int):Int = {
            val sstride = math.signum(stride)
            val nvals = (hi - lo)/stride
            val im = if (nvals > 600) {
                med9(a, ii, lo, hi, stride)   
            } else if (nvals > 100) {
                med3(a, ii, lo, hi, stride)
            } else {
                lo + stride*(math.floor(nvals*java.lang.Math.random()).asInstanceOf[Int])
            }
            var v = a(im)
            var iv = ii(im)
            var done = false
            var i = lo - stride
            var j = hi 
            while (! done) { 
                i += stride
                j -= stride
                while ((hi-i)*sstride > sstride*stride && ((a(i) <= v) && ((a(i) < v) || ii(i) <= iv))) {i += stride}
                while (                                   ((a(j) >= v) && ((a(j) > v) || ii(j) > iv)))  {j -= stride}
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


	def quickSort2(a:Array[Int], ii:Array[Int], lo:Int, hi:Int, stride:Int, nthreads:Int):Unit = {
			if ((hi - lo)/stride > 0) {
				if ((hi - lo)/stride <= 16) {
					isort(a, ii, lo, hi, stride)
				} else {
					val ip = partition(a, ii, lo, hi, stride)
					if (nthreads > 1 && (hi-lo)/stride > 400) {
						var done0 = false
						var done1 = false
						future { quickSort2(a, ii, lo, ip, stride, nthreads/2); done0 = true }
						future { quickSort2(a, ii, ip, hi, stride, nthreads/2); done1 = true }
						while (!done0 || !done1) {Thread.`yield`}
					} else {
						quickSort2(a, ii, lo, ip, stride, nthreads/2)
						quickSort2(a, ii, ip, hi, stride, nthreads/2)
					}
				}
			}
	}

	def isort(a:Array[Int], ii:Array[Int], lo:Int, hi:Int, stride:Int):Unit = {
			var i = lo
			while (i != hi) {
				var j = i+stride
				var imin = i
				var vmin = a(i)
				while (j != hi) {
					if (a(j) <= vmin && ((a(j) < vmin) || ii(j) < ii(imin))) {
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

	def med3(a:Array[Int], ii:Array[Int], lo:Int, hi:Int, stride:Int):Int = {
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
			if ((v2 >= v1) && ((v2 > v1) || ii2 > ii1)) {
				if ((v3 >= v2) && ((v3 > v2) || ii3 > ii2)) i2 else {
					if ((v3 >= v1) && ((v3 > v1) || ii3 > ii1)) i3 else i1
				}
			} else {
				if ((v3 >= v1) && ((v3 > v1) || ii3 > ii1)) i1 else {
					if ((v3 >= v2) && ((v3 > v2) || ii3 > ii2)) i3 else i2
				}
			}
	}

	def med9(a:Array[Int], ii:Array[Int], lo:Int, hi:Int, stride:Int):Int = {
			val i1 = med3(a, ii, lo, hi, stride)
			val i2 = med3(a, ii, lo, hi, stride)
			val i3 = med3(a, ii, lo, hi, stride)
			val v1 = a(i1)
			val v2 = a(i2)
			val v3 = a(i3)
			val ii1 = ii(i1)
			val ii2 = ii(i2)
			val ii3 = ii(i3)
			if ((v2 >= v1) && ((v2 > v1) || ii2 > ii1)) {
				if ((v3 >= v2) && ((v3 > v2) || ii3 > ii2)) i2 else {
					if ((v3 >= v1) && ((v3 > v1) || ii3 > ii1)) i3 else i1
				}
			} else {
				if ((v3 >= v1) && ((v3 > v1) || ii3 > ii1)) i1 else {
					if ((v3 >= v2) && ((v3 > v2) || ii3 > ii2)) i3 else i2
				}
			}
	}

	def partition(a:Array[Int], ii:Array[Int], lo:Int, hi:Int, stride:Int):Int = {
			val sstride = math.signum(stride)
			val nvals = (hi - lo)/stride
			val im = if (nvals > 600) {
				med9(a, ii, lo, hi, stride)	  
			} else if (nvals > 100) {
				med3(a, ii, lo, hi, stride)
			} else {
				lo + stride*(math.floor(nvals*java.lang.Math.random()).asInstanceOf[Int])
			}
			var v = a(im)
			var iv = ii(im)
			var done = false
			var i = lo - stride
			var j = hi 
			while (! done) { 
				i += stride
				j -= stride
				while ((hi-i)*sstride > sstride*stride && ((a(i) <= v) && ((a(i) < v) || ii(i) <= iv))) {i += stride}
				while (                                   ((a(j) >= v) && ((a(j) > v) || ii(j) > iv)))  {j -= stride}
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
				var m = off + (len >> 1) 
				if (len > 30) {
					var l = off
					var n = off + len - 1
					if (len > 300) { 
						val s = len / 8
						l = med3(l, l+s, l+2*s)
						m = med3(m-s, m, m+s)
						n = med3(n-2*s, n-s, n)
					}
					m = med3(l, m, n) 
				}

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

				val n = off + len
				var s = math.min(a-off, b-a)
				vecswap(off, b-s, s)
				s = math.min(d-c, n-d-1)
				vecswap(b, n-s, s)

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
/*		val xa = new myMedian[Float]
		import myOrdering._
		val v = xa.median(a.data)
		println("%f (%f %f %f)" format (v, a(0), a(1), a(2))) */
	}
}


class myMedian[@specialized(Float, Double, Int) T] {
	def median(a:Array[T])(implicit ord:myOrdering[T]):T = {
		val v1 = a(0)
		val v2 = a(1)
		val v3 = a(2)
		if (ord.gt(v2, v1)) {
			if (ord.gt(v3, v2)) v2 else {
				if (ord.gt(v3, v1)) v3 else v1
			}
		} else {
			if (ord.gt(v3, v1)) v1 else {
				if (ord.gt(v3, v2)) v3 else v2
			}
		}
	}
}


trait myOrdering[T] {    
	outer => 

	def compare(x: T, y: T): Int

	def lteq(x: T, y: T): Boolean = compare(x, y) <= 0

	def gteq(x: T, y: T): Boolean = compare(x, y) >= 0

	def lt(x: T, y: T): Boolean = compare(x, y) < 0

	def gt(x: T, y: T): Boolean = compare(x, y) > 0

	def equiv(x: T, y: T): Boolean = compare(x, y) == 0

	def max(x: T, y: T): T = if (gteq(x, y)) x else y

	def min(x: T, y: T): T = if (lteq(x, y)) x else y

	def reverse: myOrdering[T] = new myOrdering[T] {
		override def reverse = outer
		def compare(x: T, y: T) = outer.compare(y, x)
	}
}

object myOrdering {
	def apply[T](implicit ord: myOrdering[T]) = ord

	trait FloatOrdering extends myOrdering[Float] {
		outer =>

		def compare(x: Float, y: Float) = java.lang.Float.compare(x, y)

		@inline override def lteq(x: Float, y: Float): Boolean = x <= y
		@inline override def gteq(x: Float, y: Float): Boolean = x >= y
		@inline override def lt(x: Float, y: Float): Boolean = x < y
		@inline override def gt(x: Float, y: Float): Boolean = x > y
		@inline override def equiv(x: Float, y: Float): Boolean = x == y
		@inline override def max(x: Float, y: Float): Float = math.max(x, y)
		@inline override def min(x: Float, y: Float): Float = math.min(x, y)

		override def reverse: myOrdering[Float] = new FloatOrdering {
			override def reverse = outer
			@inline override def compare(x: Float, y: Float) = outer.compare(y, x)

			@inline override def lteq(x: Float, y: Float): Boolean = outer.lteq(y, x)
			@inline override def gteq(x: Float, y: Float): Boolean = outer.gteq(y, x)
			@inline override def lt(x: Float, y: Float): Boolean = outer.lt(y, x)
			@inline override def gt(x: Float, y: Float): Boolean = outer.gt(y, x)
		}
	}
	implicit object Float extends FloatOrdering
}


