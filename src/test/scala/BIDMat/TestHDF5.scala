package BIDMat


import MatFunctions._
import SciFunctions._
import CMat._



object TestHDF5 {
	def main(args: Array[String]) : Unit = {
			val n = 50000
			val k = 10
			val l = 1
			val a = rand(n,k)
			val b = IMat(l,n)
			val c = sprand(10,10,0.1)
			val d = CSMat(1,2)
			d(0,0) = "test"
			d(0,1) = "try"
			val fname = "d:\\sentiment\\tmp\\mtest.mat"
//			println(a.toString)
//			println(b.toString)
			saveAs(fname, d, "c")
/*			val fid = new java.io.FileInputStream(fname)
			var next:Int = 0
			var i = 0
			while (next >= 0) {
				var next = fid.read()
				if (next > 0) println("buf("+i+")=("+(next + (if (next > 127) -256 else 0))+")")
				i += 1
			}
				fid.close() */
			val e:CMat = load(fname, "c").asInstanceOf[CMat]
	    println(e(0,0).asInstanceOf[String])
		  println(e(0,1).asInstanceOf[String])
	}
}

