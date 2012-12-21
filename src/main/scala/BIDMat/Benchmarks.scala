package BIDMat


object AltaVistaCrawl { 
	def tcToSMat(m:IMat, n:Int) ={
	  val ioff = Mat.ioneBased
	  var mnnz = 0
	  var source = 0
	  var ncols = 0
	  var i = 0
	  while (i < m.nrows) {
	    if (i == 0 || m(i,0) > m(i-1,0)) {
	    } else {
	      if (m(i,1) > 0) {
	        mnnz += 1
	      }
	    }
	    i += 1
	  }
	  val out = SMat(n, n, mnnz)
	  i = 0
	  var jc0 = 0
	  while (i < m.nrows) {
	    if (i == 0 || m(i,0) > m(i-1,0)) {
	      source = m(i,1)
	      while (jc0 < source) {
	        out.jc(jc0) = i + ioff
	        jc0 += 1
	      }
	    } else {
	      if (m(i,1) > 0) {
	        out.data(i) = 1
	        out.ir(i) = m(i,1) + ioff
	      }
	    }
	    i += 1
	  }
	  while (jc0 <= n) {
	  	out.jc(jc0) = i + ioff
	  	jc0 += 1
	  }
	}

}
