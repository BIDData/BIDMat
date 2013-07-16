package BIDMat
import scala.collection.mutable.HashMap
import MatFunctions._
import edu.berkeley.bid.CUMAT

class Dict(val cstr:CSMat) { 

  val length = cstr.length
  
  var counts:DMat = null
  
  var hash:HashMap[String,Int] = null

  def makeHash:HashMap[String, Int] = { 
    if (hash.asInstanceOf[AnyRef] == null) { 
      hash = new HashMap[String, Int]()
      var i = 0
      while (i < cstr.length) { 
        hash(cstr.data(i)) = i
        i += 1
      } 
    }
    hash
  }
  
  def --> (b:Dict):IMat = { 
    val bhash = b.makeHash
    val out = IMat(length, 1)
    var i = 0
    while (i < length) { 
      if (bhash.contains(cstr(i))) { 
        out.data(i) = bhash(cstr(i))
      } else { 
        out.data(i) = -1
      }
      i += 1
    }
    out
  }
  
  def apply(s:String) = {
    makeHash
    hash.get(s).get
  }
  
 def apply(i:Int) = {
    cstr(i)
  }
  
  def trim(thresh:Int):Dict = {
    val ii = find(counts >= thresh.toDouble)
    Dict(cstr(ii), DMat(counts(ii)))
  }

}

object Dict {
  
  def apply(cstr:CSMat, counts:DMat):Dict = {
    val out = new Dict(cstr)
    out.counts = counts
    out
  }
  def apply(cstr:CSMat, counts:IMat):Dict = Dict(cstr, DMat(counts))
  
  def apply(cstr:CSMat, counts:DMat, thresh:Int):Dict = {
    val ii = find(counts >= thresh.toDouble)
    val out = new Dict(cstr(ii))
    out.counts = counts(ii)
    out
  }
  def apply(cstr:CSMat, counts:IMat, thresh:Int):Dict = Dict(cstr, DMat(counts), thresh)
  
  def apply(b:BMat, counts:DMat):Dict = {
    val out = new Dict(CSMat(b))
    out.counts = counts
    out
  }
  
  def apply(b:BMat):Dict = {
    val out = new Dict(CSMat(b))
    out.counts = null
    out
  }
  
  def apply(b:BMat, counts:IMat):Dict = Dict(b, DMat(counts))
  
  def apply(b:BMat, counts:DMat, thresh:Int):Dict = {
    val ii = find(counts >= thresh.toDouble)
    val out = new Dict(CSMat(b(?,ii)))
    out.counts = counts(ii)
    out
  }
  def apply(b:BMat, counts:IMat, thresh:Int):Dict = Dict(b, DMat(counts), thresh)
  
  def apply(cstr:CSMat, counts:DMat, h:HashMap[String,Int]):Dict = {
    val out = new Dict(cstr)
    out.counts = counts
    out.hash = h
    out
  }
    
  def _union(dicts:Dict*):HashMap[String, Int] = {
    val hash = new HashMap[String, Int]()
    var nelems = 0
    dicts.foreach((d:Dict) => {
    	var i = 0 
    	while (i < d.length) {
    		if (! hash.contains(d.cstr(i))) {
    			hash(d.cstr(i)) = nelems
    			nelems += 1
    		}
    		i += 1
    	}
    })
    hash
  }
  
  def getCSMat(h:HashMap[String,Int]):CSMat = {
    val out = CSMat(h.size, 1)
    h.foreach{case (s:String, i:Int) => out.data(i) = s}
    out
  }
    
  def flatten(d1:Dict):Dict = {
    val h = _union(d1)
    val d = Dict(getCSMat(h), null, h)
    val d1d = d1 --> d
    d.counts = accum(d1d, d1.counts, d.length, 1)
    d
  }
  
  def union(dd:Dict*):Dict = {
  	val h = _union(dd:_*)
  	val d = Dict(getCSMat(h), null, h) 	
  	for (i <- 0 until dd.length) {
  		val d1d = dd(i) --> d  
  		if (i == 0) {
  			d.counts = accum(d1d, dd(i).counts, d.length, 1)
  		}	else {
  			d.counts += accum(d1d, dd(i).counts, d.length, 1)
  		}
    }    
    d
  }
  
  def union3(d1:Dict, d2:Dict):(Dict, IMat, IMat) = {
  	val h = _union(d1, d2)
    val d = Dict(getCSMat(h), null, h)
    val d1d = d1 --> d
    val d2d = d2 --> d
    d.counts = accum(d1d, d1.counts, d.length, 1) +
               accum(d2d, d2.counts, d.length, 1);
    (d, d1d, d2d)
  }

  def union4(d1:Dict, d2:Dict, d3:Dict):(Dict, IMat, IMat, IMat) = {
  	val h = _union(d1, d2, d3)
    val d = Dict(getCSMat(h), null, h)
    val d1d = d1 --> d
    val d2d = d2 --> d
    val d3d = d3 --> d
    d.counts = accum(d1d, d1.counts, d.length, 1) +
               accum(d2d, d2.counts, d.length, 1) +
               accum(d3d, d3.counts, d.length, 1);
    (d, d1d, d2d, d3d)
  }
  
  def union5(d1:Dict, d2:Dict, d3:Dict, d4:Dict):(Dict, IMat, IMat, IMat, IMat) = {
    val h = _union(d1, d2, d3, d4)
    val d = Dict(getCSMat(h), null, h)
    val d1d = d1 --> d
    val d2d = d2 --> d
    val d3d = d3 --> d
    val d4d = d4 --> d
    d.counts = (accum(d1d, d1.counts, d.length, 1) +
               accum(d2d, d2.counts, d.length, 1)) +
               (accum(d3d, d3.counts, d.length, 1) +
               accum(d4d, d4.counts, d.length, 1));
    (d, d1d, d2d, d3d, d4d)
  }
  
  def union6(d1:Dict, d2:Dict, d3:Dict, d4:Dict, d5:Dict):(Dict, IMat, IMat, IMat, IMat, IMat) = {
    val h = _union(d1, d2, d3, d4, d5)
    val d = Dict(getCSMat(h), null, h)
    val d1d = d1 --> d
    val d2d = d2 --> d
    val d3d = d3 --> d
    val d4d = d4 --> d
    val d5d = d5 --> d
    d.counts = (accum(d1d, d1.counts, d.length, 1) +
               accum(d2d, d2.counts, d.length, 1)) +
               (accum(d3d, d3.counts, d.length, 1) +
               accum(d4d, d4.counts, d.length, 1)) +
               accum(d5d, d5.counts, d.length, 1);
    (d, d1d, d2d, d3d, d4d, d5d)
  }
  
  def union7(d1:Dict, d2:Dict, d3:Dict, d4:Dict, d5:Dict, d6:Dict):(Dict, IMat, IMat, IMat, IMat, IMat, IMat) = {
    val h = _union(d1, d2, d3, d4, d5, d6)
    val d = Dict(getCSMat(h), null, h)
    val d1d = d1 --> d
    val d2d = d2 --> d
    val d3d = d3 --> d
    val d4d = d4 --> d
    val d5d = d5 --> d
    val d6d = d6 --> d
    d.counts = (accum(d1d, d1.counts, d.length, 1) +
               accum(d2d, d2.counts, d.length, 1)) +
               (accum(d3d, d3.counts, d.length, 1) +
               accum(d4d, d4.counts, d.length, 1)) +
               (accum(d5d, d5.counts, d.length, 1) +                              
               accum(d6d, d6.counts, d.length, 1));
    (d, d1d, d2d, d3d, d4d, d5d, d6d)
  }

}

class IDict(val grams:IMat) {

  val length = grams.nrows
  
  var counts:DMat = null
  
  var sortedMat:IMat = null
  
  var perm:IMat = null

  def makeSorted:IMat = { 
    if (sortedMat.asInstanceOf[AnyRef] == null) {
    	sortedMat = grams.copy
    	perm = icol(0->grams.nrows)
    	IDict.sortlex2or3cols(sortedMat, perm) 
    }
    sortedMat
  }
  
  def cmp(a:IMat, b:IMat, ia:Int, ib:Int):Int = {
    var i = 0
    var retv = 0
    while (retv == 0 && i < a.ncols) {
    	if (a(ia, i) > b(ib, i)) {
    		retv = 1
    	} else if (a(ia, i) < b(ib, i)) {
    		retv = -1
    	}
    	i += 1
    }
    retv
  }
  
  def --> (b:IDict):IMat = { 
    makeSorted
    b.makeSorted
    val out = -iones(length, 1)
    var i = 0
    var ito = 0
    while (i < length && ito < b.length) { 
      val xx = cmp(sortedMat, b.sortedMat, i, ito)
      if (xx == 0) {
        out(perm(i)) = b.perm(ito)
        i += 1
        ito += 1
      } else if (xx > 0) {
        ito += 1
      } else {
        i += 1
      }
    }
    out
  }
  
  def trim(thresh:Int):IDict = {
    val ii = find(counts >= thresh.toDouble)
    IDict(grams(ii), DMat(counts(ii)))
  }

}

object IDict {
  
  var useGPUsort = true;
  
  def apply(grams:IMat, counts:DMat):IDict = {
    val out = new IDict(grams)
    out.counts = counts
    out
  }
  def apply(grams:IMat, counts:IMat):IDict = IDict(grams, DMat(counts))
  
  def apply(grams:IMat, counts:DMat, thresh:Int):IDict = {
    val ii = find(counts >= thresh.toDouble)
    val out = new IDict(grams(ii,?))
    out.counts = counts(ii)
    out
  }
  def apply(grams:IMat, counts:IMat, thresh:Int):IDict = IDict(grams, DMat(counts), thresh)
  
  def dictFromData(grams:IMat, counts:DMat, countsort:Boolean):IDict = {
    val (outy, ia, ib) = IDict.uniquerows(grams)
    val countsy = accum(ib, if (counts == null) drow(1.0) else counts, outy.nrows, 1)
    if (countsort) {    	
    	val (countsz, ip) = GMat.sortdown2(countsy)
    	IDict(outy(ip, ?), countsz)
    } else {
      IDict(outy, countsy)
    }
  }
  
  def dictFromData(grams:IMat):IDict = dictFromData(grams, null, false)
  
  def sortlex2or3cols(mat:IMat, inds:IMat) = _sortlex2or3cols(mat, inds, true) 
  
  def _sortlex2or3cols(mat:IMat, inds:IMat, asc:Boolean) {
    import MatFunctions._
  	if (if (useGPUsort && Mat.hasCUDA > 0) {
  		val (dmy, freebytes, allbytes) = SciFunctions.GPUmem
  		if ((mat.length+inds.length)*12L < freebytes) {
  			if (mat.ncols == 2) {
  				GIMat.i2lexsortGPU(mat, inds, asc)
  				false
  			} else if (mat.ncols == 3) {
  				GIMat.i3lexsortGPU(mat, inds, asc)
  				false
  			} else true
  		} else true
  	} else true) {
  		val perm = IMat.sortlex(mat, asc) 
  		val indsp = inds(perm)
  		inds <-- indsp
  		val matp = mat(perm, ?)
  		mat <-- matp
  	}
  }
  
    
  def sortlexfast(mat:IMat, asc:Boolean):IMat = {
    import MatFunctions._
  	if (useGPUsort && Mat.hasCUDA > 0 && {
  	  val (dmy, freebytes, allbytes) = SciFunctions.GPUmem; 
  	  (mat.nrows*(mat.ncols+1)*12L < freebytes)
  	  }) 
  	{
  		val inds = icol(0->mat.nrows)
  		val tmat = mat.copy
  		if (mat.ncols == 2) {
  			GIMat.i2lexsortGPU(tmat, inds, asc)
  			inds
  		} else if (mat.ncols == 3) {
  			GIMat.i3lexsortGPU(tmat, inds, asc)
  			inds
  		} else IMat.sortlex(mat, asc) 
  	} else IMat.sortlex(mat, asc)
  }
  
  def uniquerows(a:IMat):(IMat, IMat, IMat) = {
    val iss = IMat.newOrCheckIMat(a.nrows, 1, null, a.GUID, "Dict.uniquerows".hashCode)
    val sortv = IMat.newOrCheckIMat(a.nrows, a.ncols, null, a.GUID, "Dict.uniquerows_1".hashCode)
    sortv <-- a
    var i = 0; while (i < iss.nrows) {iss(i) = i; i += 1}
    sortlex2or3cols(sortv, iss)
    def compeq(i:Int, j:Int):Boolean = {
      var k:Int = 0;
      while (k < a.ncols && (a.data(i+k*a.nrows) == a.data(j+k*a.nrows))) {
        k += 1
      }
      if (k == a.ncols) true
      else false
    }
    val iptrs = IMat.newOrCheckIMat(a.nrows, 1, null, a.GUID, "Dict.uniquerows_2".hashCode)
    var lastpos = 0
    iptrs.data(iss.data(0)) = lastpos
    i = 1
    while (i < iss.length) {
      if (!compeq(iss.data(i-1), iss.data(i))) {
        lastpos += 1
      }
      iptrs.data(iss.data(i)) = lastpos
      i += 1
    }
    val bptrs = IMat.newOrCheckIMat(lastpos+1, 1, null, a.GUID, "Dict.uniquerows_3".hashCode)
    val outv = IMat.newOrCheckIMat(lastpos+1, a.ncols, null, a.GUID, "Dict.uniquerows_4".hashCode)
    i = iss.length
    while (i > 0) {
      bptrs.data(iptrs.data(i-1)) = i-1
      i = i - 1
    }
    (a(bptrs, ?), bptrs, iptrs)    
  }  
  
  def lexcomp(a:IMat, b:IMat):(Int, Int) => Int = {
  	val aa = a.data
  	val bb = b.data
  	val nra = a.nrows
  	val nrb = b.nrows
  	val nc = a.ncols
  	(i:Int, j:Int) => {
  		var k = 0
  		while (k < nc && aa(i+k*nra) == bb(j+k*nrb)) {
  			k += 1
  		}
  		if (k == nc) {
  			0
  		} else {
  			if (aa(i+k*nra) < bb(j+k*nrb)) {
  				-1
  			} else {
  				1
  			}
  		}
  	}
  }
  
  def copyrow(a:IMat, i:Int, b:IMat, j:Int) = {
    var k = 0 
    while (k < a.ncols) {
      b.data(j + k*b.nrows) = a.data(i + k*a.nrows)
      k += 1
    }
  }
  
  def treeMerge(aa:Array[IDict]):IDict = {
    var n = 1
    var ll = 0
    while (n <= aa.length) {n *= 2; ll += 1}
    val tree = new Array[IDict](ll)
    var i = 0
    while (i < aa.length) {
      var dd = aa(i)
      var j = 0 
      while (tree(j) != null) {
        dd = merge2(tree(j), dd)
        tree(j) = null
        j += 1
      }
      tree(j) = dd
      i += 1
    }
    var j = ll-1
    var dd = tree(j)
    while (j > 0) {
    	j -= 1
    	if (tree(j) != null) dd = merge2(tree(j), dd)
    }
    dd
  }
  
  def treeAdd(x:IDict, tree:Array[IDict]) = {
    if (x != null) {
    	var dd = x
    	var j = 0 
    	while (tree(j) != null) {
    		dd = merge2(tree(j), dd)
    		tree(j) = null
    		j += 1
    	}
    	tree(j) = dd
    }
  }
  
  def treeFlush(tree:Array[IDict]):IDict = {
    var j = 0
    var dd:IDict = null
    while (j < tree.length) {
    	if (tree(j) != null) {
    	  if (dd != null) {
    	  	dd = merge2(tree(j), dd)
    	  } else {
    	    dd = tree(j)
    	  }
    	  tree(j) = null
    	}
    	j += 1
    }
    dd
  }
  
  def merge2(ad:IDict, bd:IDict):(IDict) = {
    val a = ad.grams
    val ac = ad.counts
    val b = bd.grams
    val bc = bd.counts
    var i = 0
    var j = 0
    var nout = 0
    val ccomp = lexcomp(a, b)
    while (i < a.nrows && j < b.nrows) {
      val c = ccomp(i,j) 
      if (c <= 0) {
        i += 1
      } 
      if (c >= 0) {
        j += 1
      } 
      nout += 1
    }
    if (i < a.nrows) {
      nout += a.nrows - i
    }
    if (j < b.nrows) {
      nout += b.nrows - j
    }
    val out = IMat.newOrCheckIMat(nout, a.ncols, null, a.GUID, "Dict.union2".hashCode)
    val cout = DMat.newOrCheckDMat(nout, 1, null, a.GUID, "Dict.union2_1".hashCode)
    i = 0
    j = 0
    nout = 0
    while (i < a.nrows && j < b.nrows) {
      val c = ccomp(i,j) 
      if (c <= 0) {
      	copyrow(a, i, out, nout)
      	cout(nout) = ac(i)
        i += 1
      }
      if (c >= 0) {
      	if (c > 0) {
      		copyrow(b, j, out, nout)
      		cout(nout) = bc(j)
      	} else {
      		cout(nout) += bc(j)
      	}
      	j += 1
      } 
      nout += 1
    }
    while (i < a.nrows) {
    	copyrow(a, i, out, nout)
    	cout(nout) = ac(i)
      i += 1
      nout += 1
    }
    while (j < b.nrows) {
    	copyrow(b, j, out, nout)
    	cout(nout) = bc(j)
      j += 1
      nout += 1
    }
    IDict(out, cout)
  }
  
  def union(dicts:Array[IDict]):IDict = {
    var totl = 0
    var somed:IDict = null
    dicts.foreach((d:IDict) => {
      if (d != null) {
        totl += d.length
        somed = d
      }
    })
    val outx = IMat(totl, somed.grams.ncols)
    val countsx = DMat(totl, 1)
    var totx = 0
    dicts.foreach((d:IDict) => {
      if (d != null) {
      	outx(totx->(totx+d.length),?) = d.grams
      	countsx(totx->(totx+d.length),0) = d.counts
      	totx += d.length
      }
      })
    dictFromData(outx, countsx, true)
  }
  
  def union(dicts:IDict*):IDict = union(dicts.toArray)
  
  def union2(d1:IDict, d2:IDict):(IDict, IMat, IMat) = {
  	val d = union(d1, d2)
    val d1d = d1 --> d
    val d2d = d2 --> d
    (d, d1d, d2d)
  }

  def union3(d1:IDict, d2:IDict, d3:IDict):(IDict, IMat, IMat, IMat) = {
  	val d = union(d1, d2, d3)
    val d1d = d1 --> d
    val d2d = d2 --> d
    val d3d = d3 --> d
    (d, d1d, d2d, d3d)
  }
  
  def union4(d1:IDict, d2:IDict, d3:IDict, d4:IDict):(IDict, IMat, IMat, IMat, IMat) = {
    val d = union(d1, d2, d3, d4)
    val d1d = d1 --> d
    val d2d = d2 --> d
    val d3d = d3 --> d
    val d4d = d4 --> d
    (d, d1d, d2d, d3d, d4d)
  }
  
  def union5(d1:IDict, d2:IDict, d3:IDict, d4:IDict, d5:IDict):(IDict, IMat, IMat, IMat, IMat, IMat) = {
    val d = union(d1, d2, d3, d4, d5)
    val d1d = d1 --> d
    val d2d = d2 --> d
    val d3d = d3 --> d
    val d4d = d4 --> d
    val d5d = d5 --> d
    (d, d1d, d2d, d3d, d4d, d5d)
  }
  
  def union6(d1:IDict, d2:IDict, d3:IDict, d4:IDict, d5:IDict, d6:IDict):(IDict, IMat, IMat, IMat, IMat, IMat, IMat) = {
    val d = union(d1, d2, d3, d4, d5, d6)
    val d1d = d1 --> d
    val d2d = d2 --> d
    val d3d = d3 --> d
    val d4d = d4 --> d
    val d5d = d5 --> d
    val d6d = d6 --> d
    (d, d1d, d2d, d3d, d4d, d5d, d6d)
  }

}
