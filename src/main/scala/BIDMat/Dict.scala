package BIDMat
import scala.collection.mutable.{Map,SynchronizedMap,HashMap}
import MatFunctions._
import edu.berkeley.bid.CUMAT

class Dict(val cstr:CSMat) { 

  val length = cstr.length
  
  var counts:DMat = null
  
  var hash:HashMap[String,Int] = null

  def makeHash:HashMap[String, Int] = { 
    this.synchronized {
    	if (hash.asInstanceOf[AnyRef] == null) { 
    		hash = new HashMap[String, Int] 
    		var i = 0
    		while (i < cstr.length) { 
    			hash(cstr.data(i)) = i
    			i += 1
    		} 
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
  
  def apply(s:String):Int = {
    makeHash
    hash.getOrElse(s, -1)
  }
  
 def apply(i:Int) = {
    cstr(i)
  }
 
  def apply(x:Range) = {
    cstr(x)
  }
 
 def count(s:String):Double = {
   makeHash
   val v = hash.getOrElse(s, -1)
   if (v >= 0) counts(v) else 0.0
 }
 
 def count(i:Int):Double = {
   counts(i) 
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
  
  def treeAdd(x:Dict, tree:Array[Dict]) = {
    if (x != null) {
    	var dd = x
    	var j = 0 
    	while (tree(j) != null) {
    		dd = union(tree(j), dd)
    		tree(j) = null
    		j += 1
    	}
    	tree(j) = dd
    }
  }
  
  def treeFlush(tree:Array[Dict]):Dict = {
    var j = 0
    var dd:Dict = null
    while (j < tree.length) {
    	if (tree(j) != null) {
    	  if (dd != null) {
    	  	dd = union(tree(j), dd)
    	  } else {
    	    dd = tree(j)
    	  }
    	  tree(j) = null
    	}
    	j += 1
    }
    dd
  }

}

class IDict(val grams:IMat) {

  val length = grams.nrows
  
  var counts:DMat = null
  
  var sortedMat:IMat = null
  
  var perm:IMat = null

  def makeSorted:IMat = { 
    this.synchronized {
    	if (sortedMat.asInstanceOf[AnyRef] == null) {
    		sortedMat = grams.copy
    		perm = icol(0->grams.nrows)
    		sortlexInds(sortedMat, perm) 
    	}
    	sortedMat
    }
  }
  
  @inline def cmp(a:IMat, b:IMat, ia:Int, ib:Int):Int = {
    var i = 0
    var retv = 0
    val anr = a.nrows
    val bnr = b.nrows
    while (retv == 0 && i < a.ncols) {
    	if (a.data(ia + i*anr) > b.data(ib + i*bnr)) {
    		retv = 1
    	} else if (a.data(ia + i*anr) < b.data(ib + i*bnr)) {
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
        out.data(perm(i)) = b.perm.data(ito)
        i += 1
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
  def apply(grams:IMat):IDict = IDict(grams, null:DMat)
  
  def apply(grams:IMat, counts:DMat, thresh:Int):IDict = {
    val ii = find(counts >= thresh.toDouble)
    val out = new IDict(grams(ii,?))
    out.counts = counts(ii)
    out
  }
  def apply(grams:IMat, counts:IMat, thresh:Int):IDict = IDict(grams, DMat(counts), thresh)
  
  def dictFromData(grams:IMat, counts:DMat, countsort:Boolean):IDict = {
    val (outy, ia, ib) = uniquerows(grams)
    val countsy = accum(ib, if (counts == null) drow(1.0) else counts, outy.nrows, 1)
    if (countsort) {    	
    	val (countsz, ip) = GMat.sortdown2(countsy)
    	IDict(outy(ip, ?), countsz)
    } else {
      IDict(outy, countsy)
    }
  }
  
  def dictFromData(grams:IMat):IDict = dictFromData(grams, null, false)
  
 
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
  
  def gramDict(nuni:Int, nbi:Int, ntri:Int, du:Dict, db:IDict, dt:IDict) = {
    val c1 = du.cstr(0->nuni,0)
    val n1 = du.counts(0->nuni, 0)
    val c2 = du.cstr(db.grams(0->nbi, 0), 0) + csrow(" ") + du.cstr(db.grams(0->nbi, 1), 0)
    val n2 = db.counts(0->nbi, 0)
    val c3 = du.cstr(dt.grams(0->ntri, 0), 0) + csrow(" ") + du.cstr(dt.grams(0->ntri, 1), 0) + csrow(" ") + du.cstr(dt.grams(0->ntri, 2), 0)
    val n3 = dt.counts(0->ntri, 0)
    Dict(c1 on c2 on c3, n1 on n2 on n3)
  }

}
