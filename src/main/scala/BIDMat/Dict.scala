package BIDMat
import scala.collection.mutable.HashMap
import MatFunctions._

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
  
  var iperm:IMat = null

  def makeSorted:IMat = { 
    if (sortedMat.asInstanceOf[AnyRef] == null) { 
      perm = sortlex(grams) 
      sortedMat = grams(perm, ?)
      iperm = accum(perm, icol(0->perm.length), perm.length, 1)
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
    val bhash = b.makeSorted
    val out = -iones(length, 1)
    var i = 0
    var ito = 0
    while (i < length && ito < b.length) { 
      val xx = cmp(sortedMat, b.sortedMat, i, ito)
      if (xx == 0) {
        out(iperm(i)) = b.perm(ito)
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
  
  def apply(grams:IMat, counts:DMat):IDict = {
    val out = new IDict(grams)
    out.counts = counts
    out
  }
  def apply(grams:IMat, counts:IMat):IDict = IDict(grams, DMat(counts))
  
  def apply(grams:IMat, counts:DMat, thresh:Int):IDict = {
    val ii = find(counts >= thresh.toDouble)
    val out = new IDict(grams(ii))
    out.counts = counts(ii)
    out
  }
  def apply(grams:IMat, counts:IMat, thresh:Int):IDict = IDict(grams, DMat(counts), thresh)
  
  def dictFromData(grams:IMat, counts:DMat):IDict = {
    val (outy, ia, ib) = uniquerows(grams)
    val countsy = accum(ib, if (counts == null) drow(1.0) else counts, outy.nrows, 1)
    val (countsz, ip) = sortdown2(countsy)
    IDict(outy(ip, ?), countsz)
  }
  
  def dictFromData(grams:IMat):IDict = dictFromData(grams, null)
  
  def _union(dicts:IDict*):IDict = {
    var totl = 0
    dicts.foreach((d:IDict) => {totl += d.length})
    val outx = IMat(totl, dicts(0).grams.ncols)
    val countsx = DMat(totl, 1)
    var totx = 0
    dicts.foreach((d:IDict) => {
      outx(totx->(totx+d.length),?) = d.grams
      countsx(totx->(totx+d.length),0) = d.counts
      totx += d.length
      })
    dictFromData(outx, countsx)
  }

  def union(dd:IDict*):IDict = {
  	val d = _union(dd:_*) 	
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
  
  def union2(d1:IDict, d2:IDict):(IDict, IMat, IMat) = {
  	val d = _union(d1, d2)
    val d1d = d1 --> d
    val d2d = d2 --> d
    d.counts = accum(d1d, d1.counts, d.length, 1) +
               accum(d2d, d2.counts, d.length, 1);
    (d, d1d, d2d)
  }

  def union3(d1:IDict, d2:IDict, d3:IDict):(IDict, IMat, IMat, IMat) = {
  	val d = _union(d1, d2, d3)
    val d1d = d1 --> d
    val d2d = d2 --> d
    val d3d = d3 --> d
    d.counts = accum(d1d, d1.counts, d.length, 1) +
               accum(d2d, d2.counts, d.length, 1) +
               accum(d3d, d3.counts, d.length, 1);
    (d, d1d, d2d, d3d)
  }
  
  def union4(d1:IDict, d2:IDict, d3:IDict, d4:IDict):(IDict, IMat, IMat, IMat, IMat) = {
    val d = _union(d1, d2, d3, d4)
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
  
  def union5(d1:IDict, d2:IDict, d3:IDict, d4:IDict, d5:IDict):(IDict, IMat, IMat, IMat, IMat, IMat) = {
    val d = _union(d1, d2, d3, d4, d5)
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
  
  def union6(d1:IDict, d2:IDict, d3:IDict, d4:IDict, d5:IDict, d6:IDict):(IDict, IMat, IMat, IMat, IMat, IMat, IMat) = {
    val d = _union(d1, d2, d3, d4, d5, d6)
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
