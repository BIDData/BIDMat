package BIDMat
import scala.collection.mutable.{Map,SynchronizedMap,HashMap}
import MatFunctions._
//import edu.berkeley.bid.CUMAT

@SerialVersionUID(100L)
class Dict(val cstr:CSMat) extends Serializable { 

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
  
  def flatten:Dict = {
    val h = Dict._union(this)
    val d = Dict(Dict.getCSMat(h), null, h)
    val d1d = this --> d
    d.counts = accum(d1d, counts, d.length, 1)
    d
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
    
  def apply(m:IMat):CSMat = {
    val out = CSMat(m.nrows,m.ncols);
    var i = 0
    while (i < m.ncols) {
      var j = 0
      while (j < m.nrows) {
        out.data(j + m.nrows * i) = cstr(m.data(j + m.nrows * i))
        j += 1
      }
      i += 1
    }
    out
  }
  
  def apply(m:CSMat):IMat = {
    val out = IMat(m.nrows,m.ncols);
    makeHash
    var i = 0
    while (i < m.ncols) {
      var j = 0
      while (j < m.nrows) {
        out.data(j + m.nrows * i) = hash.getOrElse(m.data(j + m.nrows * i), -1)
        j += 1
      }
      i += 1
    }
    out
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
    out.counts = if (counts.ncols == 1) counts else counts.t
    out
  }
  def apply(cstr:CSMat, counts:IMat):Dict = Dict(cstr, DMat(counts))
  def apply(cstr:CSMat):Dict = Dict(cstr, dones(cstr.length,1))
  
  def apply(cstr:CSMat, counts:DMat, thresh:Int):Dict = {
    val ii = find(counts >= thresh.toDouble)
    val out = new Dict(cstr(ii))
    out.counts = counts(ii)
    out
  }
  def apply(cstr:CSMat, counts:IMat, thresh:Int):Dict = Dict(cstr, DMat(counts), thresh)
  
  def apply(b:SBMat, counts:DMat):Dict = {
    val out = new Dict(CSMat(b))
    out.counts = counts
    out
  }
  
  def apply(b:SBMat):Dict = {
    val out = new Dict(CSMat(b))
    out.counts = null
    out
  }
  
  def apply(b:SBMat, counts:IMat):Dict = Dict(b, DMat(counts))
  
  def apply(b:SBMat, counts:DMat, thresh:Int):Dict = {
    val ii = find(counts >= thresh.toDouble)
    val out = new Dict(CSMat(b(?,ii)))
    out.counts = counts(ii)
    out
  }
  def apply(b:SBMat, counts:IMat, thresh:Int):Dict = Dict(b, DMat(counts), thresh)
  
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
    
  def flatten(d1:Dict):Dict = d1.flatten
  
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
  
  def treeAdd(x:Dict, tree:Array[Dict], thresh:IMat = null) = {
    if (x != null) {
    	var dd = x
    	var j = 0 
    	while (tree(j) != null) {
    		dd = union(tree(j), dd);
    		if (thresh.asInstanceOf[AnyRef] != null) dd = Dict(dd.cstr, dd.counts, thresh(j));
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
