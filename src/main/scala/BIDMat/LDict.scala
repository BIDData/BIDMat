package BIDMat
import scala.collection.mutable.{Map,SynchronizedMap,HashMap}
import MatFunctions._
import edu.berkeley.bid.CUMAT

/**
 * Note: this is still alpha code. Needs native code infill for long sort routines. 
 */

class LDict(val grams:LMat) extends Serializable {

  val length = grams.nrows
  
  var counts:DMat = null
  
  var sortedMat:LMat = null
  
  var perm:IMat = null

  def makeSorted:LMat = { 
    this.synchronized {
    	if (sortedMat.asInstanceOf[AnyRef] == null) {
    		sortedMat = grams.copy
    		perm = icol(0->grams.nrows)
    		sortlexInds(sortedMat, perm) 
    	}
    	sortedMat
    }
  }
  
  @inline def cmp(a:LMat, b:LMat, ia:Int, ib:Int):Int = {
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
  
  def --> (b:LDict):IMat = { 
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
  
  def flatten:LDict = {
    val (fl, ii, jj) = uniquerows(grams)
    if (counts.asInstanceOf[AnyRef] == null) counts = dones(grams.nrows,1)
    val newcounts = accum(jj, counts, fl.nrows, 1)
    LDict(fl, newcounts)
  }
  
  def trim(thresh:Int):LDict = {
    val ii = find(counts >= thresh.toDouble)
    LDict(grams(ii), DMat(counts(ii)))
  }

}

object LDict {
  
  var useGPUsort = true;
  
  def apply(grams:LMat, counts:DMat):LDict = {
    val out = new LDict(grams)
    out.counts = counts
    out
  }
  def apply(grams:LMat, counts:IMat):LDict = LDict(grams, DMat(counts))
  def apply(grams:LMat):LDict = LDict(grams, null:DMat)
  
  def apply(grams:LMat, counts:DMat, thresh:Int):LDict = {
    val ii = find(counts >= thresh.toDouble)
    val out = new LDict(grams(ii,?))
    out.counts = counts(ii)
    out
  }
  def apply(grams:LMat, counts:IMat, thresh:Int):LDict = LDict(grams, DMat(counts), thresh)
  
  def dictFromData(grams:LMat, counts:DMat, countsort:Boolean):LDict = {
    val (outy, ia, ib) = uniquerows(grams)
    val countsy = accum(ib, if (counts == null) drow(1.0) else counts, outy.nrows, 1)
    if (countsort) {    	
    	val (countsz, ip) = GFunctions.sortdown2(countsy)
    	LDict(outy(ip, ?), countsz)
    } else {
      LDict(outy, countsy)
    }
  }
  
  def dictFromData(grams:LMat):LDict = dictFromData(grams, null, false)
  
 
  def lexcomp(a:LMat, b:LMat):(Int, Int) => Int = {
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
  
  def treeMerge(aa:Array[LDict], thresholds:IMat = null):LDict = {
    var n = 1
    var ll = 0
    while (n <= aa.length) {n *= 2; ll += 1}
    val tree = new Array[LDict](ll)
    var i = 0
    while (i < aa.length) {
      var dd = aa(i)
      var j = 0 
      while (tree(j) != null) {
	  val thresh = if (thresholds == null) 0 else if (thresholds.length > 1) thresholds(j) else thresholds(0);
	  dd = merge2(tree(j), dd, thresh);
	  tree(j) = null;
	  j += 1;
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
  
  def treeAdd(x:LDict, tree:Array[LDict], thresholds:IMat=null) = {
    if (x != null) {
    	var dd = x
    	var j = 0 
    	while (tree(j) != null) {
	    val thresh = if (thresholds == null) 0 else if (thresholds.length > 1) thresholds(j) else thresholds(0);
	    dd = merge2(tree(j), dd, thresh);
	    tree(j) = null;
	    j += 1;
    	}
    	tree(j) = dd
    }
  }
  
  def treeFlush(tree:Array[LDict], thresholds:IMat = null):LDict = {
    var j = 0;
    var dd:LDict = null;
    while (j < tree.length) {
    	if (tree(j) != null) {
	    if (dd != null) {
		val thresh = if (thresholds == null) 0 else if (thresholds.length > 1) thresholds(j) else thresholds(0);
    	  	dd = merge2(tree(j), dd, thresh);
	    } else {
		dd = tree(j);
	    }
	    tree(j) = null;
    	}
    	j += 1;
    }
    dd;
  }
  
  def merge2(ad:LDict, bd:LDict, threshold:Int = 0):(LDict) = {
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
    val out = LMat.newOrCheckLMat(nout, a.ncols, null, a.GUID, "Dict.union2".hashCode)
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
    LDict(out, cout, threshold)
  }
  
  def union(dicts:Array[LDict]):LDict = {
    var totl = 0
    var somed:LDict = null
    dicts.foreach((d:LDict) => {
      if (d != null) {
        totl += d.length
        somed = d
      }
    })
    val outx = LMat(totl, somed.grams.ncols)
    val countsx = DMat(totl, 1)
    var totx = 0
    dicts.foreach((d:LDict) => {
      if (d != null) {
      	outx(totx->(totx+d.length),?) = d.grams
      	countsx(totx->(totx+d.length),0) = d.counts
      	totx += d.length
      }
      })
    dictFromData(outx, countsx, true)
  }
  
  def union(dicts:LDict*):LDict = union(dicts.toArray)
  
  def union2(d1:LDict, d2:LDict):(LDict, IMat, IMat) = {
  	val d = union(d1, d2)
    val d1d = d1 --> d
    val d2d = d2 --> d
    (d, d1d, d2d)
  }

  def union3(d1:LDict, d2:LDict, d3:LDict):(LDict, IMat, IMat, IMat) = {
  	val d = union(d1, d2, d3)
    val d1d = d1 --> d
    val d2d = d2 --> d
    val d3d = d3 --> d
    (d, d1d, d2d, d3d)
  }
  
  def union4(d1:LDict, d2:LDict, d3:LDict, d4:LDict):(LDict, IMat, IMat, IMat, IMat) = {
    val d = union(d1, d2, d3, d4)
    val d1d = d1 --> d
    val d2d = d2 --> d
    val d3d = d3 --> d
    val d4d = d4 --> d
    (d, d1d, d2d, d3d, d4d)
  }
  
  def union5(d1:LDict, d2:LDict, d3:LDict, d4:LDict, d5:LDict):(LDict, IMat, IMat, IMat, IMat, IMat) = {
    val d = union(d1, d2, d3, d4, d5)
    val d1d = d1 --> d
    val d2d = d2 --> d
    val d3d = d3 --> d
    val d4d = d4 --> d
    val d5d = d5 --> d
    (d, d1d, d2d, d3d, d4d, d5d)
  }
  
  def union6(d1:LDict, d2:LDict, d3:LDict, d4:LDict, d5:LDict, d6:LDict):(LDict, IMat, IMat, IMat, IMat, IMat, IMat) = {
    val d = union(d1, d2, d3, d4, d5, d6)
    val d1d = d1 --> d
    val d2d = d2 --> d
    val d3d = d3 --> d
    val d4d = d4 --> d
    val d5d = d5 --> d
    val d6d = d6 --> d
    (d, d1d, d2d, d3d, d4d, d5d, d6d)
  }
  
  def gramDict(nuni:Int, nbi:Int, ntri:Int, du:Dict, db:LDict, dt:LDict) = {
    val c1 = du.cstr(0->nuni,0)
    val n1 = du.counts(0->nuni, 0)
    val c2 = du.cstr(IMat(db.grams(0->nbi, 0)), 0) + csrow(" ") + du.cstr(IMat(db.grams(0->nbi, 1)), 0);
    val n2 = db.counts(0->nbi, 0)
    val c3 = du.cstr(IMat(dt.grams(0->ntri, 0)), 0) + csrow(" ") + du.cstr(IMat(dt.grams(0->ntri, 1)), 0) + csrow(" ") + du.cstr(IMat(dt.grams(0->ntri, 2)), 0)
    val n3 = dt.counts(0->ntri, 0)
    Dict(c1 on c2 on c3, n1 on n2 on n3)
  }

}
