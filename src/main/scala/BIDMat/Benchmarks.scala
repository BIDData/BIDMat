package BIDMat
import SciFunctions._
import MatFunctions._
import java.io.File


object AltaVista { 

	
	def partition(dirname:String, fname:String, outname:String, nfiles:Int, nparts:Int,
	    dopart:Boolean, docount:Boolean, doperm:Boolean, dotranspose:Boolean, randperm:Boolean) = {
	  
	  // Randomly partition edges ready for butterfly mixing
	  val seed = 453154334
	  var maxv = 0 // 1413511394
	  val blksize = 1000000
	  var ninpart:IMat = null
	  var ff = fblank
	  
	  // Get the max node index and count edges in each partition
	  if (dopart) {
	  	println("Computing number of indices and random partition sizes")
	  	ninpart = izeros(nparts, 1)
	  	setseed(seed)
	  	for (i <- 0 until nfiles) {
	  		val mm = HMat.loadIMat((dirname + fname) format i)
	  		var j = 0 
	  		while (j < mm.nrows) {	  		  
	  		  val rr = math.floor(myrand.nextFloat * (0.9999f*nparts)).toInt
	  		  if (mm(j, 1) > 0) {
	  		    ninpart(rr) += 1
	  		    maxv = math.max(maxv, 1+math.max(mm(j,0), mm(j,1)))
	  		  }
	  		  j += 1
	  		  
	  		  /* val jend = math.min(mm.nrows, j+blksize)
	  		  ff = rand(jend-j, 1, ff)
	  			var rr = IMat(floor(ff * (0.9999f*nparts)))
	  			val iigood = find(mm(j->jend,1) > 0)
	  			ninpart += accum(rr(iigood), 1, nparts)
	  			maxv = math.max(maxv, 1+maxi(maxi(mm)).v) 
	  			j += blksize */
	  		}
	  		if (i % 3 == 0) printf(".")
	  	}
	  	printf("\nmaxv=%d\n" format maxv)	
	  	HMat.saveIMat(dirname + "ninpart.gz", ninpart)
	  } else {
	    println("loading partition sizes")
	    ninpart = HMat.loadIMat(dirname + "ninpart.gz")
	    maxv = 1413511394
	  }

	  
	  // Now count nodes
	  var nodecounts:IMat = null
	  if (docount) {
	  	nodecounts = izeros(maxv, 1)
	  	println("Computing node frequencies")
	  	for (i <- 0 until nfiles) {
	  		val mm = HMat.loadIMat((dirname + fname) format i)
	  		var j = 0 
	  		while (j < mm.nrows && mm(j,1) > 0) {
	  			nodecounts(mm(j,0)) += 1
	  			nodecounts(mm(j,1)) += 1
	  			j += 1
	  		}
	  		if (i % 3 == 0) printf(".")
	  	}	  
	  	println("\nSaving") 
	  	HMat.saveIMat(dirname + "nodecounts.gz", nodecounts)
	  } else {
	    println("loading node counts")
	    nodecounts = HMat.loadIMat(dirname + "nodecounts.gz")
	  }
	  
	  var isinv:IMat = null
	  if (doperm) {
	  	printf("Creating arrays...")	  
	  	val inds = icol(0->maxv)
	  	if (!randperm) {
	  		nodecounts ~ nodecounts * -1
	  		printf("sorting...")
	  		Sorting.quickSort2(nodecounts.data, inds.data, 0, maxv, 1)
	  	} else {
	  		printf("permuting...")
	  	  Sorting.quickSort2(rand(maxv,1).data, inds.data, 0, maxv, 1)
	  	}
	  	printf("inverting...")
	  	var i = 0
	  	while (i < maxv) {
	  		nodecounts(inds(i)) = i
	  		i += 1
	  	}	 
	  	isinv = nodecounts  
	  	printf("saving permutation\n")
	  	HMat.saveIMat(dirname + "isinv.gz", isinv)
	  } else {
	    println("loading permutation")
	    isinv = HMat.loadIMat(dirname + "isinv.gz")
	  }

	  // Now build each sampled set in turn
	  for (i <- 0 until nparts) {
	    printf("output set %d" format i)
	  	val ii = IMat(ninpart(i), 1)
	  	val jj = IMat(ninpart(i), 1)
	  	var sofar = 0
	  	setseed(seed)
	  	for (j <- 0 until nfiles) {
	  		val mm = HMat.loadIMat((dirname + fname) format j)
	  		val ici = if (dotranspose) 1 else 0
	  		var k = 0
	  		while (k < mm.nrows) {
	  			val rr = math.floor(myrand.nextFloat * (0.9999f*nparts)).toInt
	  		  if (mm(k, 1) > 0 && rr == i) {
	  		    ii(sofar) = isinv(mm(k, ici))
	  		    jj(sofar) = isinv(mm(k, 1-ici))
	  		    sofar += 1
	  		  }
	  			k += 1
	  		  
	  		  /*val jend = math.min(mm.nrows, j+blksize)
	  			val rr = floor(rand(jend-j, 1) * (0.9999f*nparts))
	  			val iigood = find((rr === i) *@ (mm(j->jend, 1) > 0))
	  			val iin = iigood + k
	  			val iout = icol(sofar->(sofar+iigood.length))
	  			ii(iout,0) = isinv(mm(iin, ici))
	  			jj(iout,0) = isinv(mm(iin, 1-ici))
	  			sofar += iigood.length	  			
	  			k += blksize */
	  		}
	  		if (j % 3 == 0) printf(".")
	  	}
	    Mat.ilexsort2(ii.data, jj.data)
	    saveAs("/big/Yahoo/G2/test.mat",ii,"ii",jj,"jj")
	    val avals = ones(ii.length, 1).data
	    sofar = SparseMat.remdups(ii.data, jj.data, avals)
	  	val ss = new SMat(maxv, maxv, sofar, SparseMat.incInds(jj.data), new Array[Int](maxv+1), avals)
	    SparseMat.compressInds(ii.data, maxv, ss.jc, sofar)
	    ss.sparseTrim
	  	HMat.saveSMat((dirname + outname) format i, ss)
	  	ss.check
	  	println("")
	  }
	}
	
	def pagerank_setup(fpath:String, size:Int, nfiles:Int, scalepath:String, iterpath:String) = {
	  printf("setting up")
	  val scale = zeros(size,1)
	  val iter = zeros(1, size)
	  val tmp1 = zeros(size, 1)
	  val tmp2 = new FMat(1, size, tmp1.data)
	  for (i <- 0 until nfiles) {
	    val ss = HMat.loadSMat(fpath format i, false)
	    scale ~ scale + sum(ss, 2, tmp1)  // add up outdegree
	    iter ~ iter + sum(ss, 1, tmp2)   // add up indegree
	    printf(".")
	  }
	  println("")
	  max(scale, 1.0f, scale)
	  scale ~ 1.0f / scale
	  val scalet = new FMat(1, size, scale.data)
	  iter ~ iter * (1.0f / sum(iter).v)
	  HMat.saveFMat(scalepath, scalet, false)
	  HMat.saveFMat(iterpath format 0, iter, false)
	}
	
  def pagerank_iter(fpath:String, size:Int, nfiles:Int, iiter:Int, scalepath:String, iterpath:String, alpha:Float) = {
		printf("iteration %d" format (iiter+1))
	  val iter = HMat.loadFMat(iterpath format iiter, false)
	  val scale = HMat.loadFMat(scalepath, false)
	  iter ~ iter *@ scale
	  val newiter = zeros(1, size)
	  var tmp = scale
	  for (i <- 0 until nfiles) {
	    val ss = HMat.loadSMat(fpath format i, false)	    
	    newiter ~ newiter + (tmp ~ iter * ss)
	    printf(".")
	  }
	  newiter ~ newiter * (1-alpha)
	  newiter ~ newiter + (alpha/size)
	  tmp ~ newiter - iter
	  val v = tmp ddot tmp
	  println("resid = %g, time = %f, gf = %f" format (math.sqrt(v/size), gflop._2, gflop._1))
	  HMat.saveFMat(iterpath format (iiter+1), newiter, false)
	}
  
  def pagerank_run(dirname:String, fname:String, nparts:Int, niter:Int) = {
    val alpha = 0.1f
    val scalename = dirname + "pagerank/scale.fmat"
    val itername = dirname + "pagerank/iter%03d.fmat"
    flip
    pagerank_setup(dirname + fname, 1413511394, nparts, scalename, itername)
    for (i <- 0 until niter) {
      pagerank_iter(dirname + fname, 1413511394, nparts, i, scalename, itername, alpha)      
    }
    val gf = gflop
    println("time=%f, gflops=%f" format (gf._2, gf._1))
  }
	
	def main(args:Array[String]):Unit = {
	  val dirname = if (args != null && args.length > 0) args(0) else "/big/Yahoo/G2/"
	  val fname =   if (args != null && args.length > 1) args(1) else "parts/out%03d.mtab.gz"
	  val outname = if (args != null && args.length > 2) args(2) else "tparts/part%03dp4.smat"
	  val nfiles =  if (args != null && args.length > 3) args(3).toInt else 179
	  val nparts =  if (args != null && args.length > 4) args(4).toInt else 8

	  // partition(dirname, fname, outname, nfiles, nparts, true, false, true, true, false)
	  pagerank_run(dirname, outname, nparts, 10)
	}

}

object Twitter { 
  
  	def main(args:Array[String]):Unit = {
	  val infname =  if (args != null && args.length > 0) args(0) else "/disk%02d/twitter/tokenized/"
	  val outfname = if (args != null && args.length > 1) args(1) else "/big/twitter/tokenized/"
	  val threshold = if (args != null && args.length > 2) args(2).toInt else 10

	  mergedicts(2011, 2013, infname, outfname, threshold)
	}
  
	def mergedicts(year1:Int, year2:Int, infname:String, outfname:String, threshold:Int) = {
  	val dd = new Array[Dict](6)
  	val md = new Array[Dict](6)
	  for (yy <- year1 to year2) {
	  	for (mm <- 1 to 12) {
	  		print("\n%d/%02d" format (yy, mm))
	  		val ff = new File(outfname + "%04d/%02d/wcount.gz" format (yy, mm))
	  		if (! ff.exists) {
	  			var here = 0
	  			var ndone = 0
	  			for (id <- 1 to 31) {
	  				var there = (here + 1) % 16
	  				var foundit = false
	  				while (!foundit && here != there) {
	  					val fname = (infname + "%04d/%02d/%02d/" format (there, yy, mm, id))
	  					val ff = new File(fname + "wcount.gz")
	  					if (ff.exists) {
	  						val bb = HMat.loadBMat(fname + "dict.gz")
	  						val cc = HMat.loadIMat(fname + "wcount.gz")
	  						dd(ndone % 6) = Dict(bb, cc, threshold)
	  						ndone = ndone + 1
	  						print("-")
	  						if (ndone % 6 == 0) {
	  							md(ndone / 6 - 1) = Dict.union(dd:_*)
	  							print("+")
	  						}
	  						foundit = true
	  						here = there
	  					}
	  					there = (there + 1) % 16
	  				}
	  			}
	  			if (ndone % 6 != 0) {
	  				md(ndone / 6) = Dict.union(dd.slice(0, ndone % 6):_*)
	  				print("+")
	  			}
	  			if (ndone > 0) {
	  				val dx = Dict.union(md.slice(0, (ndone-1)/6+1):_*)
	  				val (sv, iv) = sortdown2(dx.counts)
	  				val dxx = Dict(dx.cstr(iv), sv)
	  				HMat.saveBMat(outfname + "%04d/%02d/dict.gz" format (yy, mm), BMat(dxx.cstr))
	  				HMat.saveDMat(outfname + "%04d/%02d/wcount.gz" format (yy, mm), dxx.counts)
	  			}
//	  			println("")
	  		}
	  	}
	  }
	}
}
