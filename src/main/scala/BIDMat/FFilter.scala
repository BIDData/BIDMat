package BIDMat;
import BIDMat.MatFunctions._;
import edu.berkeley.bid.CBLAS._
import SciFunctions._

//
// Basic CPU convolutional Filter class.
//
// Works in any dimension, but the last three (minor) dimensions should be HWC. 
// Matrix blocks for each Filter dimension are stored in IO order (input-major, output-minor).
//
// e.g. a 4D data tensor would be NHWC (Channel minor, N(minibatch number) major).
// a 4D filter block would be HWIO
// 
// There are three algorithms:
// * Element-wise - pure Scala implementation
// * im2col - copy input elements into a im2col matrix and do a single multiply.
// * 1x1 convolutions. For a kxk filter, perform k^2 1x1 convolutions each with a single matrix multiply. 


@SerialVersionUID(100L)
class FFilter(inDims0:IMat, outDims0:IMat, stride0:IMat, pad0:IMat, data0:Array[Float]) extends
FND((inDims0 *@ outDims0).data, data0) with Filter {

	override val inDims = inDims0;
	override val outDims = outDims0;
	override val stride = if (stride0.asInstanceOf[AnyRef] != null) stride0 else iones(1, inDims.length);
	override val pad = if (pad0.asInstanceOf[AnyRef] != null) pad0 else izeros(1,inDims.length);

	def convolve(a:FND, omat:ND, doclear:Boolean):FND = {
		val bdims = Filter.getOutputDims(a.dims, inDims, outDims, stride, pad);
    val hmm = Filter.hashIMat(stride, Filter.hashIMat(pad));
		val b = FND.newOrCheckFND(bdims, omat, a.GUID, GUID, hmm, "convout".##);
		val apadmat = if (pad.data.exists(_ != 0)) {
      val m = FND.newOrCheckFND(a.dims + pad * 2, null, a.GUID, GUID, hmm, "convinpad".##);
      m.clear;
      _copy_padded(a, m, pad, inDims.length-1, 0, 0, true);
      m
    } else a;
		if (Mat.useMKL && inDims(0) == a.dims(0) && outDims(0) == bdims(0)) {             // Use gemm acceleration
		  if (inDims(0) < FFilter.im2colThreshold) {
		    val cellsize = SciFunctions.prod(inDims).v;
		    val outlength = SciFunctions.prod(bdims(1 -> bdims.length)).v;
		  	val i2cmat = FMat.newOrCheckFMat(cellsize, outlength, null, a.GUID, GUID, hmm, "convIm2col".##);
		  	_im2col(apadmat, b, i2cmat, inDims.length-1, cellsize, 0, 0);
		  	if (doclear) b.clear;
		  	sgemm(ORDER.ColMajor, TRANSPOSE.NoTrans, TRANSPOSE.NoTrans, bdims(0), outlength, cellsize, 1f,
		  	    data, bdims(0), i2cmat.data, cellsize, 1f, b.data, bdims(0));
		  } else {
		  	val ddiff = (inDims - 1)/stride - outDims + 1;
		  	ddiff(0) = 0;
		  	val bpadmat = if (ddiff.data.exists(_ != 0)) {
		  		FND.newOrCheckFND(b.dims + ddiff, null, a.GUID, GUID, hmm, "convoutpad".##);
		  	} else {
		  		b;
		  	}
		  	val firststride = 2 + find((stride(2->stride.length) > 1) \ 1)(0); 
		  	if (doclear) { 
		  		bpadmat.clear;
		  	} else {
		  		_copy_padded(b, bpadmat, ddiff/2, inDims.length-1, 0, 0, true);
		  	}
		  	_fast_convolve(apadmat, bpadmat, inDims.length-1, 0, 0, 0, firststride, Filter.forward);
		  	_copy_padded(bpadmat, b, ddiff/2, inDims.length-1, 0, 0, false);
		  }
		} else {
			if (doclear) b.clear; 
			_convolve(apadmat, b, inDims.length-1, 0, 0, 0, Filter.forward);
		}
		Mat.nflops += computeFlops(a, stride, pad);
		b;
	};
  
	def convolve(a:FND):FND = convolve(a, null, true);
  
  def convolveT(b:FND, omat:ND, doclear:Boolean):FND = {
    val bdims = b.dims;
    val adims = Filter.getInputDims(b.dims, inDims, outDims, stride, pad);
    val hmm = Filter.hashIMat(stride, Filter.hashIMat(pad));
    val a = FND.newOrCheckFND(adims, omat, b.GUID, GUID, hmm, "convTin".##);
    val apadmat = if (pad.data.exists(_ != 0)) {
      val m = FND.newOrCheckFND(a.dims + pad * 2, null, b.GUID, GUID, hmm, "convTinpad".##);
      if (!doclear) _copy_padded(a, m, pad, inDims.length-1, 0, 0, true);
      m
    } else a;
    if (doclear) apadmat.clear;
    val ddiff = (inDims - 1)/stride - outDims + 1;
    ddiff(0) = 0;
    if (Mat.useMKL && inDims(0) == adims(0) && outDims(0) == b.dims(0)) {             // Use gemm acceleration
    	if (inDims(0) < FFilter.im2colThreshold) {
    		val cellsize = SciFunctions.prod(inDims).v;
    		val outlength = SciFunctions.prod(bdims(1 -> bdims.length)).v;
    		val i2cmat = FMat.newOrCheckFMat(cellsize, outlength, null, a.GUID, GUID, hmm, "convTIm2col".##);
    		sgemm(ORDER.ColMajor, TRANSPOSE.Trans, TRANSPOSE.NoTrans, cellsize, outlength, bdims(0), 1f,
    				data, bdims(0), b.data, bdims(0), 0f, i2cmat.data, cellsize);
    		_col2im(apadmat, b, i2cmat, inDims.length-1, cellsize, 0, 0);
    	} else {
    		val bpadmat = if (ddiff.data.exists(_ != 0)) {
    			val m = FND.newOrCheckFND(b.dims + ddiff, null, b.GUID, GUID, hmm, "convToutpad".##);
    			m.clear;
    			_copy_padded(b, m, ddiff/2, inDims.length-1, 0, 0, true);
    			m;
    		} else {
    			b;
    		}
    		val firststride = 2 + find((stride(2->stride.length) > 1) \ 1)(0); 
    		_fast_convolve(apadmat, bpadmat, inDims.length-1, 0, 0, 0, firststride, Filter.backwardGradient);
    	} 
    } else {
    	_convolve(apadmat, b, inDims.length-1, 0, 0, 0, Filter.backwardGradient);
    }
    if (apadmat.GUID != a.GUID) {
    	_copy_padded(apadmat, a, pad, inDims.length-1, 0, 0, false);
    }
    Mat.nflops += computeFlops(a, stride, pad);
    a;
  };
  
  def convolveT(a:FND):FND = convolveT(a, null, true);
  
  def convolveM(a:FND, b:FND, doclear:Boolean):FND = {
    val bdims = b.dims;
		val outdims = Filter.getOutputDims(a.dims, inDims, outDims, stride, pad);
		if ((bdims - outdims).data.exists(_ != 0)) {
		  throw new RuntimeException("Output dimensions mismatch in convolveM")
		}
    val hmm = Filter.hashIMat(stride, Filter.hashIMat(pad));
		val apadmat = if (pad.data.exists(_ != 0)) {
      val m = FND.newOrCheckFND(a.dims + pad * 2, null, a.GUID, b.GUID, hmm, "convMinpad".##);
      m.clear;
      _copy_padded(a, m, pad, inDims.length-1, 0, 0, true);
      m
    } else a;
		if (doclear) clear;
		if (Mat.useMKL && inDims(0) == a.dims(0) && outDims(0) == outdims(0)) {             // Use gemm acceleration
			if (inDims(0) < FFilter.im2colThreshold) {
				val cellsize = SciFunctions.prod(inDims).v;
				val outlength = SciFunctions.prod(bdims(1 -> bdims.length)).v;
				val i2cmat = FMat.newOrCheckFMat(cellsize, outlength, null, a.GUID, GUID, hmm, "convMIm2col".##);
				_im2col(apadmat, b, i2cmat, inDims.length-1, cellsize, 0, 0);
				sgemm(ORDER.ColMajor, TRANSPOSE.NoTrans, TRANSPOSE.Trans, bdims(0), cellsize, outlength, 1f,
						b.data, bdims(0), i2cmat.data, cellsize, if (doclear) 0f else 1f, data, bdims(0));
			} else {
				val outdiff = (inDims - 1)/stride - outDims + 1;
				outdiff(0) = 0;
				val bpadmat = if (outdiff.data.exists(_ != 0)) {
					val m = FND.newOrCheckFND(b.dims + outdiff, null, a.GUID, b.GUID, hmm, "convMoutpad".##);
					m.clear;
					_copy_padded(b, m, outdiff/2, inDims.length-1, 0, 0, true);
					m;
				} else {
					b;
				}
				val firststride = 2 + find((stride(2->stride.length) > 1) \ 1)(0); 
				_fast_convolve(apadmat, bpadmat, inDims.length-1, 0, 0, 0, firststride, Filter.backwardModel);
			}
		} else {
			_convolve(apadmat, b, inDims.length-1, 0, 0, 0, Filter.backwardModel);
		}
		Mat.nflops += computeFlops(a, stride, pad);
		this;
	};
	
  def convolveM(a:FND, b:FND):FND = convolveM(a, b, true);

	def _fast_convolve(a:FND, b:FND, idim:Int, astart:Int, bstart:Int, fstart:Int, firststride:Int, convType:Int) {
		val adims = a.dims;
		val bdims = b.dims;
		val iwidth = inDims(idim);
		val owidth = outDims(idim);
		val kstride = stride(idim);
		if (idim > 0) {
			var astep = 1;
			var bstep = 1;
			var fstep = 1;
			var ix = 0;
			while (ix < idim) {
				astep *= adims(ix);
				bstep *= bdims(ix);
				fstep *= inDims(ix);
				fstep *= outDims(ix);
				ix += 1; 
			}
			var j = 0;
			while (j < owidth) {
				var i = 0;
				while (i < iwidth) {
				  if (idim >= firststride) {
				  	var k = pad(idim);
				  	var ks = 0;
				  	while (ks + iwidth - 1 < adims(idim)) {
				  		_fast_convolve(a, b, idim-1, 
				  				astart + (ks + i) * astep, 
				  				bstart + (k + j) * bstep,
				  				fstart + (j + i * owidth) * fstep,
				  				firststride, convType);
				  		k += 1;
				  		ks += kstride;
				  	} 
				  } else {
				  	_fast_convolve(a, b, idim-1, 
				  			astart + i * astep, 
				  			bstart + (pad(idim) + j) * bstep,
				  			fstart + (j + i * owidth) * fstep,
				  			firststride, convType);
				  }
					i += 1;
				}
				j += 1;
			}
		} else {
		  var mstep = 1;
			var ix = 1;
			while (ix < firststride) {
				mstep *= bdims(ix);
				ix += 1; 
			}
      val bstart0 = bstart + pad(0);
      val mstep0 = math.min(mstep, math.min((a.length - astart)/(iwidth*stride(1)), (b.length - bstart0)/owidth));
			convType match {			  
			  case Filter.forward => {
			  	sgemmx(ORDER.ColMajor, TRANSPOSE.NoTrans, TRANSPOSE.NoTrans, owidth, mstep0, iwidth, 1f,
			  			data, fstart, owidth, 
			  			a.data, astart, iwidth*stride(1), 1f,
			  			b.data, bstart0, owidth);	
			  }
			  case Filter.backwardGradient => {
			  	sgemmx(ORDER.ColMajor, TRANSPOSE.Trans, TRANSPOSE.NoTrans, iwidth, mstep0, owidth, 1f,
			  			data, fstart, owidth, 
			  			b.data, bstart0, owidth, 1f,
			  			a.data, astart, iwidth*stride(1));	
			  }
			  case Filter.backwardModel => {
			  	sgemmx(ORDER.ColMajor, TRANSPOSE.NoTrans, TRANSPOSE.Trans, owidth, iwidth, mstep0, 1f, 
			  			b.data, bstart0, owidth,  
			  			a.data, astart, iwidth*stride(1), 1f,
			  			data, fstart, owidth);	
			  }
			}
		}
	}
	
	// copy into padded array, or from padded array

	def _copy_padded(in:FND, out:FND, padx:IMat, idim:Int, astart:Int, bstart:Int, topadded:Boolean) {
	  val idims = in.dims;
	  val odims = out.dims;
	  val width = if (topadded) idims(idim) else odims(idim);
	  val ipad = padx(idim);
	  if (idim > 1) {
	  	var instep = 1;
	  	var outstep = 1;
	  	var ix = 0;
	  	while (ix < idim) {
	  		instep *= idims(ix);
	  		outstep *= odims(ix);
	  		ix += 1; 
	  	}
	  	var i = 0;
	  	while (i < width) {
        if (topadded) {
          _copy_padded(in, out, padx, idim-1, astart + instep * i, bstart + outstep * (ipad + i), topadded);
        } else {
        	_copy_padded(in, out, padx, idim-1, astart + instep * (ipad + i), bstart + outstep * i, topadded);
        }
	  		i += 1;
	  	}
	  } else if (idim == 1){
	    val instep = idims(0);
	    val outstep = odims(0);
	    val ipad0 = padx(0);
	    val width0 = if (topadded) idims(0) else odims(0);
	    var i = 0;
	    if (topadded) {
	    	while (i < width) {
	    		val astart0 = astart + instep * i;
	    		val bstart0 = bstart + outstep * (ipad + i);
	    		var j = 0;
	    		while (j < width0) {
	    			out.data(bstart0 + ipad0 + j) = in.data(astart0 + j);
	    			j += 1;
	    		}
	    		i += 1;
        }
	    } else {
	      while (i < width) {
	    		val astart0 = astart + instep * (ipad + i);
	    		val bstart0 = bstart + outstep * i;
	    		var j = 0;
	    		while (j < width0) {
	    			out.data(bstart0 + j) = in.data(astart0 + ipad0 + j);
	    			j += 1;
	    		}
	    		i += 1;
	      }
	    }
	  } else {
	  	var i = 0;
      if (topadded) {
        while (i < width) {
          out.data(bstart + ipad + i) = in.data(astart + i);
          i += 1;
        }       
      } else {
    	  while (i < width) {
    		  out.data(bstart + i) = in.data(astart + ipad + i);
    		  i += 1;
    	  }
      }
	  }
	}

	def _convolve(a:FND, b:FND, idim:Int, astart:Int, bstart:Int, fstart:Int, convType:Int) {
		val adims = a.dims;
		val bdims = b.dims;
		val iwidth = inDims(idim);
		val owidth = outDims(idim);
		val kstride = stride(idim);
		if (idim > 0) {
			var astep = 1;
			var bstep = 1;
			var fstep = 1;
			var ix = 0;
			while (ix < idim) {
				astep *= adims(ix);
				bstep *= bdims(ix);
				fstep *= inDims(ix);
				fstep *= outDims(ix);
				ix += 1; 
			}
			var k = 0;
			var ks = 0;
			while (ks + iwidth - 1 < adims(idim)) {
				var j = 0;
				while (j < owidth) {
					var i = 0;
					while (i < iwidth) {
						_convolve(a, b, idim-1, 
								astart + (ks + i) * astep, 
								bstart + (k + j) * bstep,
								fstart + (j + i * owidth) * fstep,
								convType);
						i += 1;
					}
					j += 1;
				}
				k += 1;
				ks += kstride;
			}
		} else {
			var k = 0;
			var ks = 0;
			while (ks + iwidth - 1 < adims(0)) {              // Move forward over input+output tensors              
				convType match {
				  case Filter.forward => {
				  	var i = 0;
				  	var ifilt = 0;
				  	while (i < iwidth) {                      // Move over input tensor
				  	  val adat = a.data(astart + ks + i);
				  		var j = 0;
				  		while (j < owidth) {                        // Move over output tensor
				  			b.data(bstart + k + j) += data0(fstart + ifilt + j) * adat;
				  			j += 1;
				  		}
				  		i += 1;
				  		ifilt += owidth;
				  	}
				  }
				  case Filter.backwardGradient => {
				  	var i = 0;
				  	var ifilt = 0;
				  	while (i < iwidth) {
				  	  var ss = 0f;
				  		var j = 0;
				  		while (j < owidth) {
				  			ss += b.data(bstart + k + j) * data0(fstart + j + ifilt);
				  			j += 1;
				  		}
				  		a.data(astart + ks + i) += ss;
				  		i += 1;
				  		ifilt += owidth;
				  	}
				  }
				  case Filter.backwardModel => {
				  	var i = 0;
				  	var ifilt = 0;
				  	while (i < iwidth) {
				  	  val adat = a.data(astart + ks + i);
				  		var j = 0;
				  		while (j < owidth) {
				  			data0(fstart + j + ifilt) += b.data(bstart + k + j) * adat;
				  			j += 1;
				  		}
				  		i += 1;
				  		ifilt += owidth;
				  	}
				  }
				}
				k += 1;
				ks += kstride;
			}
		}
	};
	
	def _im2col(a:FND, b:FND, i2c:FMat, idim:Int, celldim:Int, astart:Int, bstart:Int):FMat = {
	  val adims = a.dims;
	  val bdims = b.dims;
		val iwidth = inDims(idim);
		val owidth = outDims(idim);
	  if (idim > 1) {
	  	var astep = 1;
	  	var bstep = 1;
	  	var cellstep = 1;
	  	var ix = 0;
	  	while (ix < idim) {
	  		astep *= adims(ix);
	  		if (ix > 0) bstep *= bdims(ix);
	  		cellstep *= inDims(ix);
	  		ix += 1; 
	  	}
	  	var iin = 0;
	  	var iout = 0;
	  	while (iout < bdims(idim)) {
	  	  var j = 0;
	  	  while (j < iwidth) {
	  	  	_im2col(a, b, i2c, idim - 1, celldim, astart + astep * (iin + j), bstart + celldim * bstep * iout + cellstep * j);
	  	  	j += 1;
	  	  }
	  	  iout += 1;
	  		iin += stride(idim);
	  	}
	  } else if (idim == 1) {
	  	var	astep = adims(0);
	  	var	cellstep = inDims(0);
	  	var iin = 0;
	  	var iout = 0;
	  	while (iout < bdims(idim)) {
	  	  var j = 0;
	  	  while (j < iwidth) {
	  	  	val astart0 = astart + astep * (iin + j);
	  	  	val bstart0 = bstart + celldim * iout + cellstep * j;
	  	  	if (cellstep > 16) {
	  	  		System.arraycopy(a.data, astart0, i2c.data, bstart0, cellstep);
	  	  	} else {
	  	  		var k = 0;
	  	  		while (k < cellstep) {
	  	  			i2c.data(bstart0 + k) = a.data(astart0 + k);
	  	  			k += 1;
	  	  		}
	  	    }
	  	  	j += 1;
	  	  }
	  	  iout += 1;
	  		iin += stride(idim);
	  	}	    
	  } else {
	    if (iwidth > 16) {
	      System.arraycopy(a.data, astart, i2c.data, bstart, iwidth);
	    } else {
	      var i = 0;
	      while (i < iwidth) {
	        i2c.data(bstart + i) = a.data(astart + i);
	        i += 1;
	      }
	    }
	  }
	  i2c
	}
	
	
		
	def _col2im(a:FND, b:FND, i2c:FMat, idim:Int, celldim:Int, astart:Int, bstart:Int):FMat = {
	  val adims = a.dims;
	  val bdims = b.dims;
		val iwidth = inDims(idim);
		val owidth = outDims(idim);
	  if (idim > 1) {
	  	var astep = 1;
	  	var bstep = 1;
	  	var cellstep = 1;
	  	var ix = 0;
	  	while (ix < idim) {
	  		astep *= adims(ix);
	  		if (ix > 0) bstep *= bdims(ix);
	  		cellstep *= inDims(ix);
	  		ix += 1; 
	  	}
	  	var iin = 0;
	  	var iout = 0;
	  	while (iout < bdims(idim)) {
	  	  var j = 0;
	  	  while (j < iwidth) {
	  	  	_col2im(a, b, i2c, idim - 1, celldim, astart + astep * (iin + j), bstart + celldim * bstep * iout + cellstep * j);
	  	  	j += 1;
	  	  }
	  	  iout += 1;
	  		iin += stride(idim);
	  	}
	  } else if (idim == 1) {
	  	var	instep = adims(0);
	  	var	cellstep = inDims(0);
	  	var iin = 0;
	  	var iout = 0;
	  	while (iout < bdims(idim)) {
	  	  var j = 0;
	  	  while (j < iwidth) {
	  	  	val astart0 = astart + instep * (iin + j);
	  	  	val bstart0 = bstart + celldim * iout + cellstep * j;
	  	  	var k = 0;
	  	  	while (k < cellstep) {
	  	  		a.data(astart0 + k) += i2c.data(bstart0 + k);
	  	  		k += 1;
	  	    }
	  	  	j += 1;
	  	  }
	  	  iout += 1;
	  		iin += stride(idim);
	  	}	    
	  } else {
	  	var i = 0;
	  	while (i < iwidth) {
	  		a.data(astart + i) += i2c.data(bstart + i);
	  		i += 1;
	  	}
	  }
	  i2c;
	}
	
  override def * (a:FND):FND = {
    convolve(a);
  }
  
  def ^* (a:FND):FND = {
    convolveT(a);
  }
  
  override def * (a:ND):ND = {
    a match {
      case aa:FND => convolve(aa);
    }
  };
  
  override def ^* (a:ND):ND = {
  	a match {
  	case aa:FND => convolveT(aa);
  	}
  }
}


class FFiltPair(val omat:Filter, val a:FND) extends Pair {
  def *^ (b:FND):FND = {
    omat.asInstanceOf[FFilter].convolveM(a, b);
  }
}

object FFilter {
  
  var im2colThreshold = 16;
  
	def FFilter1D(w:Int, nstride:Int, npad:Int) = {
		val inDims = irow(w);
		val outDims = irow(1);
		val stride = irow(nstride);
		val pad = irow(npad);
		new FFilter(inDims, outDims, stride, pad, new Array[Float](w));
	}

	def FFilter1Dd(w:Int, din:Int, dout:Int, nstride:Int, npad:Int) = {
		val inDims = irow(din, w);
		val outDims = irow(dout, 1);
		val stride = irow(1, nstride);
		val pad = irow(0, npad);
		new FFilter(inDims, outDims, stride, pad, new Array[Float](din*dout*w));
	}

	def FFilter2D(w:Int, h:Int, nstride:Int, npad:Int) = {
		val inDims = irow(w, h);
		val outDims = irow(1, 1);
		val stride = irow(nstride, nstride);
		val pad = irow(npad, npad);
		new FFilter(inDims, outDims, stride, pad, new Array[Float](w*h));
	}

	def FFilter2Dd(w:Int, h:Int, din:Int, dout:Int, nstride:Int, npad:Int) = {
		val inDims = irow(din, w, h);
		val outDims = irow(dout, 1, 1);
		val stride = irow(1, nstride, nstride);
		val pad = irow(0, npad, npad);
		new FFilter(inDims, outDims, stride, pad, new Array[Float](din*dout*w*h));
	}
  
  def FFilter2Ddn(w:Int, h:Int, din:Int, dout:Int, nstride:Int, npad:Int) = {
    val inDims = irow(din, w, h, 1);
    val outDims = irow(dout, 1, 1, 1);
    val stride = irow(1, nstride, nstride, 1);
    val pad = irow(0, npad, npad, 0);
    new FFilter(inDims, outDims, stride, pad, new Array[Float](din*dout*w*h));
  }

}
