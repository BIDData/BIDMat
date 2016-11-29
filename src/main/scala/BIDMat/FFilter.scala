package BIDMat;
import BIDMat.MatFunctions._;
import edu.berkeley.bid.CBLAS._
import scala.util.hashing.MurmurHash3

@SerialVersionUID(100L)
class FFilter(inDims0:IMat, outDims0:IMat, stride0:IMat, pad0:IMat, data0:Array[Float]) extends
FND((inDims0 *@ outDims0).data, data0) with Filter {

	override val inDims = inDims0;
	override val outDims = outDims0;
	override val stride = if (stride0.asInstanceOf[AnyRef] != null) stride0 else iones(1, inDims.length);
	override val pad = if (pad0.asInstanceOf[AnyRef] != null) pad0 else izeros(1,inDims.length);
  var outpadmat:FND = null;

	def convolve(a:FND, omat:ND, doclear:Boolean):FND = {
		val outdims = Filter.getOutputDims(a.dims, inDims, outDims, stride, pad);
    val hmm = hashIMat(stride, hashIMat(pad));
		val out = FND.newOrCheckFND(outdims, omat, a.GUID, GUID, hmm, "convout".##);
		val inpadmat = if (SciFunctions.sum(pad).dv > 0) pad(a, pad) else a;
		if (Mat.useMKL && inDims(0) == a.dims(0) && outDims(0) == outdims(0)) {             // Use gemm acceleration
			val outpaddims = if (SciFunctions.sum(pad).dv > 0) out.dims + pad*2 else out.dims;
			outpadmat = FND.newOrCheckFND(outpaddims, null, a.GUID, GUID, hmm, "convoutpad".##);
			val firststride = 2 + find((stride(2->stride.length) > 1) \ 1)(0); 
			if (doclear) { 
				outpadmat.clear;
      } else {
    	  _copy_padded(out, outpadmat, inDims.length-1, 0, 0, true);
      }
			_fast_convolve(inpadmat, outpadmat, inDims.length-1, 0, 0, 0, firststride, Filter.forward);
			_copy_padded(outpadmat, out, inDims.length-1, 0, 0, false);
		} else {
			if (doclear) { 
        out.clear;
      } 
			_convolve(inpadmat, out, inDims.length-1, 0, 0, 0, Filter.forward);
		}
		Mat.nflops += computeFlops(a, stride, pad);
		out;
	};
  
	def convolve(a:FND):FND = convolve(a, null, true);
  
  def hashIMat(a:IMat, start:Int):Int = {
    var i = 0; 
    var hv = start;
    while (i < a.length) {
      hv = MurmurHash3.mix(hv, a.data(i));
      i += 1;
    }
    hv;
  }
  
  def hashIMat(a:IMat):Int = hashIMat(a, 23412154);

	def _fast_convolve(in:FND, out:FND, idim:Int, astart:Int, bstart:Int, fstart:Int, firststride:Int, convType:Int) {
		val idims = in.dims;
		val odims = out.dims;
		val iwidth = inDims(idim);
		val owidth = outDims(idim);
		val kstride = stride(idim);
		if (idim > 0) {
			var instep = 1;
			var outstep = 1;
			var fstep = 1;
			var ix = 0;
			while (ix < idim) {
				instep *= idims(ix);
				outstep *= odims(ix);
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
				  	while (ks + iwidth - 1 < idims(idim)) {
				  		_fast_convolve(in, out, idim-1, 
				  				astart + (ks + iwidth - i - 1) * instep, 
				  				bstart + (k + owidth - j - 1) * outstep,
				  				fstart + (i + j * iwidth) * fstep,
				  				firststride, convType);
				  		k += 1;
				  		ks += kstride;
				  	} 
				  } else {
				  	_fast_convolve(in, out, idim-1, 
				  			astart + (iwidth - i - 1) * instep, 
				  			bstart + (pad(idim) + owidth - j - 1) * outstep,
				  			fstart + (i + j * iwidth) * fstep,
				  			firststride, convType);
				  }
					i += 1;
				}
				j += 1;
			}
		} else {
		  var outstep = 1;
			var ix = 1;
			while (ix < firststride) {
				outstep *= odims(ix);
				ix += 1; 
			}
      val bstart0 = bstart + pad(0);
      val outstep0 = outstep - (owidth + bstart0 - 1)/owidth;
//      println("conv dims %d %d %d" format (owidth, outstep0, iwidth))
			convType match {			  
			  case Filter.forward => {
			  	sgemmx(ORDER.ColMajor, TRANSPOSE.Trans, TRANSPOSE.NoTrans, owidth, outstep0, iwidth, 1f,
			  			data, fstart, iwidth, 
			  			in.data, astart, iwidth*stride(1), 1f,
			  			out.data, bstart0, owidth);	
			  }
			  case Filter.backwardGradient => {
			  	sgemmx(ORDER.ColMajor, TRANSPOSE.NoTrans, TRANSPOSE.NoTrans, iwidth, outstep0, owidth, 1f,
			  			data, fstart, iwidth, 
			  			out.data, bstart0, owidth, 1f,
			  			in.data, astart, iwidth*stride(1));	
			  }
			  case Filter.backwardModel => {
			  	sgemmx(ORDER.ColMajor, TRANSPOSE.NoTrans, TRANSPOSE.Trans, iwidth, owidth, outstep0, 1f, 
			  			in.data, astart, iwidth*stride(1), 
              out.data, bstart0, owidth, 1f,
			  			data, fstart, iwidth);	
			  }
			}
		}
	}

	def _copy_padded(in:FND, out:FND, idim:Int, astart:Int, bstart:Int, dopad:Boolean) {
	  val idims = in.dims;
	  val odims = out.dims;
	  val owidth = odims(idim);
	  val ipad = pad(idim);
	  if (idim > 0) {
	  	var instep = 1;
	  	var outstep = 1;
	  	var ix = 0;
	  	while (ix < idim) {
	  		instep *= idims(ix);
	  		outstep *= odims(ix);
	  		ix += 1; 
	  	}
	  	var i = 0;
	  	while (i < owidth) {
        if (dopad) {
          _copy_padded(in, out, idim-1, astart + instep * i, bstart + outstep * (ipad + i), dopad);
        } else {
        	_copy_padded(in, out, idim-1, astart + instep * (ipad + i), bstart + outstep * i, dopad);
        }
	  		i += 1;
	  	}
	  } else {
	  	var i = 0;
      if (dopad) {
        while (i < owidth) {
          out.data(bstart + ipad + i) = in.data(astart + i);
          i += 1;
        }       
      } else {
    	  while (i < owidth) {
    		  out.data(bstart + i) = in.data(astart + ipad + i);
    		  i += 1;
    	  }
      }
	  }
	}

	def _convolve(in:FND, out:FND, idim:Int, astart:Int, bstart:Int, fstart:Int, convType:Int) {
		val idims = in.dims;
		val odims = out.dims;
		val iwidth = inDims(idim);
		val owidth = outDims(idim);
		val kstride = stride(idim);
		if (idim > 0) {
			var instep = 1;
			var outstep = 1;
			var fstep = 1;
			var ix = 0;
			while (ix < idim) {
				instep *= idims(ix);
				outstep *= odims(ix);
				fstep *= inDims(ix);
				fstep *= outDims(ix);
				ix += 1; 
			}
			var k = 0;
			var ks = 0;
			while (ks + iwidth - 1 < idims(idim)) {
				var j = 0;
				while (j < owidth) {
					var i = 0;
					while (i < iwidth) {
						_convolve(in, out, idim-1, 
								astart + (ks + iwidth - i - 1) * instep, 
								bstart + (k + owidth - j - 1) * outstep,
								fstart + (i + j * iwidth) * fstep,
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
			while (ks + iwidth - 1 < idims(0)) {              // Move forward over input+output tensors              
				convType match {
				  case Filter.forward => {
				  	var j = 0;
				  	while (j < owidth) {                        // Move over output tensor
				  		val jfwidth = j * iwidth;
				  		var ss = 0f;
				  		var i = 0;
				  		while (i < iwidth) {                      // Move over input tensor
				  			ss += in.data(astart + ks + i) * data0(fstart + jfwidth + i);
				  			i += 1;
				  		}
				  		out.data(bstart + k + j) += ss;
				  		j += 1;
				  	}
				  }
				  case Filter.backwardGradient => {
				  	var j = 0;
				  	while (j < owidth) {
				  		val jfwidth = j * iwidth;
				  		val odata = out.data(bstart + k + owidth - j - 1);
				  		var i = 0;
				  		while (i < iwidth) {
				  			in.data(astart + ks + iwidth - i - 1) += odata * data0(fstart + jfwidth + i);
				  			i += 1;
				  		}
				  		j += 1;
				  	}
				  }
				  case Filter.backwardModel => {
				  	var j = 0;
				  	while (j < owidth) {
				  	  val jfwidth = j * iwidth;
				  	  val odata = out.data(bstart + k + owidth - j - 1);
				  		var i = 0;
				  		while (i < iwidth) {
				  		  data0(fstart + jfwidth + i) += odata * in.data(astart + ks + iwidth - i - 1);
				  			i += 1;
				  		}
				  		j += 1;
				  	}
				  }
				}
				k += 1;
				ks += kstride;
			}
		}
	};

	def correlate(a:FND, omat:ND, doclear:Boolean):FND = {
			val outdims = Filter.getOutputDims(a.dims, inDims, outDims, stride, pad);
			val hmm = hashIMat(stride, hashIMat(pad));
			val out = FND.newOrCheckFND(outdims, omat, a.GUID, GUID, hmm, "convout".##);
			val inpadmat = if (SciFunctions.sum(pad).dv > 0) pad(a, pad) else a;
			if (Mat.useMKL && inDims(0) == a.dims(0) && outDims(0) == outdims(0)) {             // Use gemm acceleration
				val outpaddims = if (SciFunctions.sum(pad).dv > 0) out.dims + pad*2 else out.dims;
				val outpadmat = FND.newOrCheckFND(outpaddims, null, a.GUID, GUID, hmm, "convoutpad".##);
				val firststride = 1 + find((stride(1->stride.length) > 1) \ 1)(0); 
        if (doclear) {
        	outpadmat.clear;          
        } else {
        	_copy_padded(out, outpadmat, inDims.length-1, 0, 0, true);
        }
				_fast_correlate(inpadmat, outpadmat, inDims.length-1, 0, 0, 0, firststride, Filter.forward);
				_copy_padded(outpadmat, out, inDims.length-1, 0, 0, false);
			} else {
				if (doclear) {
					out.clear;
        }
				_correlate(inpadmat, out, inDims.length-1, 0, 0, 0, Filter.forward);
			}
			Mat.nflops += computeFlops(a, stride, pad);
			out
	};
  
	def correlate(a:FND):FND = correlate(a, null, true);
	
	def _fast_correlate(in:FND, out:FND, idim:Int, astart:Int, bstart:Int, fstart:Int, firststride:Int, convType:Int) {
		val idims = in.dims;
		val odims = out.dims;
		val iwidth = inDims(idim);
		val owidth = outDims(idim);
		val kstride = stride(idim);
		if (idim > 0) {
			var instep = 1;
			var outstep = 1;
			var fstep = 1;
			var ix = 0;
			while (ix < idim) {
				instep *= idims(ix);
				outstep *= odims(ix);
				fstep *= inDims(ix);
				fstep *= outDims(ix);
				ix += 1; 
			}
			val aiwidth = if (idim == 1) 1 else iwidth;
			var j = 0;
			while (j < owidth) {
				var i = 0;
				while (i < aiwidth) {
				  if (idim > firststride) {
				  	var k = 0;
				  	var ks = 0;
				  	while (ks + iwidth - 1 < idims(idim)) {
				  		_fast_correlate(in, out, idim-1, 
				  				astart + (ks + i) * instep, 
				  				bstart + (k + j) * outstep,
				  				fstart + (i + j * iwidth) * fstep,
				  				firststride, convType);
				  		k += 1;
				  		ks += kstride;
				  	} 
				  } else {
				  	_fast_correlate(in, out, idim-1, 
				  			astart + i * instep, 
				  			bstart + j * outstep,
				  			fstart + (i + j * iwidth) * fstep,
				  			firststride, convType);
				  }
					i += 1;
				}
				j += 1;
			}
		} else {
		  var outstep = 1;
			var ix = 1;
			while (ix < firststride) {
				outstep *= odims(ix);
				ix += 1; 
			}
      val bstart0 = bstart + pad(0);
      val outstep0 = outstep - (owidth + bstart0 - 1)/owidth;
			convType match {
			  case Filter.forward => {
			  	sgemmx(ORDER.ColMajor, TRANSPOSE.Trans, TRANSPOSE.NoTrans, owidth, outstep0/inDims(1), iwidth*inDims(1), 1f,
			  			data, fstart, iwidth, 
			  			in.data, astart, iwidth*inDims(1), 1f,
			  			out.data, bstart0, owidth);	
			  }
			  case Filter.backwardGradient => {
			  	sgemmx(ORDER.ColMajor, TRANSPOSE.NoTrans, TRANSPOSE.NoTrans, iwidth*inDims(1), outstep0/inDims(1), owidth, 1f,
			  			data, fstart, iwidth,
			  			out.data, bstart0, owidth, 1f,
			  			in.data, astart, iwidth*inDims(1));	
			  }
			  case Filter.backwardModel => {
			  	sgemmx(ORDER.ColMajor, TRANSPOSE.NoTrans, TRANSPOSE.Trans, iwidth*inDims(1), owidth, outstep0/inDims(1), 1f,		  			 
			  			in.data, astart, iwidth*inDims(1), 
              out.data, bstart0, owidth, 1f,
			  			data, fstart, iwidth);	
			  }
			}
		}
	}

	def _correlate(in:FND, out:FND, idim:Int, astart:Int, bstart:Int, fstart:Int, convType:Int) {
		val idims = in.dims;
		val odims = out.dims;
		val iwidth = inDims(idim);
		val owidth = outDims(idim);
		val kstride = stride(idim);
		if (idim > 0) {
			var instep = 1;
			var outstep = 1;
			var fstep = 1;
			var ix = 0;
			while (ix < idim) {
				instep *= idims(ix);
				outstep *= odims(ix);
				fstep *= inDims(ix);
				fstep *= outDims(ix);
				ix += 1; 
			}
			var k = 0;
			var ks = 0;
			while (ks + iwidth - 1 < idims(idim)) {
				var j = 0;
				while (j < owidth) {
					var i = 0;
					while (i < iwidth) {
						_correlate(in, out, idim-1, 
								astart + (ks + i) * instep, 
								bstart + (k + j) * outstep,
								fstart + (i + j * iwidth) * fstep,
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
			while (ks + iwidth - 1 < idims(0)) {           // Move forward over input+output tensors  
			  convType match {
			    case Filter.forward => {
			    	var j = 0;
			    	while (j < owidth) {                         // Move over output tensor
			    		var ss = 0f;
			    		val jfwidth = j * iwidth;
			    		var i = 0;
			    		while (i < iwidth) {                       // Move over input tensor
			    			ss += in.data(astart + ks + i) * data0(fstart + jfwidth + i);
			    			i += 1;
			    		}
			    		out.data(bstart + k + j) += ss;
			    		j += 1;
			    	}
			    }
			    case Filter.backwardGradient => {
			    	var j = 0;
			    	while (j < owidth) { 		
			    		val jfwidth = j * iwidth;
			    		val odata = out.data(bstart + k + j);
			    		var i = 0;
			    		while (i < iwidth) {
			    			in.data(astart + ks + i) += odata * data0(fstart + jfwidth + i);
			    			i += 1;
			    		}
			    		j += 1;
			    	}
			    }
			    case Filter.backwardModel => {
			    	var j = 0;
			    	while (j < owidth) {                         // Move over output tensor
			    	  val jfwidth = j * iwidth;
			    	  val odata = out.data(bstart + k + j);
			    		var i = 0;
			    		while (i < iwidth) {
			    			data0(fstart + jfwidth + i) += odata* in.data(astart + ks + i);
			    			i += 1;
			    		}
			    		j += 1;
			    	}
			    }
			  }
				k += 1;
				ks += kstride;
			}
		}
	};


	def pad(a:FND, pad:IMat):FND = { 
			val inpaddims = if (pad.asInstanceOf[AnyRef] != null) a.dims + pad*2 else a.dims;
			val out = FND.newOrCheckFND(inpaddims, null, a.GUID, GUID, "convpad".##);
			_pad(a, pad, out, inpaddims.length-1, 0, 0);
			out
	};

	def _zerofill(out:FND, idim:Int, bstart0:Int) {
		val odims = out.dims;
		var bstart = bstart0;
		var i = 0;
		if (idim > 0) {
			var outcol = 1;
			var ix = 0;
			while (ix < idim) {
				outcol *= odims(ix);
				ix += 1; 
			}    
			while (i < odims(idim)) {
				_zerofill(out, idim-1, bstart);
				bstart += outcol;
				i += 1;
			}
		} else {
			while (i < odims(idim)) {
				out.data(bstart) = 0;
				bstart += 1;
				i += 1;
			}
		}
	};

	def _pad(a:FND, pad:IMat, out:FND, idim:Int, astart0:Int, bstart0:Int) {
		val adims = a.dims;
		val odims = out.dims;
		var astart = astart0;
		var bstart = bstart0;
		var i = 0;
		if (idim > 0) {
			var incol = 1;
			var outcol = 1;
			var ix = 0;
			while (ix < idim) {
				incol *= adims(ix);
				outcol *= odims(ix);
				ix += 1; 
			}
			while (i < pad(idim)) {
				_zerofill(out, idim-1, bstart);
				bstart += outcol;
				i += 1;
			}
			while (i - pad(idim) < adims(idim)) {
				_pad(a, pad, out, idim-1, astart, bstart);
				astart += incol;
				bstart += outcol;
				i += 1;
			}
			while (i - pad(idim) - adims(idim) < pad(idim)) {
				_zerofill(out, idim-1, bstart);
				bstart += outcol;
				i += 1;
			}
		} else {
			while (i < pad(idim)) {
				out.data(bstart) = 0;
				bstart += 1;
				i += 1;
			}
			while (i - pad(idim) < adims(idim)) {
				out.data(bstart) = a.data(astart);
				astart += 1;
				bstart += 1;
				i += 1;
			}
			while (i - pad(idim) - adims(idim) < pad(idim)) {
				out.data(bstart) = 0;
				bstart += 1;
				i += 1;
			}   
		}
	}
    
  def computeFlops(in:ND, stride:IMat, pad:IMat):Long = {
    var i = 0;
    var flops = 2L;
    while (i < stride.length) {
      flops *= 1L * inDims(i) * outDims(i) * (in.dims(i) - inDims(i) + 1 + 2*pad(i))/stride(i)
      i += 1;
    }
    flops;
  }
}

object FFilter {
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
