package BIDMat;
import BIDMat.MatFunctions._;
import edu.berkeley.bid.CBLAS._

@SerialVersionUID(100L)
class FFilter(inDims0:IMat, outDims0:IMat, stride0:IMat, pad0:IMat, data0:Array[Float]) extends
FND((inDims0 *@ outDims0).data, data0) with Filter {

	override val inDims = inDims0;
	override val outDims = outDims0;
	override val pad = pad0;
	override val stride = stride0;

	def convolve(a:FND):FND = {
		val outdims = Filter.getOutputDims(a.dims, inDims, outDims, stride, pad);
		val out = FND.newOrCheckFND(outdims, null, a.GUID, GUID, "convout".##);
		val inpadmat = if (pad.asInstanceOf[AnyRef] != null && SciFunctions.sum(pad).dv > 0) pad(a, pad) else a;
		out.clear;
		if (inDims(0) == a.dims(0) && outDims(0) == outdims(0)) {             // Use gemm acceleration
			val outpaddims = if (pad.asInstanceOf[AnyRef] != null) out.dims + pad*2 else out.dims;
			val outpadmat = FND.newOrCheckFND(outpaddims, null, a.GUID, GUID, "convoutpad".##);
			val firststride = find((stride(1->stride.length) > 1) \ 1)(0); 
			outpadmat.clear;
			_fast_convolve(inpadmat, outpadmat, inDims.length-1, 0, 0, 0, firststride);
			_copy_padded(outpadmat, out, inDims.length-1, 0, 0);
		} else {
			_convolve(inpadmat, out, inDims.length-1, 0, 0, 0);
		}
		out;
	};

	def _fast_convolve(in:FND, out:FND, idim:Int, astart:Int, bstart:Int, fstart:Int, firststride:Int) {
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
				  if (idim > firststride) {
				  	var k = 0;
				  	var ks = 0;
				  	while (ks + iwidth - 1 < idims(idim)) {
				  		_fast_convolve(in, out, idim-1, 
				  				astart + (ks + iwidth - i - 1) * instep, 
				  				bstart + (k + owidth - j - 1) * outstep,
				  				fstart + (i + j * iwidth) * fstep,
				  				firststride);
				  		k += 1;
				  		ks += kstride;
				  	} 
				  } else {
				  	_fast_convolve(in, out, idim-1, 
				  			astart + (iwidth - i - 1) * instep, 
				  			bstart + (owidth - j - 1) * outstep,
				  			fstart + (i + j * iwidth) * fstep,
				  			firststride);
				  }
					i += 1;
				}
				j += 1;
			}
		} else {
		  var instep = 1;
			var ix = 1;
			while (ix < firststride) {
				instep *= idims(ix);
				ix += 1; 
			}
		  sgemmx(ORDER.ColMajor, TRANSPOSE.NoTrans, TRANSPOSE.NoTrans, outDims(0), instep, inDims(0), 1f,
		      data, fstart, outDims(0), 
		      in.data, astart, outDims(0), 1f,
		      out.data, bstart, outDims(0));	  
		}
	}

	def _copy_padded(in:FND, out:FND, idim:Int, astart:Int, bstart:Int) {
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
	  		_copy_padded(in, out, idim-1, astart + instep * (ipad + i), bstart + outstep * i);
	  		i += 1;
	  	}
	  } else {
	  	var i = 0;
	  	while (i < owidth) {
	  	  out.data(bstart + i) = in.data(astart + ipad + i);
	  	  i += 1;
	  	}
	  }
	}

	def _convolve(in:FND, out:FND, idim:Int, astart:Int, bstart:Int, fstart:Int) {
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
								fstart + (i + j * iwidth) * fstep);
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
				var j = 0;
				while (j < owidth) {                    // Move over output tensor
					val jfwidth = j * iwidth;
					var ss = 0f;
					var i = 0;
					while (i < iwidth) {                      // Move over input tensor
						ss += in.data(astart + ks + iwidth - i - 1) * data0(fstart + jfwidth + i);
						i += 1;
					}
					out.data(bstart + k + outDims(0) - j - 1) += ss;
					j += 1;
				}
				k += 1;
				ks += kstride;
			}
		}
	};

	def correlate(a:FND):FND = {
			val outdims = Filter.getOutputDims(a.dims, inDims, outDims, stride, pad);
			val out = FND.newOrCheckFND(outdims, null, a.GUID, GUID, "convout".##);
			val inpadmat = if (pad.asInstanceOf[AnyRef] != null && SciFunctions.sum(pad).dv > 0) pad(a, pad) else a;
			if (inDims(0) == a.dims(0) && outDims(0) == outdims(0)) {             // Use gemm acceleration
				val outpaddims = if (pad.asInstanceOf[AnyRef] != null) out.dims + pad*2 else out.dims;
				val outpadmat = FND.newOrCheckFND(outpaddims, null, a.GUID, GUID, "convoutpad".##);
				val firststride = find((stride(1->stride.length) > 1) \ 1)(0); 
				outpadmat.clear;
				_fast_correlate(inpadmat, outpadmat, inDims.length-1, 0, 0, 0, firststride);
				_copy_padded(outpadmat, out, inDims.length-1, 0, 0);
			} else {
				out.clear;
				_correlate(inpadmat, out, inDims.length-1, 0, 0, 0);
			}
			out
	};
	
	def _fast_correlate(in:FND, out:FND, idim:Int, astart:Int, bstart:Int, fstart:Int, firststride:Int) {
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
				  				firststride);
				  		k += 1;
				  		ks += kstride;
				  	} 
				  } else {
				  	_fast_correlate(in, out, idim-1, 
				  			astart + i * instep, 
				  			bstart + j * outstep,
				  			fstart + (i + j * iwidth) * fstep,
				  			firststride);
				  }
					i += 1;
				}
				j += 1;
			}
		} else {
		  var instep = 1;
			var ix = 1;
			while (ix < firststride) {
				instep *= idims(ix);
				ix += 1; 
			}
		  sgemmx(ORDER.ColMajor, TRANSPOSE.NoTrans, TRANSPOSE.NoTrans, outDims(0), instep/inDims(1), inDims(0)*inDims(1), 1f,
		      data, fstart, outDims(0), 
		      in.data, astart, outDims(0), 1f,
		      out.data, bstart, outDims(0));	  
		}
	}

	def _correlate(in:FND, out:FND, idim:Int, astart:Int, bstart:Int, fstart:Int) {
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
								fstart + (i + j * iwidth) * fstep);
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
				var j = 0;
				while (j < owidth) {                    // Move over output tensor
					var ss = 0f;
					val jfwidth = j * iwidth;
					var i = 0;
					while (i < iwidth) {                      // Move over input tensor
						ss += in.data(astart + ks + i) * data0(fstart + jfwidth + i);
						i += 1;
					}
					out.data(bstart + k + j) += ss;
					j += 1;
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
}

object FFilter {
	def FFilter1D(w:Int, nstride:Int, npad:Int, contents:ND) = {
		val inDims = irow(w);
		val outDims = irow(1);
		val stride = irow(nstride);
		val pad = irow(npad);
		if (contents.length != w) {
			throw new RuntimeException("FFilter1D bad initialization matrix")
		}
		val filt = FND.newOrCheckFND(irow(w), null, contents.GUID, "FFilter1D()".##);
		filt <-- contents;
		new FFilter(inDims, outDims, stride, pad, filt.data);
	}

	def FFilter1Dd(w:Int, din:Int, dout:Int, nstride:Int, npad:Int, contents:ND) = {
		val inDims = irow(din, w);
		val outDims = irow(dout, 1);
		val stride = irow(1, nstride);
		val pad = irow(1, npad);
		if (contents.length != w * din * dout) {
			throw new RuntimeException("FFilter1D bad initialization matrix")
		}
		val filt = FND.newOrCheckFND(irow(din*dout, w), null, contents.GUID, "FFilter1D()".##);
		filt <-- contents;
		new FFilter(inDims, outDims, stride, pad, filt.data);
	}

	def FFilter2D(w:Int, h:Int, nstride:Int, npad:Int, contents:ND) = {
		val inDims = irow(w, h);
		val outDims = irow(1, 1);
		val stride = irow(nstride, nstride);
		val pad = irow(npad, npad);
		if (contents.length != w*h) {
			throw new RuntimeException("FFilter2D bad initialization matrix")
		}
		val filt = FND.newOrCheckFND(irow(w, h), null, contents.GUID, "FFilter2D()".##);
		filt <-- contents;
		new FFilter(inDims, outDims, stride, pad, filt.data);
	}

	def FFilter2Dd(w:Int, h:Int, din:Int, dout:Int, nstride:Int, npad:Int, contents:ND) = {
		val inDims = irow(din, w, h);
		val outDims = irow(dout, 1, 1);
		val stride = irow(1, nstride, nstride);
		val pad = irow(1, npad, npad);
		if (contents.length != w*h*din*dout) {
			throw new RuntimeException("FFilter2D bad initialization matrix")
		}
		val filt = FND.newOrCheckFND(irow(din*dout, w, h), null, contents.GUID, "FFilter2D()".##);
		filt <-- contents;
		new FFilter(inDims, outDims, stride, pad, filt.data);
	}

}
