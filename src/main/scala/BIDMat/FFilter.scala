package BIDMat
import BIDMat.MatFunctions._

@SerialVersionUID(100L)
class FFilter(inDims0:IMat, outDims0:IMat, stride0:IMat, pad0:IMat, data0:Array[Float]) extends FND((inDims0 *@ outDims0).data, data0) with Filter {
	override val inDims = inDims0;
	override val outDims = outDims0;
	override val pad = pad0;
	override val stride = stride0;
  
  def convolve(a:FND):FND = {
    val outdims = Filter.getOutputDims(a.dims, inDims, outDims, stride, pad);
    val out = FND.newOrCheckFND(outdims, null, a.GUID, GUID, "convout".##);
    val inpadded = if (pad.asInstanceOf[AnyRef] != null && SciFunctions.sum(pad).dv > 0) pad(a, pad) else a;
    out.clear;
    _convolve(inpadded, out, inDims.length-1, 0, 0, 0);
    out
  }
  
  def _convolve(in:FND, out:FND, idim:Int, astart:Int, bstart:Int, fstart:Int) {
	  val idims = in.dims;
	  val odims = out.dims;
	  if (idim > 0) {
		  var instep = 1;
		  var outstep = 1;
		  var fstep = 1;
		  val fwidth = inDims(idim);
		  var ix = 0;
		  while (ix < idim) {
			  instep *= idims(ix);
			  outstep *= odims(ix);
			  fstep *= inDims(ix);
			  fstep *= outDims(ix);
			  ix += 1; 
		  }
		  var k = 0;
		  while (k + fwidth - 1 < idims(idim)) {
			  var j = 0;
			  while (j < outDims(idim)) {
				  var i = 0;
				  while (i < fwidth) {
					  _convolve(in, out, idim-1, 
							  astart + (k + fwidth - i - 1) * instep, 
							  bstart + (k + outDims(idim) - j - 1) * outstep,
							  fstart + (i + j * fwidth) * fstep);
					  i += 1;
				  }
				  j + 1;
			  }
			  k += 1;
		  }
	  } else {
		  val fwidth = inDims(0); 
		  var k = 0;
		  while (k + fwidth - 1 < idims(0)) {           // Move forward over input+output tensors              
			  var j = 0;
			  while (j < outDims(0)) {                    // Move over output tensor
				  var i = 0;
				  var ss = 0f;
				  while (i < fwidth) {                      // Move over input tensor
					  ss += in.data(astart + k + fwidth - i - 1) * data0(fstart + j * fwidth + i);
					  i += 1;
				  }
				  out.data(bstart + k + outDims(0) - j - 1) += ss;
				  j += 1;
			  }
			  k += 1;
		  }
	  }
  }

  
  def pad(a:FND, pad:IMat):FND = { 
		val inpaddims = if (pad.asInstanceOf[AnyRef] != null) a.dims + pad*2 else a.dims;
    val out = FND.newOrCheckFND(inpaddims, null, a.GUID, GUID, "convpad".##);
    _pad(a, pad, out, inpaddims.length-1, 0, 0);
    out
  }
  
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
  }

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
  def FFilter1D(w:Int, nstride:Int, npad:Int, contents:FND) = {
    val inDims = irow(w);
    val outDims = irow(1);
    val stride = irow(nstride);
    val pad = irow(npad);
    if (contents.length != w) {
      throw new RuntimeException("FFilter1D bad initialization matrix")
    }
    val filt = FND.newOrCheckFND(irow(w), null, contents.GUID, "FFilter1D()".##);
    new FFilter(inDims, outDims, stride, pad, filt.data);
  }
  
  def FFilter2D(w:Int, h:Int, nstride:Int, npad:Int, contents:FND) = {
    val inDims = irow(w, h);
    val outDims = irow(1, 1);
    val stride = irow(nstride, nstride);
    val pad = irow(npad, npad);
    if (contents.length != w*h) {
      throw new RuntimeException("FFilter2D bad initialization matrix")
    }
    val filt = FND.newOrCheckFND(irow(w, h), null, contents.GUID, "FFilter2D()".##);
    new FFilter(inDims, outDims, stride, pad, filt.data);
  }

}