package BIDMat
import BIDMat.MatFunctions._

@SerialVersionUID(100L)
class FFilter(inDims0:IMat, outDims0:IMat, pad0:IMat, stride0:IMat, data0:Array[Float]) extends FND((inDims0 *@ outDims0).data, data0) with Filter {
	override val inDims = inDims0;
	override val outDims = outDims0;
	override val pad = pad0;
	override val stride = stride0;
  def convolve(a:FND):FND = {
    val outdims = Filter.getOutputDims(a.dims, inDims, outDims, stride, pad);
    val inpad = if (pad.asInstanceOf[AnyRef] != null) inDims + pad else inDims;
    val out = FND.newOrCheckFND(outdims, null, a.GUID, GUID, "convout".##);
    val in = if (pad.asInstanceOf[AnyRef] != null && SciFunctions.sum(pad).dv > 0) pad(a, pad, inpad) else a;
    out.clear;
    _convolve(in, out, inDims.length-1, 0, 0, 0);
    out
  }
  
  def _convolve(in:FND, out:FND, idim:Int, astart0:Int, bstart0:Int, fstart0:Int) {
    if (idim > 0) {
    	val idims = in.dims;
    	val odims = out.dims;
    	var astart = astart0;
    	var bstart = bstart0;
      var fstart = fstart0;
    	var i = 0;
    	if (idim > 0) {
    		var incol = 1;
    		var outcol = 1;
        var fcol = 1; 
    		var ix = 0;
    		while (ix < idim) {
    			incol *= idims(ix);
    			outcol *= odims(ix);
          fcol *= inDims(ix);
          fcol *= outDims(ix);
    			ix += 1; 
    		}
    		while (i < inDims(idim)) {
    			var j = 0;
    			while (j < outDims(idim)) {
    				_convolve(in, out, idim-1, astart, bstart + j*outcol, fstart);
    				j += 1;
    			}
    			astart += incol;
    			i + 1;
    		}
    	}
    }
  }
  
  def pad(a:FND, pad:IMat, inpad:IMat):FND = { 
    val out = FND.newOrCheckFND(inpad, null, a.GUID, GUID, "convpad".##);
    _pad(a, pad, out, inpad.length-1, 0, 0);
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

}