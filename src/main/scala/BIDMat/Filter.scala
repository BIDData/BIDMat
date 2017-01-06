
package BIDMat
import scala.util.hashing.MurmurHash3


trait Filter {
  val inDims:IMat;
  val outDims:IMat;
  val pad:IMat;
  val stride:IMat;
  val outPad:IMat;
  
  def computeFlops(in:ND, stride:IMat, pad:IMat):Long = {
    var i = 0;
    var flops = 2L;
    while (i < stride.length) {
      flops *= 1L * inDims(i) * outDims(i) * (in.dims(i) - inDims(i) + 1 + 2*pad(i))/stride(i)
      i += 1;
    }
    flops;
  }
  
  def notImplemented0ND(s:String):ND = { 
    throw new RuntimeException("operator "+s+" not implemented for Filter")
  }
  
  def notImplemented1ND(s:String, a:ND):ND = { 
    throw new RuntimeException("operator "+s+" not implemented for Filter and "+a.mytype)
  }
  
  def notImplemented2ND(s:String, a:ND, b:ND):ND = { 
    throw new RuntimeException("operator "+s+" not implemented for Filter and "+a.mytype+" and "+b.mytype)
  }
  
  def ^* (a:ND):ND = {
    notImplemented1ND("^*", a);
  }

}

object Filter {
  
  def getOutputDims(imageDims:IMat, finDims:IMat, foutDims:IMat, stride:IMat, pad:IMat, outPad:IMat, compress:Boolean = false):IMat = {
    val ilen = imageDims.length;
    if (finDims.length != ilen) throw new RuntimeException("getOutputDims Image and Filter number of dimensions mismatch %d %d" format (ilen, finDims.length))
    if (stride.length != ilen) throw new RuntimeException("getOutputDims Image and Filter Stride number of dimensions mismatch %d %d" format (ilen, stride.length))
    if (pad.length != ilen) throw new RuntimeException("getOutputDims Image and Filter Pad number of dimensions mismatch %d %d" format (ilen, pad.length))
    if (outPad.length != ilen) throw new RuntimeException("getOutputDims Image and Filter Pad number of dimensions mismatch %d %d" format (ilen, outPad.length))
    val odims = (imageDims + pad * 2 - finDims) / stride + foutDims + outPad * 2;
    (if (compress) ND.trimDims(odims) else odims);   
  }
  
  def getInputDims(imageDims:IMat, finDims:IMat, foutDims:IMat, stride:IMat, pad:IMat, outPad:IMat, compress:Boolean = false):IMat = {
    val ilen = imageDims.length;
    if (finDims.length != ilen) throw new RuntimeException("getInputDims Image and Filter number of dimensions mismatch %d %d" format (ilen, finDims.length))
    if (stride.length != ilen) throw new RuntimeException("getInputDims Image and Filter Stride number of dimensions mismatch %d %d" format (ilen, stride.length))
    if (pad.length != ilen) throw new RuntimeException("getInputDims Image and Filter Pad number of dimensions mismatch %d %d" format (ilen, pad.length))
    if (outPad.length != ilen) throw new RuntimeException("getInputDims Image and Filter Pad number of dimensions mismatch %d %d" format (ilen, outPad.length))
    val indims = (imageDims - foutDims - outPad*2 + 1) *@ stride - pad*2 + finDims - 1;
    (if (compress) ND.trimDims(indims) else indims);   
  }
  
  final val forward = 1;
  final val backwardGradient = 2;
  final val backwardModel = 3;

}