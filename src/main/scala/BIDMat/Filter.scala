
package BIDMat


trait Filter {
  val inDims:IMat = null;
  val outDims:IMat = null;
  val pad:IMat = null;
  val stride:IMat = null;
  
}

object Filter {
  
  def getOutputDims(imageDims:IMat, filtInDims:IMat, filtOutDims:IMat, filtStride:IMat, filtPad:IMat, compress:Boolean = false):IMat = {
    val ilen = imageDims.length;
    if (filtInDims.length != ilen) throw new RuntimeException("getOutputDims Image and Filter number of dimensions mismatch %d %d" format (ilen, filtInDims.length))
    if (filtStride.asInstanceOf[AnyRef] != null && filtStride.length != ilen) throw new RuntimeException("getOutputDims Image and Filter Stride number of dimensions mismatch %d %d" format (ilen, filtInDims.length))
    if (filtPad.asInstanceOf[AnyRef] != null && filtPad.length != ilen) throw new RuntimeException("getOutputDims Image and Filter Pad number of dimensions mismatch %d %d" format (ilen, filtInDims.length))
    val fstride = if (filtStride.asInstanceOf[AnyRef] != null) filtStride else MatFunctions.iones(1, ilen);
    val fpad = if (filtPad.asInstanceOf[AnyRef] != null) filtPad else MatFunctions.izeros(1, ilen);
    val odims0 = (imageDims + fpad * 2 - filtInDims) / fstride + filtOutDims;
    (if (compress) ND.trimDims(odims0) else odims0);   
  }

}