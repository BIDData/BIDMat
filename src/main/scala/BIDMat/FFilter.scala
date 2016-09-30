package BIDMat

@SerialVersionUID(100L)
class FFilter(inDims0:IMat, outDims0:IMat, pad0:IMat, stride0:IMat, data0:Array[Float]) extends FND((inDims0 *@ outDims0).data, data0) with Filter {
	override val inDims = inDims0;
	override val outDims = outDims0;
	override val pad = pad0;
	override val stride = stride0;
  def convolve(a:FND):FND = {
    val outdims = Filter.getOutputDims(a.dims, inDims, outDims, stride, pad);
    val out = FND.newOrCheckFND(outdims, null, a.GUID, GUID, 12345);
    out
  }
}

object FFilter {

}