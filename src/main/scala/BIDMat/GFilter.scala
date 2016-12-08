package BIDMat
import MatFunctions._
import jcuda._
import jcuda.runtime._
import jcuda.runtime.JCuda._
import jcuda.jcudnn._
import jcuda.jcudnn.JCudnn._

@SerialVersionUID(100L)
class GFilter(inDims0:IMat, outDims0:IMat, stride0:IMat, pad0:IMat, data0:Pointer) extends
  GND((inDims0 *@ outDims0).data, data0) with Filter {

	override val inDims = inDims0;
	override val outDims = outDims0;
	override val stride = if (stride0.asInstanceOf[AnyRef] != null) stride0 else iones(1, inDims.length);
	override val pad = if (pad0.asInstanceOf[AnyRef] != null) pad0 else izeros(1,inDims.length);
	var dataType = cudnnDataType.CUDNN_DATA_FLOAT;
	var tensorFormat = cudnnTensorFormat.CUDNN_TENSOR_NCHW;
  var convType = cudnnConvolutionMode.CUDNN_CROSS_CORRELATION;
  var fwdAlgo = cudnnConvolutionFwdAlgo.CUDNN_CONVOLUTION_FWD_ALGO_GEMM;
  var bwdFilterAlgo = cudnnConvolutionBwdFilterAlgo.CUDNN_CONVOLUTION_BWD_FILTER_ALGO_1;
  var bwdDataAlgo = cudnnConvolutionBwdDataAlgo.CUDNN_CONVOLUTION_BWD_DATA_ALGO_1;
  
	def convolve(a:GND, omat:ND, doclear:Boolean):GND = {
    val bdims = Filter.getOutputDims(a.dims, inDims, outDims, stride, pad);
    val hmm = Filter.hashIMat(stride, Filter.hashIMat(pad));
    val b = GND.newOrCheckGND(bdims, omat, a.GUID, GUID, hmm, "convout".##);
    if (dims.length == 4) {
      val adesc = new cudnnTensorDescriptor;
      cudnnCreateTensorDescriptor(adesc);
      val astatus = cudnnSetTensor4dDescriptor(adesc, tensorFormat, dataType, a.dims(3), a.dims(0), a.dims(2), a.dims(1));
      if (astatus > 0) throw new RuntimeException("Error creating A tensor for forward convolution %d, bad stride?" format astatus)
      
      val bdesc = new cudnnTensorDescriptor;
      cudnnCreateTensorDescriptor(bdesc);
      val bstatus = cudnnSetTensor4dDescriptor(bdesc, tensorFormat, dataType, b.dims(3), b.dims(0), b.dims(2), b.dims(1));
      if (bstatus > 0) throw new RuntimeException("Error creating B tensor for forward convolution %d, bad stride?" format bstatus)
      
      val fdesc = new cudnnFilterDescriptor;
      cudnnCreateFilterDescriptor(fdesc);
      val fstatus = cudnnSetFilter4dDescriptor(fdesc, dataType, tensorFormat, outDims(0), inDims(0), inDims(2), inDims(1));
      if (fstatus > 0) throw new RuntimeException("Error creating filter tensor for forward convolution %d" format bstatus)
      
      val convdesc = new cudnnConvolutionDescriptor;
      cudnnCreateConvolutionDescriptor(convdesc);
      val cstatus = cudnnSetConvolution2dDescriptor(convdesc, pad(2), pad(1), stride(2), stride(1), 1, 1, convType);
      if (cstatus > 0) throw new RuntimeException("Error setting convolution descriptor for forward convolution %d" format bstatus)
      
      val _workspaceSizeInBytes = new Array[Long](1);
      var wserr = cudnnGetConvolutionForwardWorkspaceSize(GFilter.getHandle, adesc, fdesc, convdesc, bdesc, fwdAlgo, _workspaceSizeInBytes);
      val workspaceSizeInBytes = _workspaceSizeInBytes(0);
      val workspace = GMat.newOrCheckGMat((workspaceSizeInBytes/4).toInt, 1, null, GUID, a.GUID, hmm, "ConvFwdWS".##);
      
//      println("workspace size = %d" format workspaceSizeInBytes)

      var err = cudnnConvolutionForward(GFilter.getHandle, GFilter.myone, adesc, a.data, fdesc, data, convdesc, 
          fwdAlgo, workspace.data, workspaceSizeInBytes, if (doclear) GFilter.myzero else GFilter.myone, bdesc, b.data);
      
      cudaDeviceSynchronize;
      if (err == 0) err = cudaGetLastError();
      if (err > 0) throw new RuntimeException("Error in CUDNN forward convolution %s" format cudaGetErrorString(err));
      
      cudnnDestroyConvolutionDescriptor(convdesc);
      cudnnDestroyFilterDescriptor(fdesc);
      cudnnDestroyTensorDescriptor(bdesc);
      cudnnDestroyTensorDescriptor(adesc);
    }
    Mat.nflops += computeFlops(a, stride, pad);
    b
  }
  
  def convolve(a:GND):GND = convolve(a, null, true);
    
  def convolveT(b:GND, omat:ND, doclear:Boolean):GND = {
    val adims = Filter.getInputDims(b.dims, inDims, outDims, stride, pad);
    val hmm = Filter.hashIMat(stride, Filter.hashIMat(pad));
    val a = GND.newOrCheckGND(adims, omat, b.GUID, GUID, hmm, "convoutT".##);
    if (dims.length == 4) {
      val adesc = new cudnnTensorDescriptor;
      cudnnCreateTensorDescriptor(adesc);
      val astatus = cudnnSetTensor4dDescriptor(adesc, tensorFormat, dataType, a.dims(3), a.dims(0), a.dims(2), a.dims(1));
      if (astatus > 0) throw new RuntimeException("Error creating A tensor for backward data convolution %d, bad stride?" format astatus)
      
      val bdesc = new cudnnTensorDescriptor;
      cudnnCreateTensorDescriptor(bdesc);
      val bstatus = cudnnSetTensor4dDescriptor(bdesc, tensorFormat, dataType, b.dims(3), b.dims(0), b.dims(2), b.dims(1));
      if (bstatus > 0) throw new RuntimeException("Error creating B tensor for backward data convolution %d, bad stride?" format bstatus)
      
      val fdesc = new cudnnFilterDescriptor;
      cudnnCreateFilterDescriptor(fdesc);
      val fstatus = cudnnSetFilter4dDescriptor(fdesc, dataType, tensorFormat, outDims(0), inDims(0), inDims(2), inDims(1));
      if (fstatus > 0) throw new RuntimeException("Error creating filter tensor for backward data convolution %d" format bstatus)
      
      val convdesc = new cudnnConvolutionDescriptor;
      cudnnCreateConvolutionDescriptor(convdesc);
      val cstatus = cudnnSetConvolution2dDescriptor(convdesc, pad(2), pad(1), stride(2), stride(1), 1, 1, convType);
      if (cstatus > 0) throw new RuntimeException("Error setting convolution descriptor for backward data convolution %d" format bstatus)
      
      val _workspaceSizeInBytes = new Array[Long](1);
      var wserr = cudnnGetConvolutionBackwardDataWorkspaceSize(GFilter.getHandle, fdesc, bdesc, convdesc, adesc, bwdDataAlgo, _workspaceSizeInBytes);
      val workspaceSizeInBytes = _workspaceSizeInBytes(0);
      val workspace = GMat.newOrCheckGMat((workspaceSizeInBytes/4).toInt, 1, null, GUID, a.GUID, hmm, "ConvBwdWS".##);
      
//      println("workspace size = %d" format workspaceSizeInBytes)

      var err = cudnnConvolutionBackwardData(GFilter.getHandle, GFilter.myone, fdesc, data, bdesc, b.data, convdesc, 
          bwdDataAlgo, workspace.data, workspaceSizeInBytes, if (doclear) GFilter.myzero else GFilter.myone, adesc, a.data);
      
      cudaDeviceSynchronize;
      if (err == 0) err = cudaGetLastError();
      if (err > 0) throw new RuntimeException("Error in CUDNN backward data convolution %s" format cudaGetErrorString(err));
      
      cudnnDestroyConvolutionDescriptor(convdesc);
      cudnnDestroyFilterDescriptor(fdesc);
      cudnnDestroyTensorDescriptor(bdesc);
      cudnnDestroyTensorDescriptor(adesc);
    }
    Mat.nflops += computeFlops(a, stride, pad);
    a
	}
  
  def convolveT(a:GND):GND = convolveT(a, null, true);
}

object GFilter {
  var cudnnContexts:Array[cudnnHandle] = null
  var cudnnContextsInitialized = false
  val _myone = new Array[Float](1);
  val _myzero = new Array[Float](1);
  _myone(0) = 1f;
  _myzero(0) = 0f;
  val myone = Pointer.to(_myone);
  val myzero = Pointer.to(_myzero);
  
  
  def apply(a:FFilter):GFilter = {
    val outnd = GND.newOrCheckGND(a.dims, null, a.GUID, "GFilter".##);
    outnd <-- a;
    val out = new GFilter(a.inDims, a.outDims, a.stride, a.pad, outnd.data);
    out;
  }
  
  def initHandles = {
		GFilter.synchronized { 
			if (!cudnnContextsInitialized) {
				val thisGPU = SciFunctions.getGPU
						val nGPUs = Mat.hasCUDA
						cudnnContexts = new Array[cudnnHandle](nGPUs)
						for (i <- 0 until nGPUs) {
							SciFunctions.setGPU(i)
							cudnnContexts(i) = new cudnnHandle()
							val err = cudnnCreate(cudnnContexts(i));
							if (err != 0) throw new RuntimeException("Cudnn initialization error %d on GPU %d" format (err, i));
						}  
				SciFunctions.setGPU(thisGPU)
				cudnnContextsInitialized = true
			}
		}
  }


  def getHandle = {
		if (!cudnnContextsInitialized) initHandles
		cudnnContexts(SciFunctions.getGPU)
  }
  
  def GFilter1D(w:Int, nstride:Int, npad:Int) = {
    val inDims = irow(w);
    val outDims = irow(1);
    val stride = irow(nstride);
    val pad = irow(npad);
    val out = new GFilter(inDims, outDims, stride, pad, new Pointer);
    cudaMalloc(out.data, 1L*w*Sizeof.FLOAT);
    out
  }

  def GFilter1Dd(w:Int, din:Int, dout:Int, nstride:Int, npad:Int) = {
    val inDims = irow(din, w);
    val outDims = irow(dout, 1);
    val stride = irow(1, nstride);
    val pad = irow(0, npad);
    val out = new GFilter(inDims, outDims, stride, pad, new Pointer);
    cudaMalloc(out.data, 1L*w*din*dout*Sizeof.FLOAT);
    out    
  }

  def GFilter2D(w:Int, h:Int, nstride:Int, npad:Int) = {
    val inDims = irow(w, h);
    val outDims = irow(1, 1);
    val stride = irow(nstride, nstride);
    val pad = irow(npad, npad);
    val out = new GFilter(inDims, outDims, stride, pad, new Pointer);
    cudaMalloc(out.data, 1L*w*h*Sizeof.FLOAT);
    out
  }

  def GFilter2Dd(w:Int, h:Int, din:Int, dout:Int, nstride:Int, npad:Int) = {
    val inDims = irow(din, w, h);
    val outDims = irow(dout, 1, 1);
    val stride = irow(1, nstride, nstride);
    val pad = irow(0, npad, npad);
    val out = new GFilter(inDims, outDims, stride, pad, new Pointer);
    cudaMalloc(out.data, 1L*din*dout*w*h*Sizeof.FLOAT);
    out
  }
  
  def GFilter2Ddn(w:Int, h:Int, din:Int, dout:Int, nstride:Int, npad:Int) = {
    val inDims = irow(din, w, h, 1);
    val outDims = irow(dout, 1, 1, 1);
    val stride = irow(1, nstride, nstride, 1);
    val pad = irow(0, npad, npad, 0);
    val out = new GFilter(inDims, outDims, stride, pad, new Pointer);
    cudaMalloc(out.data, 1L*din*dout*w*h*Sizeof.FLOAT);
    out
  }

}