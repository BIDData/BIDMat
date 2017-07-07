package BIDMat
import MatFunctions._
import jcuda._
import jcuda.runtime._
import jcuda.runtime.JCuda._
import jcuda.jcudnn._
import jcuda.jcudnn.JCudnn._
import jcuda.runtime.cudaMemcpyKind._
import edu.berkeley.bid.MurmurHash3.MurmurHash3_x64_64

//
// Basic GPU convolutional Filter class.
//
// Currently a basic wrapper for NVIDIA CUDNN 4d filters (2 spatial dimensions). 
// Supported tensor orders are NHWC and NCHW. The NHWC order matches BIDMat's CPU implementation. 
//
// The filter is stored with the output dimension major, and the other dimensions in the same order as a data tensor.
// so a 4D filter block would be OHWC. This matches the ordering in NVIDIA CUDNN, and in BIDMat's CPU filter class. 
// 


@SerialVersionUID(100L)
class GFilter(inDims0:IMat, outDims0:IMat, stride0:IMat, pad0:IMat, outPad0:IMat, dataDims0:IMat, data0:Pointer) extends
  GMat(dataDims0.data, data0, dataDims0.data.reduce(_*_)) with Filter {
  
	override def mytype = "GFilter";

	val inDims = inDims0;
	val outDims = outDims0;
	val stride = if (stride0.asInstanceOf[AnyRef] != null) stride0 else iones(1, inDims.length);
	val pad = if (pad0.asInstanceOf[AnyRef] != null) pad0 else izeros(1,inDims.length);
	val outPad = if (outPad0.asInstanceOf[AnyRef] != null) outPad0 else izeros(1,inDims.length);
	val dataDims = dataDims0;
	var dataType = cudnnDataType.CUDNN_DATA_FLOAT;
	var tensorFormat = cudnnTensorFormat.CUDNN_TENSOR_NHWC;
  var convType = cudnnConvolutionMode.CUDNN_CROSS_CORRELATION;
//  var convType = cudnnConvolutionMode.CUDNN_CONVOLUTION;

  val fwdAlgo = Array(0);
  val bwdDataAlgo = Array(0);
  val bwdFilterAlgo = Array(0);
  var fwdTrained = false;
  var bwdDataTrained = false;
  var bwdFilterTrained = false;
  var adesc:cudnnTensorDescriptor = null;
  var bdesc:cudnnTensorDescriptor = null;
  var fdesc:cudnnFilterDescriptor = null;
  var convdesc:cudnnConvolutionDescriptor = null;
  
  var a:GMat = null;
  @volatile var workspaceFWD:GMat = null;
  @volatile var workspaceBWDdata:GMat = null;
  @volatile var workspaceBWDfilter:GMat = null;
  
  var cudnnMainHandle:cudnnHandle = null;
  var cudnn2ndHandle:cudnnHandle = null;
  var cudnnMainStream:cudaStream_t = null;
  var cudnn2ndStream:cudaStream_t = null;
  
  def initHandles() = {
  		cudnnMainHandle = new cudnnHandle;
  		cudnn2ndHandle = new cudnnHandle;
  		cudnnMainStream = new cudaStream_t;
  		cudnn2ndStream = new cudaStream_t;

  		var err = cudnnCreate(cudnnMainHandle);
  		if (err == 0) err = cudnnCreate(cudnn2ndHandle);
  		if (err == 0) err = cudaStreamCreate(cudnnMainStream);
  		if (err == 0) err = cudaStreamCreate(cudnn2ndStream);
  		if (err == 0) err = cudnnSetStream(cudnnMainHandle, cudnnMainStream);
  		if (err == 0) err = cudnnSetStream(cudnn2ndHandle, cudnn2ndStream);

  		if (err != 0) throw new RuntimeException("Error in CUDNN filter creation %s" format cudaGetErrorString(err));
  }
  
  initHandles();
  
  def setNHWC = {
		  tensorFormat = cudnnTensorFormat.CUDNN_TENSOR_NHWC;
  }
  
  def setNCHW = {
		  tensorFormat = cudnnTensorFormat.CUDNN_TENSOR_NCHW;
  }
  
  def setTensorFormat(tformat:Int) = {
    tformat match {
      case cudnnTensorFormat.CUDNN_TENSOR_NCHW => setNCHW;
      case cudnnTensorFormat.CUDNN_TENSOR_NHWC => setNHWC;
    }
  }
  
	def convolve(a:GMat, omat:Mat, doclear:Boolean, workspace:Mat):GMat = {
    val bdims = Filter.getOutputDims(a.dims, inDims, outDims, stride, pad, outPad);
    val hmm = ND.hashIMat(stride, ND.hashIMat(pad));
    val b = GMat.newOrCheckGMat(bdims, omat, a.GUID, GUID, hmm, "convout".##);
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
      if (fstatus > 0) throw new RuntimeException("Error creating filter tensor for forward convolution %d" format fstatus)
      
      val convdesc = new cudnnConvolutionDescriptor;
      cudnnCreateConvolutionDescriptor(convdesc);
      val cstatus = cudnnSetConvolution2dDescriptor(convdesc, pad(2), pad(1), stride(2), stride(1), 1, 1, convType);
      if (cstatus > 0) throw new RuntimeException("Error setting convolution descriptor for forward convolution %d" format cstatus);
      
      if (!fwdTrained) {
        val (preference, limit) = if (workspace.asInstanceOf[AnyRef] != null) {
          (cudnnConvolutionFwdPreference.CUDNN_CONVOLUTION_FWD_SPECIFY_WORKSPACE_LIMIT, workspace.length*4L);
        } else {
          (cudnnConvolutionFwdPreference.CUDNN_CONVOLUTION_FWD_PREFER_FASTEST, 0L);
        }
      	val gstatus = cudnnGetConvolutionForwardAlgorithm(cudnnMainHandle, adesc, fdesc, convdesc, bdesc, preference, limit, fwdAlgo);
      	if (gstatus > 0) throw new RuntimeException("Error getting best algorithm for forward convolution %d" format gstatus);
      	fwdTrained = true;
      }
       
      val _workspaceSizeInBytes = new Array[Long](1);
      var wserr = cudnnGetConvolutionForwardWorkspaceSize(cudnnMainHandle, adesc, fdesc, convdesc, bdesc, fwdAlgo(0), _workspaceSizeInBytes);
      val workspaceSizeInBytes = _workspaceSizeInBytes(0);
      workspaceFWD = if (workspace.asInstanceOf[AnyRef] != null) {
        workspace.asInstanceOf[GMat];
      } else {
        GMat.newOrCheckGMat((workspaceSizeInBytes/4).toInt, 1, null, GUID, a.GUID, hmm, "ConvFwdWS".##);
      }

      var err = cudnnConvolutionForward(cudnnMainHandle, GFilter.ONE, adesc, a.pdata, fdesc, pdata, convdesc, 
          fwdAlgo(0), workspaceFWD.pdata, workspaceSizeInBytes, if (doclear) GFilter.ZERO else GFilter.ONE, bdesc, b.pdata);
      
      cudaStreamSynchronize(cudnnMainStream);
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
	
	def convolve(a:GMat, omat:Mat, doclear:Boolean):GMat = convolve(a, omat, doclear, null);
  
  def convolve(a:GMat):GMat = convolve(a, null, true, null);
    
  def convolveT(deriv:GMat, inderiv:Mat, doclear:Boolean, workspace:Mat):GMat = {
  	val ddims = if (inderiv.asInstanceOf[AnyRef] != null) {                       // The in-out map with stride/pad is many to one, so not invertible always
  	  Filter.getOutputDims(inderiv.dims, inDims, outDims, stride, pad, outPad);   // Here check if the deriv/inputDeriv pair are compatible
  	} else {
  	  null;
  	}                                                                            // Otherwise compute an input size from the params
    val adims = if (ddims.asInstanceOf[AnyRef] != null && ND.compareDims(ddims.data, deriv.dims.data)){
      inderiv.dims;
    } else {
      Filter.getInputDims(deriv.dims, inDims, outDims, stride, pad, outPad);
    }
    val hmm = ND.hashIMat(stride, ND.hashIMat(pad));
    a = GMat.newOrCheckGMat(adims, inderiv, deriv.GUID, GUID, hmm, "convoutT".##);
    if (dims.length == 4) {
      val adesc = new cudnnTensorDescriptor;
      cudnnCreateTensorDescriptor(adesc);
      val astatus = cudnnSetTensor4dDescriptor(adesc, tensorFormat, dataType, a.dims(3), a.dims(0), a.dims(2), a.dims(1));
      if (astatus > 0) throw new RuntimeException("Error creating A tensor for backward data convolution %d, bad stride?" format astatus);
      
      val bdesc = new cudnnTensorDescriptor;
      cudnnCreateTensorDescriptor(bdesc);
      val bstatus = cudnnSetTensor4dDescriptor(bdesc, tensorFormat, dataType, deriv.dims(3), deriv.dims(0), deriv.dims(2), deriv.dims(1));
      if (bstatus > 0) throw new RuntimeException("Error creating B tensor for backward data convolution %d, bad stride?" format bstatus);
      
      val fdesc = new cudnnFilterDescriptor;
      cudnnCreateFilterDescriptor(fdesc);
      val fstatus = cudnnSetFilter4dDescriptor(fdesc, dataType, tensorFormat, outDims(0), inDims(0), inDims(2), inDims(1));
      if (fstatus > 0) throw new RuntimeException("Error creating filter tensor for backward data convolution %d" format fstatus);
      
      val convdesc = new cudnnConvolutionDescriptor;
      cudnnCreateConvolutionDescriptor(convdesc);
      val cstatus = cudnnSetConvolution2dDescriptor(convdesc, pad(2), pad(1), stride(2), stride(1), 1, 1, convType);
      if (cstatus > 0) throw new RuntimeException("Error setting convolution descriptor for backward data convolution %d" format cstatus);
      
      if (!bwdDataTrained) {
        val (preference, limit) = if (workspace.asInstanceOf[AnyRef] != null) {
          (cudnnConvolutionBwdDataPreference.CUDNN_CONVOLUTION_BWD_DATA_SPECIFY_WORKSPACE_LIMIT, workspace.length*4L);
        } else {
          (cudnnConvolutionBwdDataPreference.CUDNN_CONVOLUTION_BWD_DATA_PREFER_FASTEST, 0L);
        }
      	val gstatus = cudnnGetConvolutionBackwardDataAlgorithm(cudnnMainHandle, fdesc, bdesc, convdesc, adesc, preference, limit, bwdDataAlgo);
      	if (gstatus > 0) throw new RuntimeException("Error getting best algorithm for backward data convolution %d" format gstatus);
      	bwdDataTrained = true;
      }
 
      val _workspaceSizeInBytes = new Array[Long](1);
      var wserr = cudnnGetConvolutionBackwardDataWorkspaceSize(cudnnMainHandle, fdesc, bdesc, convdesc, adesc, bwdDataAlgo(0), _workspaceSizeInBytes);
      val workspaceSizeInBytes = _workspaceSizeInBytes(0);
      workspaceBWDdata =  if (workspace.asInstanceOf[AnyRef] != null) {
        workspace.asInstanceOf[GMat];
      } else {
        GMat.newOrCheckGMat((workspaceSizeInBytes/4).toInt, 1, null, GUID, a.GUID, hmm, "ConvBwdWS".##);
      }

      var err = cudnnConvolutionBackwardData(cudnnMainHandle, GFilter.ONE, fdesc, pdata, bdesc, deriv.pdata, convdesc, 
          bwdDataAlgo(0), workspaceBWDdata.pdata, workspaceSizeInBytes, if (doclear) GFilter.ZERO else GFilter.ONE, adesc, a.pdata);
      
      cudaStreamSynchronize(cudnnMainStream);
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

  def convolveT(a:GMat, b:Mat, doclear:Boolean):GMat = convolveT(a, b, doclear, null);
  
  def convolveT(a:GMat):GMat = convolveT(a, null, true, null);
  
  def convolveMfork(a:GMat, b:GMat, doclear:Boolean):GFilter= {
		val bdims = b.dims;
    val outdims = Filter.getOutputDims(a.dims, inDims, outDims, stride, pad, outPad);
    if ((bdims - outdims).data.exists(_ != 0)) {
      throw new RuntimeException("Output dimensions mismatch in convolveM")
    }
    val hmm = ND.hashIMat(stride, ND.hashIMat(pad));
    if (dims.length == 4) {
      adesc = new cudnnTensorDescriptor;
      cudnnCreateTensorDescriptor(adesc);
      val astatus = cudnnSetTensor4dDescriptor(adesc, tensorFormat, dataType, a.dims(3), a.dims(0), a.dims(2), a.dims(1));
      if (astatus > 0) throw new RuntimeException("Error creating A tensor for backward filter convolution %d, bad stride?" format astatus)
      
      bdesc = new cudnnTensorDescriptor;
      cudnnCreateTensorDescriptor(bdesc);
      val bstatus = cudnnSetTensor4dDescriptor(bdesc, tensorFormat, dataType, b.dims(3), b.dims(0), b.dims(2), b.dims(1));
      if (bstatus > 0) throw new RuntimeException("Error creating B tensor for backward filter convolution %d, bad stride?" format bstatus)
      
      fdesc = new cudnnFilterDescriptor;
      cudnnCreateFilterDescriptor(fdesc);
      val fstatus = cudnnSetFilter4dDescriptor(fdesc, dataType, tensorFormat, outDims(0), inDims(0), inDims(2), inDims(1));
      if (fstatus > 0) throw new RuntimeException("Error creating filter tensor for backward filter convolution %d" format fstatus)
      
      convdesc = new cudnnConvolutionDescriptor;
      cudnnCreateConvolutionDescriptor(convdesc);
      val cstatus = cudnnSetConvolution2dDescriptor(convdesc, pad(2), pad(1), stride(2), stride(1), 1, 1, convType);
      if (cstatus > 0) throw new RuntimeException("Error setting convolution descriptor for backward filter convolution %d" format cstatus);
      
      if (!bwdFilterTrained) {
      	val gstatus = cudnnGetConvolutionBackwardFilterAlgorithm(cudnn2ndHandle, adesc, bdesc, convdesc, fdesc, cudnnConvolutionBwdFilterPreference.CUDNN_CONVOLUTION_BWD_FILTER_PREFER_FASTEST, 0, bwdFilterAlgo);
      	if (gstatus > 0) throw new RuntimeException("Error getting best algorithm for backward filter convolution %d" format gstatus);
      	bwdFilterTrained = true;
      }
      
      val _workspaceSizeInBytes = new Array[Long](1);
      var wserr = cudnnGetConvolutionBackwardFilterWorkspaceSize(cudnn2ndHandle, adesc, bdesc, convdesc, fdesc, bwdFilterAlgo(0), _workspaceSizeInBytes);
      val workspaceSizeInBytes = _workspaceSizeInBytes(0);
      workspaceBWDfilter = GMat.newOrCheckGMat((workspaceSizeInBytes/4).toInt, 1, null, GUID, a.GUID, hmm, "ConvBwdFilterWS".##);

      var err = cudnnConvolutionBackwardFilter(cudnn2ndHandle, GFilter.ONE, adesc, a.pdata, bdesc, b.pdata, convdesc, 
          bwdFilterAlgo(0), workspaceBWDfilter.pdata, workspaceSizeInBytes, if (doclear) GFilter.ZERO else GFilter.ONE, fdesc, pdata);

      if (err > 0) throw new RuntimeException("Error in CUDNN backward filter convolution %s" format cudaGetErrorString(err));
      
    }
    Mat.nflops += computeFlops(a, stride, pad);
    this
  }
  
  def convolveMfork(a:GMat, b:GMat):GFilter = convolveMfork(a, b, true);
  
  override def convolveMjoin:GFilter = {
  	cudaStreamSynchronize(cudnn2ndStream);
  	val  err = cudaGetLastError();
    if (err > 0) throw new RuntimeException("Error in CUDNN backward data convolution %s" format cudaGetErrorString(err));
    cudnnDestroyConvolutionDescriptor(convdesc);
    cudnnDestroyFilterDescriptor(fdesc);
    cudnnDestroyTensorDescriptor(bdesc);
    cudnnDestroyTensorDescriptor(adesc);
  	this;
  }
  
  def convolveM(a:GMat, b:GMat, doclear:Boolean):GFilter = {
    convolveMfork(a, b, doclear);
    convolveMjoin;
  }
  
  def convolveM(a:GMat, b:GMat):GFilter = {
    convolveMfork(a, b, true);
    convolveMjoin;
  }
  
  override def convolve(b:Mat, omat:Mat, doclear:Boolean):Mat = {
			b match {
			case (bb:GMat) => convolve(bb, omat, doclear);
			}
	}
  
  override def convolve(b:Mat, omat:Mat, doclear:Boolean, workspace:Mat):Mat = {
			b match {
			case (bb:GMat) => convolve(bb, omat, doclear, workspace);
			}
	}

	override def convolveT(b:Mat, omat:Mat, doclear:Boolean):Mat = {
			b match {
			case (bb:GMat) => convolveT(bb, omat, doclear);
			}
	}
	
  override def convolveT(b:Mat, omat:Mat, doclear:Boolean, workspace:Mat):Mat = {
			b match {
			case (bb:GMat) => convolveT(bb, omat, doclear, workspace);
			}
	}

	override def convolveM(a:Mat, b:Mat, doclear:Boolean):Filter = {
			(a, b) match {
			case (aa:GMat, bb:GMat) => convolveM(aa, bb, doclear);
			}
	}
	
  override def convolveMfork(a:Mat, b:Mat, doclear:Boolean):Filter = {
			(a, b) match {
			case (aa:GMat, bb:GMat) => convolveMfork(aa, bb, doclear);
			}
	}
  
  override def copy:GFilter = {
		val out = new GFilter(inDims.copy, outDims.copy, stride.copy, pad.copy, outPad.copy, dataDims.copy, new Pointer);
		val len = 1L*length*Sizeof.FLOAT;
		if (Mat.debugMem) {
			 println("GFilter %d %d, %d %f" format (nrows, ncols, SciFunctions.getGPU, SciFunctions.GPUmem._1))
			 if (length > Mat.debugMemThreshold) throw new RuntimeException("GFilter alloc too large");
    }
		out.convType = convType;
		out.tensorFormat = tensorFormat;
		cudaMalloc(out.pdata, len);
		cudaStreamSynchronize(Mat.SyncMethod);
		cudaMemcpy(out.pdata, pdata, len, cudaMemcpyDeviceToDevice);
		cudaStreamSynchronize(Mat.SyncMethod);
		out;
	}
  
  def toNCHW:GFilter = {
		if (dims.length != 4) throw new RuntimeException("fromNHWCtoNCHW ndims must be 4");
		if (tensorFormat == cudnnTensorFormat.CUDNN_TENSOR_NCHW) {
			return this;
		} else {
			val tmp = transpose(MatFunctions.irow(1,2,0,3));
			val out = new GFilter(inDims.copy, outDims.copy, stride.copy, pad.copy, outPad.copy, dataDims.copy, tmp.pdata);
			out.setNCHW;
			out.setGUID(MurmurHash3_x64_64(Array(GUID), "fromNHWCtoNCHW".##));
			out
		}
	}
  
  def toNHWC:GFilter = {
    if (dims.length != 4) throw new RuntimeException("fromNCHWtoNHWC ndims must be 4");
    if (tensorFormat == cudnnTensorFormat.CUDNN_TENSOR_NHWC) {
    	return this;
    } else {
    	val tmp = reshapeView(MatFunctions.irow(dims(1), dims(2), dims(0), dims(3))).transpose(MatFunctions.irow(2,0,1,3));
    	val out = new GFilter(inDims.copy, outDims.copy, stride.copy, pad.copy, outPad.copy, dataDims.copy, tmp.pdata);
    	out.setNHWC;
    	out.setGUID(MurmurHash3_x64_64(Array(GUID), "fromNHWCtoNCHW".##));
    	out;
    }
  }
  
  override def * (a:GMat):GMat = {
			convolve(a);
	}

	override def ^* (a:GMat):GMat = {
			convolveT(a);
	}
		
	override def xavier(scale:Float):GFilter = GFilter.xavier(this, scale);
	
	override def xavier:GFilter = GFilter.xavier(this, 1f);

	override def transpose(p:IMat):GFilter = {
	  new GFilter(inDims, outDims, stride, pad, outPad, dataDims(dataDims.length-1) \ dataDims(0->(dataDims.length-1)), _transpose(p).pdata);
	}
}

object GFilter {

  var cudnnContextsInitialized = false;
  val ONE = Pointer.to(Array(1f));
  val ZERO = Pointer.to(Array(0f));
  
  def apply(a:FFilter):GFilter = {
    val outnd = GMat.newOrCheckGMat(a.dims, null, a.GUID, "GFilter".##);
    val out = new GFilter(a.inDims, a.outDims, a.stride, a.pad, a.outPad, a.dataDims, outnd.pdata);
    GMat.CPUtoGPUarraycopy(a.data, 0, out.pdata, 0, a.length, "GFilter apply");
    out.tensorFormat = a.tensorFormat;
    out.convType = a.convType;
    out.setGUID(MurmurHash3_x64_64(Array(a.GUID), "GFilter apply".##));
    out;
  }
  
  def apply(inDims:IMat,outDims:IMat, stride:IMat, pad:IMat, outPad:IMat, dataDims:IMat) = {
		  val len = dataDims.data.reduce(_*_);
		  if (Mat.debugMem) {
			  println("GFilter %d %d %f" format (len, SciFunctions.getGPU, SciFunctions.GPUmem._1))
			  if (len > Mat.debugMemThreshold) throw new RuntimeException("GFilter alloc too large");
		  }
		  val out = new GFilter(inDims, outDims, stride, pad, outPad, dataDims,new Pointer);
		  cudaMalloc(out.pdata, 1L*len*Sizeof.FLOAT);
		  cudaStreamSynchronize(Mat.SyncMethod);
		  out;
  }
  
  def GFilter1D(w:Int, nstride:Int, npad:Int, noutpad:Int):GFilter = {
    val inDims = irow(w);
    val outDims = irow(1);
    val stride = irow(nstride);
    val pad = irow(npad);
    val outPad = irow(noutpad);
    val out = GFilter(inDims, outDims, stride, pad, outPad, irow(w, 1));
    out
  }
  
  def GFilter1D(w:Int, nstride:Int, npad:Int):GFilter = GFilter1D(w, nstride, npad, 0);

  def GFilter1Dd(w:Int, din:Int, dout:Int, nstride:Int, npad:Int, noutpad:Int):GFilter = {
    val inDims = irow(din, w);
    val outDims = irow(dout, 1);
    val stride = irow(1, nstride);
    val pad = irow(0, npad);
    val outPad = irow(0, noutpad);
    val out = GFilter(inDims, outDims, stride, pad, outPad, irow(din, w, dout));
    out    
  }
  
  def GFilter1Dd(w:Int, din:Int, dout:Int, nstride:Int, npad:Int):GFilter = GFilter1Dd(w, din, dout, nstride, npad, 0);

  def GFilter2D(w:Int, h:Int, nstride:Int, npad:Int, noutpad:Int):GFilter = {
    val inDims = irow(w, h);
    val outDims = irow(1, 1);
    val stride = irow(nstride, nstride);
    val pad = irow(npad, npad);
    val outPad = irow(noutpad, noutpad);
    val out = GFilter(inDims, outDims, stride, pad, outPad, irow(w, h));
    out
  }
  
  def GFilter2D(w:Int, h:Int, nstride:Int, npad:Int):GFilter = GFilter2D(w, h, nstride, npad, 0);

  def GFilter2Dd(w:Int, h:Int, din:Int, dout:Int, nstride:Int, npad:Int, noutpad:Int):GFilter = {
    val inDims = irow(din, w, h);
    val outDims = irow(dout, 1, 1);
    val stride = irow(1, nstride, nstride);
    val pad = irow(0, npad, npad);
    val outPad = irow(0, noutpad, noutpad);
    val out = GFilter(inDims, outDims, stride, pad, outPad, irow(din, w, h, dout));
    out
  }
  
  def GFilter2Dd(w:Int, h:Int, din:Int, dout:Int, nstride:Int, npad:Int):GFilter = GFilter2Dd(w, h, dout, nstride, npad, 0);
  
  def GFilter2Ddn(w:Int, h:Int, din:Int, dout:Int, nstride:Int, npad:Int, noutpad:Int):GFilter = {
    val inDims = irow(din, w, h, 1);
    val outDims = irow(dout, 1, 1, 1);
    val stride = irow(1, nstride, nstride, 1);
    val pad = irow(0, npad, npad, 0);
    val outPad = irow(0, noutpad, noutpad, 0);
    val out = GFilter(inDims, outDims, stride, pad, outPad, irow(din, w, h, dout));
    out
  }
  
  def GFilter2Ddn(w:Int, h:Int, din:Int, dout:Int, nstride:Int, npad:Int):GFilter = GFilter2Ddn(w, h, din, dout, nstride, npad, 0);
  
  def xavier(f:GFilter, fscale:Float):GFilter = {
	  val scale = f.inDims.data.reduce(_*_);
	  GFunctions.normrnd(0, fscale/math.sqrt(scale).toFloat, f);
	  f;
	}

}
