package BIDMat;
import BIDMat.MatFunctions._;
import edu.berkeley.bid.CBLAS._;
import edu.berkeley.bid.UTILS
import SciFunctions._;
import edu.berkeley.bid.MurmurHash3.MurmurHash3_x64_64
import jcuda.jcudnn._
import jcuda.jcudnn.JCudnn._

//
// Basic CPU convolutional Filter class.
//
// Works in any dimension, but the last three (minor) dimensions should be HWC in descending stride order, so C has unit stride. 
// e.g. A 4D data tensor would be NHWC (Channel minor, N(minibatch number) major), which agrees with Tensorflow, and is one of the
// supported formats in NVIDIA CUDNN.
//
// The filter is stored with the output dimension major, and the other dimensions in the same order as a data tensor.
// so a 4D filter block would be OHWC. This matches the ordering in NVIDIA CUDNN.
// 
// There are four algorithms:
// * Element-wise - pure Scala implementation
// * im2col - copy input elements into a im2col matrix and do a single multiply.
// * 1x1 convolutions. For a kxk filter, perform k^2 1x1 convolutions each with a single matrix multiply. Strides must all be 1.
// * kx1 convolutions. For a kxk filter, perform multiple kx1 gemm convolutions over the W dimension only. Works directly on a
//   padded input array and can handle strides > 1. 


@SerialVersionUID(100L)
class FFilter(inDims0:IMat, outDims0:IMat, stride0:IMat, pad0:IMat, outPad0:IMat, dataDims0:IMat, data0:Array[Float]) extends
FMat(dataDims0.data, data0) with Filter {
  
	override def mytype = "FFilter";

	val inDims = inDims0;
	val outDims = outDims0;
	val stride = if (stride0.asInstanceOf[AnyRef] != null) stride0 else iones(1, inDims.length);
	val pad = if (pad0.asInstanceOf[AnyRef] != null) pad0 else izeros(1,inDims.length);
	val outPad = if (outPad0.asInstanceOf[AnyRef] != null) outPad0 else izeros(1,inDims.length);
	val dataDims = dataDims0;
	var tensorFormat = cudnnTensorFormat.CUDNN_TENSOR_NHWC;
	var convType = cudnnConvolutionMode.CUDNN_CROSS_CORRELATION;

	var timer = 0f;

	def convolve(a:FMat, omat:Mat, doclear:Boolean):FMat = {
			val bdims = Filter.getOutputDims(a.dims, inDims, outDims, stride, pad, outPad);
			val hmm = ND.hashIMat(stride, ND.hashIMat(pad));
			val b = FMat.newOrCheckFMat(bdims, omat, a.GUID, GUID, hmm, "convout".##);
			val apadmat = if (pad.data.exists(_ != 0)) {
				val m = FMat.newOrCheckFMat(a.dims + pad * 2, null, a.GUID, GUID, hmm, "convinpad".##);
				m.clear;
				_copy_padded(a, m, pad, inDims.length-1, 0, 0, true);
				m
			} else a;
			if (Mat.useMKL && inDims(0) == a.dims(0) && outDims(0) == bdims(0)) {             // Use gemm acceleration
				if (inDims(0) < FFilter.im2colThreshold) {                                      // shallow input matrix, do kx1 convolutions
					if (doclear) b.clear; 
					_fast_convolve2(apadmat, b, inDims.length-1, 0, 0, 0, Filter.forward);       
				} else if (stride.data.exists(_ > 1)) {                                         // input is deep but strided, use im2col
					val cellsize = SciFunctions.prod(inDims).v;
					val outlength = SciFunctions.prod(bdims(1 -> bdims.length)).v;
					val i2cmat = FMat.newOrCheckFMat(cellsize, outlength, null, a.GUID, GUID, hmm, "convIm2col".##);
//					_im2col(apadmat, b, i2cmat, inDims.length-1, cellsize, 0, 0);
					UTILS.im2col(inDims.length, inDims.data, outDims.data, apadmat.dims.data, apadmat.data, b.dims.data, b.data, i2cmat.data, stride.data);
					if (doclear) b.clear;
					//        println("convolve %d x %d by %d x %d " format (bdims(0), cellsize, cellsize, outlength));
					sgemm(ORDER.ColMajor, TRANSPOSE.Trans, TRANSPOSE.NoTrans, bdims(0), outlength, cellsize, 1f,
							data, length/bdims(0), i2cmat.data, cellsize, 1f, b.data, bdims(0));
				} else {                                                                         // do 1x1 convolutions
					val ddiff = (inDims - 1)/stride - outDims + 1;
					ddiff(0) = 0;
					val bpadmat = if (ddiff.data.exists(_ != 0)) {
						FMat.newOrCheckFMat(b.dims + ddiff, null, a.GUID, GUID, hmm, "convoutpad".##);
					} else {
						b;
					} 
					if (doclear) { 
						bpadmat.clear;
					} else {
						_copy_padded(b, bpadmat, ddiff/2, inDims.length-1, 0, 0, true);
					}
					_fast_convolve(apadmat, bpadmat, inDims.length-1, 0, 0, 0, Filter.forward);
					_copy_padded(bpadmat, b, ddiff/2, inDims.length-1, 0, 0, false);
				}
			} else {
				if (doclear) b.clear; 
				_convolve(apadmat, b, inDims.length-1, 0, 0, 0, Filter.forward);
			}
			Mat.nflops += computeFlops(a, stride, pad);
			b;
	};

	def convolve(a:FMat):FMat = convolve(a, null, true);

	def convolveT(b:FMat, omat:Mat, doclear:Boolean):FMat = {
			val bdims = b.dims;
			val adims = Filter.getInputDims(b.dims, inDims, outDims, stride, pad, outPad);
			val hmm = ND.hashIMat(stride, ND.hashIMat(pad));
			val a = FMat.newOrCheckFMat(adims, omat, b.GUID, GUID, hmm, "convTin".##);
			val apadmat = if (pad.data.exists(_ != 0)) {
				val m = FMat.newOrCheckFMat(a.dims + pad * 2, null, b.GUID, GUID, hmm, "convTinpad".##);
				if (!doclear) _copy_padded(a, m, pad, inDims.length-1, 0, 0, true);
				m
			} else a;
			if (doclear) apadmat.clear;
			val ddiff = (inDims - 1)/stride - outDims + 1;
			ddiff(0) = 0;
			if (Mat.useMKL && inDims(0) == adims(0) && outDims(0) == b.dims(0)) {             // Use gemm acceleration
				if (inDims(0) < FFilter.im2colThreshold) {                                      // shallow input matrix, do kx1 convolutions 
					_fast_convolve2(apadmat, b, inDims.length-1, 0, 0, 0, Filter.backwardGradient); 
				} else if (stride.data.exists(_ > 1)) {                                         // if strided use im2col
					val cellsize = SciFunctions.prod(inDims).v;
					val outlength = SciFunctions.prod(bdims(1 -> bdims.length)).v;
					val i2cmat = FMat.newOrCheckFMat(cellsize, outlength, null, a.GUID, GUID, hmm, "convTIm2col".##);
					sgemm(ORDER.ColMajor, TRANSPOSE.NoTrans, TRANSPOSE.NoTrans, cellsize, outlength, bdims(0), 1f,
							data, length/bdims(0), b.data, bdims(0), 0f, i2cmat.data, cellsize);
					UTILS.col2im(inDims.length, inDims.data, outDims.data, apadmat.dims.data, apadmat.data, b.dims.data, b.data, i2cmat.data, stride.data);
					//		_col2im(apadmat, b, i2cmat, inDims.length-1, cellsize, 0, 0);
				} else {                                                                        // use 1x1 convolutions
					val bpadmat = if (ddiff.data.exists(_ != 0)) {
						val m = FMat.newOrCheckFMat(b.dims + ddiff, null, b.GUID, GUID, hmm, "convToutpad".##);
						m.clear;
						_copy_padded(b, m, ddiff/2, inDims.length-1, 0, 0, true);
						m;
					} else {
						b;
					}
					_fast_convolve(apadmat, bpadmat, inDims.length-1, 0, 0, 0, Filter.backwardGradient);
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

	def convolveT(a:FMat):FMat = convolveT(a, null, true);

	def convolveM(a:FMat, b:FMat, doclear:Boolean):FMat = {
			val bdims = b.dims;
			val outdims = Filter.getOutputDims(a.dims, inDims, outDims, stride, pad, outPad);
			if ((bdims - outdims).data.exists(_ != 0)) {
				throw new RuntimeException("Output dimensions mismatch in convolveM")
			}
			val hmm = ND.hashIMat(stride, ND.hashIMat(pad));
			val apadmat = if (pad.data.exists(_ != 0)) {
				val m = FMat.newOrCheckFMat(a.dims + pad * 2, null, a.GUID, b.GUID, hmm, "convMinpad".##);
				m.clear;
				_copy_padded(a, m, pad, inDims.length-1, 0, 0, true);
				m
			} else a;
			if (doclear) clear;
			if (Mat.useMKL && inDims(0) == a.dims(0) && outDims(0) == outdims(0)) {             // Use gemm acceleration
				if (inDims(0) < FFilter.im2colThreshold) {                                      // shallow input matrix, do 0-1 dim gemm 
					_fast_convolve2(apadmat, b, inDims.length-1, 0, 0, 0, Filter.backwardModel); 
				} else if (stride.data.exists(_ > 1)) {
					val cellsize = SciFunctions.prod(inDims).v;
					val outlength = SciFunctions.prod(bdims(1 -> bdims.length)).v;
					val i2cmat = FMat.newOrCheckFMat(cellsize, outlength, null, a.GUID, GUID, hmm, "convMIm2col".##);
//					_im2col(apadmat, b, i2cmat, inDims.length-1, cellsize, 0, 0);
					UTILS.im2col(inDims.length, inDims.data, outDims.data, apadmat.dims.data, apadmat.data, b.dims.data, b.data, i2cmat.data, stride.data);
					sgemm(ORDER.ColMajor, TRANSPOSE.NoTrans, TRANSPOSE.Trans, cellsize, bdims(0), outlength, 1f,
							i2cmat.data, cellsize, b.data, bdims(0), if (doclear) 0f else 1f, data, length/bdims(0));
				} else {
					val outdiff = (inDims - 1)/stride - outDims + 1;
					outdiff(0) = 0;
					val bpadmat = if (outdiff.data.exists(_ != 0)) {
						val m = FMat.newOrCheckFMat(b.dims + outdiff, null, a.GUID, b.GUID, hmm, "convMoutpad".##);
						m.clear;
						_copy_padded(b, m, outdiff/2, inDims.length-1, 0, 0, true);
						m;
					} else {
						b;
					}
					_fast_convolve(apadmat, bpadmat, inDims.length-1, 0, 0, 0, Filter.backwardModel);
				}
			} else {
				_convolve(apadmat, b, inDims.length-1, 0, 0, 0, Filter.backwardModel);
			}
			Mat.nflops += computeFlops(a, stride, pad);
			this;
	};

	def convolveM(a:FMat, b:FMat):FMat = convolveM(a, b, true);

	override def copy:FFilter = {
		val a = new FFilter(inDims.copy, outDims.copy, stride.copy, pad.copy, outPad.copy, dataDims.copy, new Array[Float](length));
		System.arraycopy(data, 0, a.data, 0, length)
		a;
	}
	
	override def convolve(b:Mat, omat:Mat, doclear:Boolean):Mat = {
			b match {
			case (bb:FMat) => convolve(bb, omat, doclear);
			}
	}

	override def convolveT(b:Mat, omat:Mat, doclear:Boolean):Mat = {
			b match {
			case (bb:FMat) => convolveT(bb, omat, doclear);
			}
	}

	override def convolveM(a:Mat, b:Mat, doclear:Boolean):Mat = {
			(a, b) match {
			case (aa:FMat, bb:FMat) => convolveM(aa, bb, doclear);
			}
	}
	// Convolution using global 1x1 Convolutions. 

	def _fast_convolve(a:FMat, b:FMat, idim:Int, astart:Int, bstart:Int, fstart:Int, convType:Int) {
		val adims = a.dims;
		val bdims = b.dims;
		val iwidth = inDims(idim);
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
				ix += 1; 
			}
			var i = 0;
			while (i < iwidth) {
				_fast_convolve(a, b, idim-1, 
						astart + i * astep, 
						bstart + pad(idim) * bstep,
						fstart + i * fstep,
						convType);			
				i += 1;
			}
		} else {
			var mstep = 1;
			val owidth = outDims(idim);
			var ix = 1;
			while (ix < bdims.length) {
				mstep *= bdims(ix);
				ix += 1; 
			}
			val bstart0 = bstart + pad(0);
			val mstep0 = math.min(mstep, math.min((a.length - astart)/(iwidth*stride(1)), (b.length - bstart0)/owidth));
			convType match {			  
			case Filter.forward => {
				sgemmx(ORDER.ColMajor, TRANSPOSE.Trans, TRANSPOSE.NoTrans, owidth, mstep0, iwidth, 1f,
						data, fstart, length/owidth, 
						a.data, astart, iwidth*stride(1), 1f,
						b.data, bstart0, owidth);	
			}
			case Filter.backwardGradient => {
				sgemmx(ORDER.ColMajor, TRANSPOSE.NoTrans, TRANSPOSE.NoTrans, iwidth, mstep0, owidth, 1f,
						data, fstart, length/owidth, 
						b.data, bstart0, owidth, 1f,
						a.data, astart, iwidth*stride(1));	
			}
			case Filter.backwardModel => {
				sgemmx(ORDER.ColMajor, TRANSPOSE.NoTrans, TRANSPOSE.Trans, iwidth, owidth, mstep0, 1f, 
						a.data, astart, iwidth*stride(1), 
            b.data, bstart0, owidth, 1f,
						data, fstart, length/owidth);	
			}
			}
		}
	}
	//
	// Convolve into unpadded b
	//

	def _fast_convolve2(a:FMat, b:FMat, idim:Int, astart:Int, bstart:Int, fstart:Int, convType:Int) {
		//    println("rec %d %d %d %d" format (idim, astart, bstart, fstart))
		val adims = a.dims;
		val bdims = b.dims;
		val iwidth = inDims(idim);
		val owidth = outDims(idim);
		val kstride = stride(idim);
		if (idim > 1) {
			var astep = 1;
			var bstep = 1;
			var fstep = 1;
			var ix = 0;
			while (ix < idim) {
				astep *= adims(ix);
				bstep *= bdims(ix);
				fstep *= inDims(ix);
				ix += 1; 
			}
			var i = 0;
			while (i < iwidth) {
				var k = 0;
				var ks = 0;
				while (k < bdims(idim)) {
					_fast_convolve2(a, b, idim-1, 
							astart + (i + ks) * astep, 
							bstart + k * bstep,
							fstart + i * fstep,
							convType);  
					k += 1;
					ks += kstride;
				}
				i += 1;
			}
		} else {
			val astep = adims(0);
			val bstep = bdims(0);
			val fstep = inDims(0);
			val iwidth0 = inDims(0);
			val owidth0 = outDims(0);
			val nstrides = 1 + (iwidth - 1)/kstride;
			val iwidthUp = nstrides * kstride;
			//      println("widths %d %d" format (nstrides, iwidthUp))
			var k = 0;
			var ks = 0;
			while (ks < iwidth) {
				val astart1 = astart + ks * astep;
				val bstart1 = bstart + k * bstep;
				val mstep0 = math.min((adims(1)-ks)/iwidthUp, 1 + (bdims(1)-k-1)/nstrides);
				//          val mstep0 = (adims(1)-ks)/widthUp;
				convType match {        
				case Filter.forward => {
					//            println("%d %d %d, %d %d, %d %d, %d %d" format (owidth0, mstep0, iwidth0*iwidth, 
					//            		fstart1, owidth0*iwidth, astart1, iwidth0*iwidthUp, bstart1, owidth0*nstrides));
					val t0 = toc
							sgemmx(ORDER.ColMajor, TRANSPOSE.Trans, TRANSPOSE.NoTrans, owidth0, mstep0, iwidth0*iwidth, 1f,
									data, fstart, length/owidth0, 
									a.data, astart1, iwidth0*iwidthUp, 1f,
									b.data, bstart1, owidth0*nstrides); 
					val t1 = toc
							timer += t1 - t0;
				}
				case Filter.backwardGradient => {
					sgemmx(ORDER.ColMajor, TRANSPOSE.NoTrans, TRANSPOSE.NoTrans, iwidth0*iwidth, mstep0, owidth0, 1f,
							data, fstart, length/owidth0, 
							b.data, bstart1, owidth0*nstrides, 1f,
							a.data, astart1, iwidth0*iwidthUp);  
				}
				case Filter.backwardModel => {
					sgemmx(ORDER.ColMajor, TRANSPOSE.NoTrans, TRANSPOSE.Trans, iwidth0*iwidth, owidth0, mstep0, 1f, 
							a.data, astart1, iwidth0*iwidthUp, 
							b.data, bstart1, owidth0*nstrides, 1f, 
							data, fstart, length/owidth0);  
				}
				}
				k += 1;
				ks += kstride;
			}
		}
	}

	// copy into padded array, or from padded array

	def _copy_padded(in:FMat, out:FMat, padx:IMat, idim:Int, astart:Int, bstart:Int, topadded:Boolean) {
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

	def _convolve(a:FMat, b:FMat, idim:Int, astart:Int, bstart:Int, fstart:Int, convType:Int) {
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
				ix += 1; 
			}
			var k = 0;
			var ks = 0;
			while (k < bdims(idim)) {
				var i = 0;
				while (i < iwidth) {
					_convolve(a, b, idim-1, 
							astart + (ks + i) * astep, 
							bstart + k * bstep,
							fstart + i * fstep,
							convType);
					i += 1;
				}
				k += 1;
				ks += kstride;
			}
		} else {
			var k = 0;
			var ks = 0;
			val fstep = length/outDims(0);
			while (k + owidth - 1 < bdims(0)) {              // Move forward over input+output tensors              
				convType match {
				case Filter.forward => {
					var i = 0;
					while (i < iwidth) {                      // Move over input tensor
						val adat = a.data(astart + ks + i);
						var j = 0;
						var js = 0;
						while (j < owidth) {                        // Move over output tensor
							b.data(bstart + k + j) += data0(fstart + i + js) * adat;
							j += 1;
							js += fstep;
						}
						i += 1;
					}
				}
				case Filter.backwardGradient => {
					var i = 0;
					while (i < iwidth) {
						var ss = 0f;
						var j = 0;
						var js = 0;
						while (j < owidth) {
							ss += b.data(bstart + k + j) * data0(fstart + i + js);
							j += 1;
							js += fstep;
						}
						a.data(astart + ks + i) += ss;
						i += 1;
					}
				}
				case Filter.backwardModel => {
					var i = 0;
					while (i < iwidth) {
						val adat = a.data(astart + ks + i);
						var j = 0;
						var js = 0;
						while (j < owidth) {
							data0(fstart + i + js) += b.data(bstart + k + j) * adat;
							j += 1;
							js += fstep;
						}
						i += 1;
					}
				}
				}
				k += 1;
				ks += kstride;
			}
		}
	};
  
	def _im2col(a:FMat, b:FMat, i2c:FMat, idim:Int, celldim:Int, astart:Int, bstart:Int):FMat = {
			val adims = a.dims;
			val bdims = b.dims;
			val iwidth = inDims(idim);
			val owidth = outDims(idim);
			if (idim > 2) {
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
			} else if (idim == 2) {
				val astep1 = adims(0)*adims(1);
				val bstep1 = bdims(1);
				val cellstep1 = inDims(0)*inDims(1);
				var iin1 = 0;
				var iout1 = 0;
				while (iout1 < bdims(2)) {     
					var astart1 = astart + astep1 * iin1;
					var bstart1 = bstart + celldim * bstep1 * iout1;
					var j1 = 0;
					while (j1 < iwidth) {     
						var	astep0 = adims(0);
						var	cellstep0 = inDims(0);
						var iin0 = 0;
						var iout0 = 0;
						while (iout0 < bdims(1)) {
							var j0 = 0;
							var astart0 = astart1 + astep0 * iin0;
							var bstart0 = bstart1 + celldim * iout0;
							while (j0 < inDims(1)) {
								var k = 0;
								while (k < cellstep0) {
									i2c.data(bstart0 + k) = a.data(astart0 + k);
									k += 1;
								}
								astart0 += astep0;
								bstart0 += cellstep0;
								j0 += 1;
							}
							iout0 += 1;
							iin0 += stride(1);
						}
						astart1 += astep1;
						bstart1 += cellstep1;
						j1 += 1;
					}
					iout1 += 1;
					iin1 += stride(2);
				}
			} else if (idim == 1) {
				var astep = adims(0);
				var cellstep = inDims(0);
				var iin = 0;
				var iout = 0;
				while (iout < bdims(idim)) {
					var j = 0;
					var astart0 = astart + astep * iin;
					var bstart0 = bstart + celldim * iout;
					while (j < iwidth) {
						var k = 0;
						while (k < cellstep) {
							i2c.data(bstart0 + k) = a.data(astart0 + k);
							k += 1;
						}
						astart0 += astep;
						bstart0 += cellstep;
						j += 1;
					}
					iout += 1;
					iin += stride(idim);
				}    
			} else {
				if (iwidth > FFilter.arraycopy) {
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

    def _im2colpad(a:FMat, b:FMat, i2c:FMat, idim:Int, celldim:Int, astart:Int, bstart:Int, zerofill:Boolean, inparallel:Boolean):FMat = {
      val adims = a.dims;
      val bdims = b.dims;
      val iwidth = inDims(idim);
      val owidth = outDims(idim);
      if (idim > 2) {
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
        var iin = -pad(idim);
        var iout = 0;
        if (! inparallel && bdims(idim) > 4) {
          (0 until bdims(idim)).par.foreach((iout:Int) =>{
            val iin = -pad(idim) + iout * stride(idim); 
            var j = 0;
            while (j < iwidth) {
              val zerofilln = zerofill || iin + j < 0 || iin + j >= adims(idim);
              _im2colpad(a, b, i2c, idim - 1, celldim, astart + astep * (iin + j), bstart + celldim * bstep * iout + cellstep * j, zerofilln, true);
              j += 1;
            }
          });          
        } else {
          (0 until bdims(idim)).foreach((iout:Int) =>{
            val iin = -pad(idim) + iout * stride(idim);
            var j = 0;
            while (j < iwidth) {
              val zerofilln = zerofill || iin + j < 0 || iin + j >= adims(idim);
              _im2colpad(a, b, i2c, idim - 1, celldim, astart + astep * (iin + j), bstart + celldim * bstep * iout + cellstep * j, zerofilln, false);
              j += 1;
            }
          });
        }
      } else if (idim == 2) {
        val astep1 = adims(0)*adims(1);
        val bstep1 = bdims(1);
        val cellstep1 = inDims(0)*inDims(1);
        var iin1 = -pad(2);
        var iout1 = 0;
        while (iout1 < bdims(2)) {     
          var astart1 = astart + astep1 * iin1;
          var bstart1 = bstart + celldim * bstep1 * iout1;
          var j1 = 0;
          while (j1 < iwidth) {     
            val zerofill1 = zerofill || iin1 + j1 < 0 || iin1 + j1 >= adims(2);
            var astep0 = adims(0);
            var cellstep0 = inDims(0);
            var iin0 = -pad(1);
            var iout0 = 0;
            while (iout0 < bdims(1)) {
              var j0 = 0;
              var astart0 = astart1 + astep0 * iin0;
              var bstart0 = bstart1 + celldim * iout0;
              while (j0 < inDims(1)) {
                val zerofill0 = zerofill1 || iin0 + j0 < 0 || iin0 + j0 >= adims(1);
                if (zerofill0) {
                  var k = 0;
                  while (k < cellstep0) {
                    i2c.data(bstart0 + k) = 0f;
                    k += 1;
                  } 
                } else {
                  var k = 0;
                  while (k < cellstep0) {
                    i2c.data(bstart0 + k) = a.data(astart0 + k);
                    k += 1;
                  }
                }
                astart0 += astep0;
                bstart0 += cellstep0;
                j0 += 1;
              }
              iout0 += 1;
              iin0 += stride(1);
            }
            astart1 += astep1;
            bstart1 += cellstep1;
            j1 += 1;
          }
          iout1 += 1;
          iin1 += stride(2);
        }
      } else if (idim == 1) {
        var astep = adims(0);
        var cellstep = inDims(0);
        var iin = 0;
        var iout = 0;
        while (iout < bdims(idim)) {
          var j = 0;
          var astart0 = astart + astep * iin;
          var bstart0 = bstart + celldim * iout;
          while (j < iwidth) {
            val zerofill0 = zerofill || iin + j < 0 || iin + j >= adims(1);
            var k = 0;
            if (zerofill0) {
              while (k < cellstep) {
                i2c.data(bstart0 + k) = 0f;
                k += 1;
              }
            } else {
              while (k < cellstep) {
                i2c.data(bstart0 + k) = a.data(astart0 + k);
                k += 1;
              }
            }
            astart0 += astep;
            bstart0 += cellstep;
            j += 1;
          }
          iout += 1;
          iin += stride(idim);
        }    
      } else {
        if (iwidth > FFilter.arraycopy) {
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
/*	def _im2col(a:FMat, b:FMat, i2c:FMat, idim:Int, celldim:Int, astart:Int, bstart:Int, zerofill:Boolean):FMat = {
			val adims = a.dims;
			val bdims = b.dims;
			val iwidth = inDims(idim);
			val owidth = outDims(idim);
			if (idim > 2) {
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
//				var iin = -pad(idim);
				var iin = 0;
				var iout = 0;
				while (iout < bdims(idim)) {
					var j = 0;
					while (j < iwidth) {
            val nextzerofill = zerofill || iin + j < 0 || iin + j >= adims(idim);
						_im2col(a, b, i2c, idim - 1, celldim, astart + astep * (iin + j), bstart + celldim * bstep * iout + cellstep * j, nextzerofill);
						j += 1;
					}
					iout += 1;
					iin += stride(idim);
				}
			} else if (idim == 2) {
				val astep1 = adims(0)*adims(1);
				val bstep1 = bdims(1);
				val cellstep1 = inDims(0)*inDims(1);
//				var iin1 = -pad(2);
        var iin1 = 0;
				var iout1 = 0;
				while (iout1 < bdims(2)) {     
					var astart1 = astart + astep1 * iin1;
					var bstart1 = bstart + celldim * bstep1 * iout1;
          val zerofill1 = zerofill || iin1 < 0 || iin1 >= adims(2);
					var j1 = 0;
					while (j1 < iwidth) {     
						var	astep0 = adims(0);
						var	cellstep0 = inDims(0);
//						var iin0 = -pad(1);
            var iin0 = 0;
						var iout0 = 0;
						while (iout0 < bdims(1)) {
							val astart0 = astart1 + astep0 * iin0;
							val bstart0 = bstart1 + celldim * iout0;
							val zerofill0 = zerofill1 || iin0 < 0 || iin0 >= adims(1);
 //             if (zerofill0) {
							if (false) {
            	  var j0 = 0;
            	  while (j0 < cellstep1) {
            		  i2c.data(bstart0 + j0) = 0f;
            		  j0 += 1;
            	  }            
              } else {
            	  if (cellstep0 > FFilter.arraycopy) {
            		  System.arraycopy(a.data, astart0, i2c.data, bstart0, cellstep1);
            	  } else {
            		  var j0 = 0;
            		  while (j0 < cellstep1) {
                    if (bstart0+j0 >=i2c.data.length) {
                      throw new RuntimeException("i2c %d %d %d %d %d %d %d" format (bstart1, bstart0, j0, iout1, iout0, iin1, iin0))
                    }
                    if (astart0+j0 >=a.data.length) {
                      throw new RuntimeException("a %d %d %d %d %d" format (astart1, astart0, j0, iin1, iin0))
                    }
            			  i2c.data(bstart0 + j0) = a.data(astart0 + j0);
            			  j0 += 1;
            		  }
            	  }
              }
							iout0 += 1;
							iin0 += stride(1);
						}
						astart1 += astep1;
						bstart1 += cellstep1;
						j1 += 1;
					}
					iout1 += 1;
					iin1 += stride(2);
				}
			} else if (idim == 1) {
				var astep = adims(0);
				var cellstep = inDims(0);
				var iin = -pad(idim);
				var iout = 0;
				while (iout < bdims(idim)) {
					val astart0 = astart + astep * iin;
					val bstart0 = bstart + celldim * iout;
          val zerofill0 = zerofill || iin < 0 || iin >= adims(idim);
          if (zerofill0) {
        	  var j = 0;
        	  while (j < iwidth * cellstep) {
        		  i2c.data(bstart0 + j) = 0f;
        		  j += 1;
        	  }            
          } else {
        	  if (cellstep*iwidth > FFilter.arraycopy) {
        		  System.arraycopy(a.data, astart0, i2c.data, bstart0, cellstep*iwidth);
        	  } else {
        		  var j = 0;
        		  while (j < iwidth * cellstep) {
        			  i2c.data(bstart0 + j) = a.data(astart0 + j);
        			  j += 1;
        		  }
        	  }
          }
					iout += 1;
					iin += stride(idim);
				}    
			} else {
        if (zerofill) {
        	var i = 0;
        	while (i < iwidth) {
        		i2c.data(bstart + i) = 0f;
        		i += 1;
        	}         
        } else {
        	if (iwidth > FFilter.arraycopy) {
        		System.arraycopy(a.data, astart, i2c.data, bstart, iwidth);
        	} else {
        		var i = 0;
        		while (i < iwidth) {
			i2c.data(bstart + i) = a.data(astart + i);
        			i += 1;
        		}
        	}
        }
			}
			i2c
	}
*/

	def _col2im(a:FMat, b:FMat, i2c:FMat, idim:Int, celldim:Int, astart:Int, bstart:Int):FMat = {
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

	override def * (a:FMat):FMat = {
			convolve(a);
	}

	override def ^* (a:FMat):FMat = {
			convolveT(a);
	}
	
	override def xavier(scale:Float):FFilter = FFilter.xavier(this, scale);
	
	override def xavier:FFilter = FFilter.xavier(this, 1f);
	
	override def transpose(p:IMat):FFilter = {
	  new FFilter(inDims, outDims, stride, pad, outPad, dataDims(dataDims.length-1) \ dataDims(0->(dataDims.length-1)), _transpose(p).data);
	}
}


class FFiltPair2(val omat:Filter, val a:FMat) extends Pair(omat.asInstanceOf[Mat], a) {
	def *^ (b:FMat):FMat = {
			omat.asInstanceOf[FFilter].convolveM(a, b);
	}
}

object FFilter {

	var im2colThreshold = 0;
	var arraycopy = 16;

	def apply(indims:IMat, stride:IMat, pad:IMat):FFilter = {
			new FFilter(indims, iones(1, indims.length), stride, pad, null, indims(0,0->(indims.length-1)) \ 1, new Array[Float](indims.data.reduce(_*_)))
	}

	def apply(indims:IMat, outdims:IMat, stride:IMat, pad:IMat, dataDims:IMat):FFilter = {
			new FFilter(indims, outdims, stride, pad, null, indims(0,0->(indims.length-1)) \ outdims(0), new Array[Float]((indims dotr outdims).v))
	}
	
  def apply(g:GFilter):FFilter = {
          val outnd = FMat.newOrCheckFMat(g.dims, null, g.GUID, "FFilter".##);
          val out = new FFilter(g.inDims, g.outDims, g.stride, g.pad, g.outPad, g.dataDims, outnd.data);
          GMat.GPUtoCPUarraycopy(g.pdata, 0, out.data, 0, g.length, "FFilter apply");
          out.tensorFormat = g.tensorFormat;
          out.convType = g.convType;
          out.setGUID(MurmurHash3_x64_64(Array(g.GUID), "FFilter apply".##));
          out;
	}

	def FFilter1D(w:Int, nstride:Int, npad:Int, noutpad:Int):FFilter = {
		val inDims = irow(w);
		val outDims = irow(1);
		val stride = irow(nstride);
		val pad = irow(npad);
		val outPad = irow(noutpad);		
		new FFilter(inDims, outDims, stride, pad, outPad, irow(w,1), new Array[Float](w));
	}
	
	def FFilter1D(w:Int, nstride:Int, npad:Int):FFilter = FFilter1D(w, nstride, npad, 0);

	def FFilter1Dd(w:Int, din:Int, dout:Int, nstride:Int, npad:Int, noutpad:Int):FFilter = {
		val inDims = irow(din, w);
		val outDims = irow(dout, 1);
		val stride = irow(1, nstride);
		val pad = irow(0, npad);
		val outPad = irow(0, noutpad);
		new FFilter(inDims, outDims, stride, pad, outPad, irow(din, w, dout), new Array[Float](din*dout*w));
	}
	
	def FFilter1Dd(w:Int, din:Int, dout:Int, nstride:Int, npad:Int):FFilter = FFilter1Dd(w, din, dout, nstride, npad, 0);

	def FFilter2D(w:Int, h:Int, nstride:Int, npad:Int, noutpad:Int):FFilter = {
		val inDims = irow(w, h);
		val outDims = irow(1, 1);
		val stride = irow(nstride, nstride);
		val pad = irow(npad, npad);
		val outPad = irow(noutpad, noutpad);
		new FFilter(inDims, outDims, stride, pad, outPad, irow(w,h), new Array[Float](w*h));
	}
	
	def FFilter2D(w:Int, h:Int, nstride:Int, npad:Int):FFilter = FFilter2D(w, h, nstride, npad, 0);

	def FFilter2Dd(w:Int, h:Int, din:Int, dout:Int, nstride:Int, npad:Int, noutpad:Int):FFilter = {
		val inDims = irow(din, w, h);
		val outDims = irow(dout, 1, 1);
		val stride = irow(1, nstride, nstride);
		val pad = irow(0, npad, npad);
		val outPad = irow(0, noutpad, noutpad);		
		new FFilter(inDims, outDims, stride, pad, outPad, irow(din, w, h, dout), new Array[Float](din*dout*w*h));
	}
	
	def FFilter2Dd(w:Int, h:Int, din:Int, dout:Int, nstride:Int, npad:Int):FFilter = FFilter2Dd(w, h, din, dout, nstride, npad, 0);

	def FFilter2Ddn(w:Int, h:Int, din:Int, dout:Int, nstride:Int, npad:Int, noutpad:Int):FFilter = {
		val inDims = irow(din, w, h, 1);
		val outDims = irow(dout, 1, 1, 1);
		val stride = irow(1, nstride, nstride, 1);
		val pad = irow(0, npad, npad, 0);
		val outPad = irow(0, noutpad, noutpad, 0);
		new FFilter(inDims, outDims, stride, pad, outPad, irow(din, w, h, dout), new Array[Float](din*dout*w*h));
	}
	
	def FFilter2Ddn(w:Int, h:Int, din:Int, dout:Int, nstride:Int, npad:Int):FFilter = FFilter2Ddn(w, h, din, dout, nstride, npad, 0);
	
	def xavier(f:FFilter, fscale:Float):FFilter = {
	  val scale = f.inDims.data.reduce(_*_);
	  FFunctions.normrnd(0, fscale/math.sqrt(scale).toFloat, f);
	  f;
	}

}
