import BIDMat.{Mat, FMat, DMat, IMat, CMat, BMat, CSMat, SMat, SDMat, GMat, GIMat, GSMat, HMat}
import BIDMat.MatFunctions._
import BIDMat.SciFunctions._
import BIDMat.Solvers._
import BIDMat.Plotting._

{
var a = new Array[Int](1)
jcuda.runtime.JCuda.cudaGetDeviceCount(a)
if (a(0) > 0) jcuda.runtime.JCuda.initialize
Mat.hasCUDA = a(0)
}
printf("%d CUDA device%s found\n", Mat.hasCUDA, if (Mat.hasCUDA == 1) "" else "s")
