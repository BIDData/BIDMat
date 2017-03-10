package BIDMat

object SDFunctions {
  /** min, max, sum, prod, maxi, mini for SMats with no output matrix*/
   /** min, max, sum, maxi, mini for SDMats with no output matrix*/
  
  def min(a:SDMat, b:SDMat) = a.ssMatOp(b, DMat.minFun, null)
  def max(a:SDMat, b:SDMat) = a.ssMatOp(b, DMat.maxFun, null)
  def sum(a:SDMat, n:Int) = a.ssReduceOp(n, DMat.idFun, DMat.sumFun, null)
  def maxi(a:SDMat, n:Int) = a.ssReduceOp(n, DMat.idFun, DMat.maxFun, null)
  def mini(a:SDMat, n:Int) = a.ssReduceOp(n, DMat.idFun, DMat.minFun, null)
  def amax(a:SDMat, n:Int) = a.ssReduceOp(n, DMat.idFun, DMat.maxFun, null)
  def amin(a:SDMat, n:Int) = a.ssReduceOp(n, DMat.idFun, DMat.minFun, null)
  
  def sum(a:SDMat) = a.ssReduceOp(0, DMat.idFun, DMat.sumFun, null)
  def maxi(a:SDMat) = a.ssReduceOp(0, DMat.idFun, DMat.maxFun, null)
  def mini(a:SDMat) = a.ssReduceOp(0, DMat.idFun, DMat.minFun, null)
  def amax(a:SDMat) = a.ssReduceOp(0, DMat.idFun, DMat.maxFun, null)
  def amin(a:SDMat) = a.ssReduceOp(0, DMat.idFun, DMat.minFun, null)
  
  /** min, max, sum, maxi, mini for SDMats with output matrix*/
  
  def sum(a:SDMat, n:Int, omat:Mat) = a.ssReduceOp(n, DMat.idFun, DMat.sumFun, omat)
  def maxi(a:SDMat, n:Int, omat:Mat) = a.ssReduceOp(n, DMat.idFun, DMat.maxFun, null)
  def mini(a:SDMat, n:Int, omat:Mat) = a.ssReduceOp(n, DMat.idFun, DMat.minFun, null)
  def amax(a:SDMat, n:Int, omat:Mat) = a.ssReduceOp(n, DMat.idFun, DMat.maxFun, null)
  def amin(a:SDMat, n:Int, omat:Mat) = a.ssReduceOp(n, DMat.idFun, DMat.minFun, null)
  
  def sum(a:SDMat, omat:Mat) = a.ssReduceOp(0, DMat.idFun, DMat.sumFun, omat)
  def maxi(a:SDMat, omat:Mat) = a.ssReduceOp(0, DMat.idFun, DMat.maxFun, omat)
  def mini(a:SDMat, omat:Mat) = a.ssReduceOp(0, DMat.idFun, DMat.minFun, omat)
  def amax(a:SDMat, omat:Mat) = a.ssReduceOp(0, DMat.idFun, DMat.maxFun, omat)
  def amin(a:SDMat, omat:Mat) = a.ssReduceOp(0, DMat.idFun, DMat.minFun, omat)
  
  def countnz(a:SDMat, omat:Mat) = a.countnz(0, omat)
  def countnz(a:SDMat, n:Int, omat:Mat) = a.countnz(n, omat)
  def countnz(a:SDMat) = a.countnz(0, null)
  def countnz(a:SDMat, n:Int) = a.countnz(n, null)
  
  
  
  def min(a:SDMat, b:Double, omat:Mat) = a.ssMatOpScalar(b, DMat.minFun, omat)
  def max(a:SDMat, b:Double, omat:Mat) = a.ssMatOpScalar(b, DMat.maxFun, omat)
  def min(b:Double, a:SDMat, omat:Mat) = a.ssMatOpScalar(b, DMat.minFun, omat)
  def max(b:Double, a:SDMat, omat:Mat) = a.ssMatOpScalar(b, DMat.maxFun, omat)
  def min(a:SDMat, b:Double) = a.ssMatOpScalar(b, DMat.minFun, null)
  def max(a:SDMat, b:Double) = a.ssMatOpScalar(b, DMat.maxFun, null)
  def min(b:Double, a:SDMat) = a.ssMatOpScalar(b, DMat.minFun, null)
  def max(b:Double, a:SDMat) = a.ssMatOpScalar(b, DMat.maxFun, null)
  
  def sum(a:CMat, n:Int) = a.ccReduceOpv(n, CMat.vecAddFun, null)
  def sum(a:CMat, n:Int, c:Mat) = a.ccReduceOpv(n, CMat.vecAddFun, c)
     
  def countnz(a:Mat, n:Int):IMat = countnz(a, n, null)
  def countnz(a:Mat):IMat = countnz(a, 0, null)
  def countnz(a:Mat, n:Int, omat:Mat):IMat = {
    a match {
      case as:SMat => as.countnz(n, omat)
      case as:SDMat => as.countnz(n, omat)
    }
  }
}
  