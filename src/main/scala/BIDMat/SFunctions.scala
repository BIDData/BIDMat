package BIDMat

object SFunctions {
  /** min, max, sum, prod, maxi, mini for SMats with no output matrix*/
  import GMat.BinOp._
  
  def min(a:SMat, b:SMat) =  a.ssMatOpS(b, FMat.minFun, op_min, null)
  def max(a:SMat, b:SMat) = a.ssMatOpS(b, FMat.maxFun, op_max, null)  
  
  def sum(a:SMat, n:Int, omat:Mat):FMat = {
    a match {
//      case aa:GSMat => aa.sum(n, omat);
      case _ => a.ssReduceOp(n, FMat.idFun, FMat.sumFun, omat);
    }
  }
  def sum(a:SMat, n:Int):FMat = sum(a, n, null);
  def sum(a:SMat, omat:Mat):FMat = sum(a, 0, omat);
  def sum(a:SMat):FMat = sum(a, 0, null);
   
  def maxi(a:SMat, n:Int) = a.ssReduceOp(n, FMat.idFun, FMat.maxFun, null);
  def mini(a:SMat, n:Int) = a.ssReduceOp(n, FMat.idFun, FMat.minFun, null);
  def amax(a:SMat, n:Int) = a.ssReduceOp(n, FMat.idFun, FMat.maxFun, null);
  def amin(a:SMat, n:Int) = a.ssReduceOp(n, FMat.idFun, FMat.minFun, null);
  
  def maxi(a:SMat) = a.ssReduceOp(0, FMat.idFun, FMat.maxFun, null)
  def mini(a:SMat) = a.ssReduceOp(0, FMat.idFun, FMat.minFun, null)
  def amax(a:SMat) = a.ssReduceOp(0, FMat.idFun, FMat.maxFun, null)
  def amin(a:SMat) = a.ssReduceOp(0, FMat.idFun, FMat.minFun, null)
  
  def min(a:SMat, b:Float, omat:Mat):SMat = {
    a match {
//      case aa:GSMat => aa.GSDop(GMat(b), omat, op_min);
      case _ => a.ssMatOpScalar(b, FMat.minFun, omat)
    }
  }
  
  def min(b:Float, a:SMat, omat:Mat):SMat = min(a, b, omat);
  def min(a:SMat, b:Float):SMat = min(a, b, null);
  def min(b:Float, a:SMat):SMat = min(a, b, null);
  
  def max(a:SMat, b:Float, omat:Mat):SMat = {
    a match {
//      case aa:GSMat => aa.GSDop(GMat(b), omat, op_max);
      case _ => a.ssMatOpScalar(b, FMat.maxFun, omat)
    }
  }
  
  def max(b:Float, a:SMat, omat:Mat):SMat = max(a, b, omat);
  def max(a:SMat, b:Float):SMat = max(a, b, null);
  def max(b:Float, a:SMat):SMat = max(a, b, null);

  
  def countnz(a:SMat, omat:Mat) = a.countnz(0, omat)
  def countnz(a:SMat, n:Int, omat:Mat) = a.countnz(n, omat)
  def countnz(a:SMat) = a.countnz(0, null)
  def countnz(a:SMat, n:Int) = a.countnz(n, null)
  
  /** min, max, sum, prod, cumsum, maxi, mini for SMats with output matrix*/

  def maxi(a:SMat, n:Int, omat:Mat) = a.ssReduceOp(n, FMat.idFun, FMat.maxFun, omat)
  def mini(a:SMat, n:Int, omat:Mat) = a.ssReduceOp(n, FMat.idFun, FMat.minFun, omat)
  def maxi(a:SMat, omat:Mat) = a.ssReduceOp(0, FMat.idFun, FMat.maxFun, omat)
  def mini(a:SMat, omat:Mat) = a.ssReduceOp(0, FMat.idFun, FMat.minFun, omat)
  def amax(a:SMat, omat:Mat) = a.ssReduceOp(0, FMat.idFun, FMat.maxFun, omat)
  def amin(a:SMat, omat:Mat) = a.ssReduceOp(0, FMat.idFun, FMat.minFun, omat)
  
}
  