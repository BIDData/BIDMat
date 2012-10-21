package BIDMat
import MatFunctions._
import SciFunctions._
import Learner._

trait Model {
  def initmodel(data:Mat, target:Mat):Mat
  def gradfun(data:Mat, target:Mat, model:Mat, diff:Mat):Double
}

trait Updater {
	def init(model:Mat):Unit
  def update(model:Mat, update:Mat, step:Int):Unit
}

case class Learner(datamat:Mat, targetmat:Mat, model:Model, updater:Updater) (
    blocksize:Int = 10000,
    npasses:Int = 1,
    memwindow:Double = 1000000,
    convwindow:Double = 10000000,
    convslope:Double = 1e-5,
    secprint:Double = 1
    ) 
{

  val modelmat = model.initmodel(datamat, targetmat)
  val (nrows, ncols) = size(modelmat)
  val diff:Mat = modelmat.zeros(nrows, ncols)
  val n = datamat.ncols
  val nw = blocksize/memwindow
  val nww = blocksize/convwindow
  updater.init(modelmat)
  
  def run() = {
  	var done:Boolean = false
  	var ipass = 0
  	var llest:Double = 0
  	var llder:Double = 0
  	var llold:Double = 0
  	var tsecs:Double = 0
  	tic
  	while (ipass < npasses && ! done) {
  		var i = 0
  		while (i < n && ! done) {
  			var iend = math.min(n, i+blocksize)
  			var dslice = datamat(?, i->iend)
  			var tslice = targetmat(?, i->iend)
  			val tll = model.gradfun(dslice, tslice, modelmat, diff)
  			updater.update(modelmat, diff, blocksize)
  			llest = (1/nw)*(tll + (nw-1)*llest)
  			llder = (1/nww)*(tll-llold + (nww-1)*llder)
  			llold = tll
  			i += blocksize
  			if (llder > 0 && llder < 0.5*convslope*blocksize*blocksize) {
  				done = true
  			}
  			if (toc > tsecs || done) {
  				println("pass="+ipass+" i="+i+" t="+toc+" secs, ll="+llest)
  				tsecs += secprint
  			}
  		}
  		ipass += 1
  	}
  }
}

abstract class RegressionModel extends Model {
    
  var tpred:FMat = null
  var fttarget:FMat = null
  var lls:FMat = null
  var gradw:FMat = null

  def regfn(targ:FMat, pred:FMat, lls:FMat, gradw:FMat):Unit
  
  def initmodel(data:Mat,target:Mat):Mat = {
    val m = size(data, 1)
    val n = size(target, 1)
    val out = 0.1f*normrnd(0,1,m,n)
    out
  }
  
  def gradfun(data:Mat, target:Mat, model:Mat, diff:Mat):Double = {
  	val sdata = data.asInstanceOf[SMat]
  	val ftarget = target.asInstanceOf[FMat]
  	val fmodel = model.asInstanceOf[FMat]
  	val fdiff = diff.asInstanceOf[FMat]
  	
  	fttarget = checkSize(fttarget, target.ncols, target.nrows)
  	lls = checkSize(lls, fttarget)
  	tpred = checkSize(tpred, fttarget)
  	gradw = checkSize(gradw, fttarget)
  	
  	fttarget ~ ftarget.t
  	tpred ~ sdata Tx fmodel 
  	regfn(fttarget, tpred, lls, gradw)
  	fdiff ~ sdata * gradw
    mean(mean(lls)).v
  }
}

class LogisticModel extends RegressionModel {

  var tfact:FMat = null
  var ptfact:FMat = null
  var epred:FMat = null
  var lle:FMat = null
  var tmp0:FMat = null
  var tmp1:FMat = null
      
  override def regfn(targ:FMat, pred:FMat, lls:FMat, gradw:FMat):Unit = {
    tfact = checkSize(tfact, targ)
    ptfact = checkSize(ptfact, targ)
    epred = checkSize(epred, targ)
    lle = checkSize(lle, targ)
    
    var i = 0
    while (i < targ.length) {
      tfact.data(i) = 1-2*targ.data(i)
      ptfact.data(i) = math.min(40f, pred.data(i) * tfact.data(i))
      i+=1
    }
    exp(ptfact, epred)
    log1p(epred, lle)
    lls ~ row(-1) * lle
    i = 0
    while (i < targ.length) {
      gradw.data(i) = - tfact.data(i) * epred.data(i) / (1 + epred.data(i))
      i+=1
    }
    Mat.nflops += 8L * targ.length
  }
  
  def logisticMap(targ:FMat, pred:FMat, lls:FMat, gradw:FMat) = {
    var i = 0
    while (i < targ.length) {
      val tfact:Double = (1-2*targ.data(i))
      val epred:Double = math.exp(math.min(40f, pred.data(i) * tfact))
      lls.data(i) = -math.log1p(epred).asInstanceOf[Float]
      gradw.data(i) = (- tfact * epred / (1 + epred)).asInstanceOf[Float]
      i += 1
    }
    Mat.nflops += 14L * targ.length
  }
    
  def logisticMap3(targ:FMat, pred:FMat, lls:FMat, gradw:FMat) = {
    tfact = checkSize(tfact, targ)
    epred = checkSize(epred, targ)
    tmp0 = checkSize(tmp0, targ)
    tmp1 = checkSize(tmp1, targ)
    
    tfact ~ row(1) - (tmp0 ~ row(2)*targ)
    min(40f, tmp0 ~ tfact *@ pred, tmp1)
    exp(tmp1, epred)
    log1p(epred, tmp0)
    lls ~ row(-1) * tmp0
    gradw ~ -1 * (tmp0 ~ tfact *@ (tmp1 ~ epred /@ (tmp0 ~ 1 + epred)))
  }
  
  def logisticMap4(targ:FMat, pred:FMat, lls:FMat, gradw:FMat) = {
    tfact = checkSize(tfact, targ)
    ptfact = checkSize(ptfact, targ)
    epred = checkSize(epred, targ)
    lle = checkSize(lle, targ)
    
    Learner.mapfun2x2((targ:Float, pred:Float) => (1-2*targ, math.min(40f, pred * (1-2*targ))), targ, pred, tfact, ptfact)
    exp(ptfact, epred)
    log1p(epred, lle)
    lls ~ row(-1) * lle
    Learner.mapfun2x1((tfact:Float, epred:Float)=>(-tfact*epred/(1+epred)), tfact, epred, gradw)
    Mat.nflops += 8L * targ.length
  }

}

class ADAGradUpdater(
    gradwindow:Float = 1000000,
    alpha:Float = 1e-12f,
    eps:Float = 1e-1f) extends Updater 
{
  var sumSq:FMat = null
  var tmp0:FMat = null
  var nsteps:Long = 1
  
  def init(model:Mat) = {
  	sumSq = zeros(size(model,1), size(model,2))
  	nsteps = 1
  }
  
	def update(model:Mat, update:Mat, step:Int):Unit = {
	  val fmodel = model.asInstanceOf[FMat]
	  val fupdate = update.asInstanceOf[FMat]
	  val nw = (step/gradwindow)
	  var i = 0 
	  var maxupdate:Float = 0
	  var maxsumsq:Float = 0
	  while (i < model.length) {
	    maxupdate = math.max(maxupdate, fupdate.data(i))
	    maxsumsq = math.max(maxsumsq, sumSq.data(i))
	    sumSq.data(i) = (1/nw)*(fupdate.data(i)*fupdate.data(i) + (nw-1)*sumSq.data(i))
	    fmodel.data(i) = fmodel.data(i) + alpha*step*fupdate.data(i)/fsqrt(math.max(eps, sumSq.data(i)))
	    i += 1
	  }
	  println("maxs="+maxupdate+" "+maxsumsq)
	  nsteps += 1
	}
}

object Learner {
  
  def checkSize(a:FMat, nr:Int, nc:Int):FMat = {
    if (a.asInstanceOf[AnyRef] != null && a.nrows == nr && a.ncols == nc) {
      a
    } else {
      zeros(nr, nc)
    }
  }
  
  def checkSize(a:FMat, b:FMat):FMat = checkSize(a, b.nrows, b.ncols)
  
  def fsqrt(v:Float):Float = math.sqrt(v).asInstanceOf[Float]
  
  def mapfun2x2(fn:(Float, Float)=>(Float, Float), in0:FMat, in1:FMat, out0:FMat, out1:FMat) = {
    if (in0.nrows != in1.nrows || in0.nrows != out0.nrows || in0.nrows != out1.nrows ||
        in0.ncols != in1.ncols || in0.ncols != out0.ncols || in0.ncols != out1.ncols) {
      throw new RuntimeException("dimensions mismatch")
    }
    var i = 0
    while (i < in0.length) {
      val (v1, v2) = fn(in0.data(i), in1.data(i))
      out0.data(i) = v1
      out1.data(i) = v2
      i += 1
    }
  }
  def mapfun2x1(fn:(Float, Float)=>Float, in0:FMat, in1:FMat, out0:FMat) = {
    if (in0.nrows != in1.nrows || in0.nrows != out0.nrows ||
        in0.ncols != in1.ncols || in0.ncols != out0.ncols) {
      throw new RuntimeException("dimensions mismatch")
    }
    var i = 0
    while (i < in0.length) {
      out0.data(i) = fn(in0.data(i), in1.data(i))
      i += 1
    }
  }
}
