package BIDMat
import ptolemy.plot._

object Plotting { 
  var ifigure:Int = 1
  
  def _plot(mats:Mat*)(xlog:Boolean=false, ylog:Boolean=false, isconnected:Boolean=true, bars:Boolean=false):Plot = {
    var p:Plot = new Plot
    p.setXLog(xlog)
    p.setYLog(ylog)
    p.setBars(bars)
    p.setConnected(isconnected)
    val dataset = 0
    if (mats.length == 1) {
      val m = mats(0)
      if (m.nrows == 1 || m.ncols == 1) { 
      	m match { 
      	case mf:FMat => for (i <- 0 until m.length) p.addPoint(dataset, i, mf(i), isconnected)
      	case md:DMat => for (i <- 0 until m.length) p.addPoint(dataset, i, md(i), isconnected)
      	case mi:IMat => for (i <- 0 until m.length) p.addPoint(dataset, i, mi(i), isconnected)
      	}
      } else {
      	for (i <- 0 until m.ncols) {
      		m match { 
      		case mf:FMat => for (j <- 0 until m.nrows) p.addPoint(i, j, mf(j,i), isconnected)
      		case md:DMat => for (j <- 0 until m.nrows) p.addPoint(i, j, md(j,i), isconnected)
      		case mi:IMat => for (j <- 0 until m.nrows) p.addPoint(i, j, mi(j,i), isconnected)
      		}     
      	}
      }
    } else {
      var i = 0
    	while (i*2 < mats.length) {      
    		(mats(2*i), mats(2*i+1)) match { 
    		case (a:FMat, b:FMat) => for (j <- 0 until a.length) p.addPoint(i, a(j), b(j), isconnected)
    		case (a:FMat, b:DMat) => for (j <- 0 until a.length) p.addPoint(i, a(j), b(j), isconnected)
    		case (a:DMat, b:FMat) => for (j <- 0 until a.length) p.addPoint(i, a(j), b(j), isconnected)
    		case (a:DMat, b:DMat) => for (j <- 0 until a.length) p.addPoint(i, a(j), b(j), isconnected)
    		case (a:FMat, b:IMat) => for (j <- 0 until a.length) p.addPoint(i, a(j), b(j), isconnected)
    		case (a:DMat, b:IMat) => for (j <- 0 until a.length) p.addPoint(i, a(j), b(j), isconnected)
    		case (a:IMat, b:FMat) => for (j <- 0 until a.length) p.addPoint(i, a(j), b(j), isconnected)
    		case (a:IMat, b:DMat) => for (j <- 0 until a.length) p.addPoint(i, a(j), b(j), isconnected)
    		case (a:IMat, b:IMat) => for (j <- 0 until a.length) p.addPoint(i, a(j), b(j), isconnected)
    		}  
    		i += 1
    	}
    }
    var pframe:PlotFrame = new PlotFrame("Figure "+ifigure, p)
    ifigure += 1
    pframe.setVisible(true)
    p
  }
  
  def plot(mats:Mat*) = _plot(mats: _*)()
  
  def loglog(mats:Mat*) = _plot(mats: _*)(xlog=true, ylog=true)
  
  def semilogx(mats:Mat*) = _plot(mats: _*)(xlog=true)
  
  def semilogy(mats:Mat*) = _plot(mats: _*)(ylog=true)

  def barplot(mats:Mat*) = _plot(mats: _*)(isconnected=false, bars=true)
  
  def barloglog(mats:Mat*) = _plot(mats: _*)(xlog=true, ylog=true, isconnected=false, bars=true)
  
  def barsemilogx(mats:Mat*) = _plot(mats: _*)(xlog=true, isconnected=false, bars=true)
  
  def barsemilogy(mats:Mat*) = _plot(mats: _*)(ylog=true, isconnected=false, bars=true)
  
  def p_plot(mats:Mat*) = _plot(mats: _*)(isconnected=false)
  
  def ploglog(mats:Mat*) = _plot(mats: _*)(xlog=true, ylog=true, isconnected=false)
  
  def psemilogx(mats:Mat*) = _plot(mats: _*)(xlog=true, isconnected=false)
  
  def psemilogy(mats:Mat*) = _plot(mats: _*)(ylog=true, isconnected=false)
   
  
  def hist(m:Mat, nbars:Int=10) = { 
    import SciFunctions._
    var p:Histogram = new Histogram
    val dataset = 0
    if (m.nrows == 1 || m.ncols == 1) { 
    	m match { 
    	case mf:FMat => {
    	  var vmax = maxi(mf,0).v
    	  var vmin = mini(mf,0).v
    	  p.setBinWidth((vmax-vmin)/nbars)
    	  for (i <- 0 until m.length) p.addPoint(dataset, mf(i))
    	}
    	case md:DMat => {
    		var vmax = maxi(md,0).v
    	  var vmin = mini(md,0).v
    	  p.setBinWidth((vmax-vmin)/nbars)
    	  for (i <- 0 until m.length) p.addPoint(dataset, md(i))
    	}
    	case mi:IMat => {
    		var vmax = maxi(mi,0).v.asInstanceOf[Double]
    	  var vmin = mini(mi,0).v
    	  p.setBinWidth((vmax-vmin)/nbars)
    	  for (i <- 0 until m.length) p.addPoint(dataset, mi(i))
    	}
      }
    }
    var pframe:PlotFrame = new PlotFrame("Figure "+ifigure, p)
    ifigure += 1
    pframe.setVisible(true)
  }
}
