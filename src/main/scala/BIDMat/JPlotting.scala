package BIDMat


object JPlotting {
	import org.jfree.chart._
	import org.jfree.chart.plot._
	import org.jfree.data.xy._
	import org.jfree.data.statistics._
	import org.jfree.chart.renderer.xy._
	import org.jfree.util._
	import java.awt.image.BufferedImage
	
  var ifigure:Int = 0;
  val marksmat = Array("points","dots","various");
  
  var plotXscale = 500;
  var plotYscale = 400;
  var plotPscale = 0.2f;

  def _addToDataset(dataset:XYSeriesCollection, mats:Array[Mat]) = {
    if (mats.length == 1) {
      val m = mats(0)
      if (m.nrows == 1 || m.ncols == 1) { 
      	val p = new XYSeries("plot");
      	m match { 
      	case mf:FMat => for (i <- 0 until m.length) p.add(i, mf(i));
      	case md:DMat => for (i <- 0 until m.length) p.add(i, md(i));
      	case mi:IMat => for (i <- 0 until m.length) p.add(i, mi(i));
      	}
      	dataset.addSeries(p);
      } else {
      	for (i <- 0 until m.ncols) {
      		val p = new XYSeries("plot %d" format i);
      		m match { 
      		case mf:FMat => for (j <- 0 until m.nrows) p.add(j, mf(j,i));
      		case md:DMat => for (j <- 0 until m.nrows) p.add(j, md(j,i));
      		case mi:IMat => for (j <- 0 until m.nrows) p.add(j, mi(j,i));
      		}
      		dataset.addSeries(p);
      	}
      }
    } else {
      var i = 0
    	while (i*2 < mats.length) { 
    		val p = new XYSeries("plot %d" format i);
    		(mats(2*i), mats(2*i+1)) match { 
    		case (a:FMat, b:FMat) => for (j <- 0 until a.length) p.add(a(j), b(j))
    		case (a:FMat, b:DMat) => for (j <- 0 until a.length) p.add(a(j), b(j))
    		case (a:DMat, b:FMat) => for (j <- 0 until a.length) p.add(a(j), b(j))
    		case (a:DMat, b:DMat) => for (j <- 0 until a.length) p.add(a(j), b(j))
    		case (a:FMat, b:IMat) => for (j <- 0 until a.length) p.add(a(j), b(j))
    		case (a:DMat, b:IMat) => for (j <- 0 until a.length) p.add(a(j), b(j))
    		case (a:IMat, b:FMat) => for (j <- 0 until a.length) p.add(a(j), b(j))
    		case (a:IMat, b:DMat) => for (j <- 0 until a.length) p.add(a(j), b(j))
    		case (a:IMat, b:IMat) => for (j <- 0 until a.length) p.add(a(j), b(j))
    		}  
    		dataset.addSeries(p);
    		i += 1;
    	}
    }
  }
  
  def _plot(mats:Mat*)(xlog:Boolean=false, ylog:Boolean=false, isconnected:Boolean=true, bars:Boolean=false, marks:Int = 0):BufferedImage = {
   	import java.awt.Color;
   	if (Mat.inline) {
   		System.setProperty("java.awt.headless", "true");
   	}
   	val fmats = mats.toArray.map(MatFunctions.cpu)
    val dataset = new XYSeriesCollection();    
    _addToDataset(dataset, fmats);
    ifigure += 1;
    val chart = if (isconnected) {
      ChartFactory.createXYLineChart("Figure %d" format ifigure, "X", "Y", dataset);
    } else if (bars) {
      ChartFactory.createXYBarChart("Figure %d" format ifigure, "X", false, "Y", dataset);
    } else {
      ChartFactory.createScatterPlot("Figure %d" format ifigure, "X", "Y", dataset);
    }
    
    val plot = chart.getXYPlot;
    if (!isconnected && !bars) {
      val dot = ShapeUtilities.createDiamond(plotPscale);
      val renderer = plot.getRenderer;
      for (i <- 0 until dataset.getSeriesCount) renderer.setSeriesShape(i, dot);
    }
    if (xlog) {
    	val xaxis = new org.jfree.chart.axis.LogarithmicAxis("X");
    	plot.setDomainAxis(xaxis);
    }
    if (ylog) {
    	val yaxis = new org.jfree.chart.axis.LogarithmicAxis("Y");
    	plot.setRangeAxis(yaxis);
    }
    plot.setBackgroundPaint(Color.white);
    plot.setDomainGridlinePaint(Color.gray);
    plot.setRangeGridlinePaint(Color.gray);
    if (!Mat.inline) {
    	val frame = new ChartFrame("Figure %d" format ifigure, chart);
    	frame.pack();
    	frame.setVisible(true);
    }
    chart.createBufferedImage(plotXscale, plotYscale);
  }
  
  def plot(mats:Mat*) = _plot(mats: _*)();
  
  def scatter(mats:Mat*) = _plot(mats: _*)(marks=1, isconnected=false);
  
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
   
  
  def hist(m:Mat, nbars:Int=10):BufferedImage = { 
  	import java.awt.Color;
  	if (Mat.inline) {
  		System.setProperty("java.awt.headless", "true");
  	}
    var dd = new HistogramDataset();
    if (m.nrows == 1 || m.ncols == 1) { 
    	m match { 
    	case mf:FMat => {
    	  dd.addSeries("H1", DMat(mf).data, nbars);
    	}
    	case md:DMat => {
    	  dd.addSeries("H1", md.data, nbars);
    	}
    	case mi:IMat => {
    	  dd.addSeries("H1", DMat(mi).data, nbars);
    	}
      }
    }
    val show = false;
    val toolTips = false;
    val urls = true;
    ifigure += 1;
    val chart = ChartFactory.createHistogram("Figure %d" format ifigure, "X", "Count", dd, PlotOrientation.VERTICAL, show, toolTips, urls);
    val plot = chart.getXYPlot();

    plot.getRenderer.asInstanceOf[XYBarRenderer].setBarPainter(new StandardXYBarPainter);
    plot.setBackgroundPaint(Color.white);
    if (!Mat.inline) {
    	val frame = new ChartFrame("Figure %d" format ifigure, chart);
    	frame.pack();
    	frame.setVisible(true);
    }
    chart.createBufferedImage(plotXscale, plotYscale);
  }
}
