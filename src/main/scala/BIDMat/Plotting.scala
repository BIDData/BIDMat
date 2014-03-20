package BIDMat
import ptolemy.plot._
import java.awt._
import java.awt.geom.AffineTransform
import java.awt.image.BufferedImage
import javax.swing._
import javax.imageio.stream.FileImageOutputStream
import javax.imageio.ImageIO
import java.io._

object Plotting { 
  var ifigure:Int = 1
  
  def _plot(mats:Mat*)(xlog:Boolean=false, ylog:Boolean=false, isconnected:Boolean=true):Plot = {
    var p:Plot = new Plot
    p.setXLog(xlog)
    p.setYLog(ylog)
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

  def heatmap(m:Mat) = {
    val hc:HeatChart = new HeatChart(m)
    val img:BufferedImage = hc.getChartImage(true)
    hc.saveToFile(new File("heat_map_"+ifigure+".jpg"))
    val jl:JLabel = new JLabel(new ImageIcon(img))
    val jp:JPanel = new JPanel
    jp.add(jl)
    val jsp:JScrollPane = new JScrollPane(jp)
    val jFrame:JFrame = new JFrame("Figure "+ifigure)        
    jFrame.getContentPane().add( jsp )
    jFrame.setSize(800, 600)
    ifigure += 1
    jFrame.setVisible(true)
  }
}

class HeatChart(mat:Mat) {
  private var xValues:Array[Double] = new Array[Double](mat.ncols)
  private var yValues:Array[Double] = new Array[Double](mat.nrows)
  setXValues(0, 1)
  setYValues(0, 1)
  private var xValuesHorizontal:Boolean = false
  private var yValuesHorizontal:Boolean = true


  private var cellSize:Dimension = new Dimension(10,10)
  private var margin:Int = 20
  private var backgroundColor = Color.WHITE
  
  private var chartSize:Dimension = new Dimension(100,100)
  
  private var highValueColor:Color = Color.BLUE
  private var lowValueColor:Color = Color.WHITE

  private var colorValueDistance:Int = 1
  private var colorScale:Double = 1.0
  
  private var heatMapSize:Dimension = new Dimension(1000,1000)  

  private var heatMapTL:Point = new Point(0, 0)
  private var heatMapBR:Point = new Point(800, 800)
  private var heatMapC:Point = new Point(400, 400)

  private var axisThickness:Int = 2
  private var axisColor:Color = Color.BLACK
  private var axisLabelsFont:Font = new Font("Sans-Serif", Font.PLAIN, 12)
  private var axisLabelColor:Color = Color.BLACK
  private var xAxisLabel:String = "X Label"
  private var yAxisLabel:String = "Y Label"
  private var axisValuesColor:Color = Color.BLACK
  private var axisValuesFont:Font = new Font("Sans-Serif", Font.PLAIN, 10)
  private var xAxisValuesFrequency:Int = 1
  private var yAxisValuesFrequency:Int = 1
  private var showXAxisValues:Boolean = true
  private var showYAxisValues:Boolean = true

  private var xAxisValuesHeight:Int = 0
  private var xAxisValuesWidthMax:Int = 0
	
  private var yAxisValuesHeight:Int = 0
  private var yAxisValuesAscent:Int = 0
  private var yAxisValuesWidthMax:Int = 0
	
  private var xAxisLabelSize:Dimension = new Dimension(0,0)
  private var xAxisLabelDescent:Int = 0
	
  private var yAxisLabelSize:Dimension = new Dimension(0,0)
  private var yAxisLabelAscent:Int = 0


  private var lowValue:Double = min(mat)
  private var highValue:Double = max(mat)


  updateColorDistance()


  //should be replaced by built-in min and max to speed up
  private def min(mat:Mat):Double = {
    mat match {
      case mi:IMat => IMin(mi).toDouble
      case mf:FMat => FMin(mf).toDouble
      case md:DMat => DMin(md)
    }
  }


  private def max(mat:Mat):Double = {
    mat match {
      case mi:IMat => IMax(mi).toDouble
      case mf:FMat => FMax(mf).toDouble
      case md:DMat => DMax(md)
    }
  }

  private def IMin(mat:IMat):Int = {
    var minV:Int = mat(0)
    for (i:Int <- 0 until mat.length) {
      if (mat(i) < minV)
        minV = mat(i)
    }
    minV
  }

  private def FMin(mat:FMat):Float = {
    var minV:Float = mat(0)
    for (i:Int <- 0 until mat.length) {
      if (mat(i) < minV)
        minV = mat(i)
    }
    minV
  }

  private def DMin(mat:DMat):Double = {
    var minV:Double = mat(0)
    for (i:Int <- 0 until mat.length) {
      if (mat(i) < minV)
        minV = mat(i)
    }
    minV
  }

  private def IMax(mat:IMat):Int = {
    var maxV:Int = mat(0)
    for (i:Int <- 0 until mat.length) {
      if (mat(i) > maxV)
        maxV = mat(i)
    }
    maxV
  }

  private def FMax(mat:FMat):Float = {
    var maxV:Float = mat(0)
    for (i:Int <- 0 until mat.length) {
      if (mat(i) > maxV)
        maxV = mat(i)
    }
    maxV
  }

  private def DMax(mat:DMat):Double = {
    var maxV:Double = mat(0)
    for (i:Int <- 0 until mat.length) {
      if (mat(i) > maxV)
        maxV = mat(i)
    }
    maxV
  }

  def getChartImage(alpha:Boolean):BufferedImage = {
    measureComponents()
    updateCoordinates()
    var imageType:Int = if ( alpha ) BufferedImage.TYPE_4BYTE_ABGR else BufferedImage.TYPE_3BYTE_BGR
		
    var chartImage:BufferedImage = new BufferedImage(chartSize.width, chartSize.height, imageType)

    var chartGraphics:Graphics2D = chartImage.createGraphics()
		
    chartGraphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
    RenderingHints.VALUE_ANTIALIAS_ON)
		
    chartGraphics.setColor(backgroundColor)
    chartGraphics.fillRect(0, 0, chartSize.width, chartSize.height)
    
    drawHeatMap(chartGraphics, mat)

    drawXLabel(chartGraphics)
    drawYLabel(chartGraphics)

    drawAxisBars(chartGraphics)

    drawXValues(chartGraphics)
    drawYValues(chartGraphics)
    chartImage
  }

  private def drawHeatMap(chartGraphics:Graphics2D, mat:Mat) = {

    var noYCells:Int = mat.nrows
    var noXCells:Int = mat.ncols
		
    var heatMapImage:BufferedImage = new BufferedImage(heatMapSize.width, heatMapSize.height, BufferedImage.TYPE_INT_ARGB)
    var heatMapGraphics:Graphics2D = heatMapImage.createGraphics()
  
  mat match {
    case mi:IMat => 
      for ( x <- 0 until noXCells ) {
        for ( y <- 0 until noYCells ) {
          heatMapGraphics.setColor( getCellColor( mi(y+x*noYCells).toDouble, lowValue, highValue ) )

          var cellX:Int = x*cellSize.width
	  var cellY:Int = y*cellSize.height

	  heatMapGraphics.fillRect(cellX, cellY, cellSize.width, cellSize.height)
        }
      }
    case mf:FMat =>
      for ( x <- 0 until noXCells ) {
        for ( y <- 0 until noYCells ) {
          heatMapGraphics.setColor( getCellColor( mf(y+x*noYCells).toDouble, lowValue, highValue ) )

          var cellX:Int = x*cellSize.width
	  var cellY:Int = y*cellSize.height

	  heatMapGraphics.fillRect(cellX, cellY, cellSize.width, cellSize.height)
        }
      }
    case md:DMat =>
      for ( x <- 0 until noXCells ) {
        for ( y <- 0 until noYCells ) {
          heatMapGraphics.setColor( getCellColor( md(y+x*noYCells), lowValue, highValue ) )

          var cellX:Int = x*cellSize.width
	  var cellY:Int = y*cellSize.height

	  heatMapGraphics.fillRect(cellX, cellY, cellSize.width, cellSize.height)
        }
      }    
  }

    chartGraphics.drawImage(heatMapImage, heatMapTL.x, heatMapTL.y, heatMapSize.width, heatMapSize.height, null)
  }



  private def getCellColor(data:Double, min:Double, max:Double):Color = {
    var range:Double = max - min
    var position:Double = data - min

    var percentPosition:Double = position / range

    var colorPosition:Int = getColorPosition(percentPosition)
    var r:Int = lowValueColor.getRed()
    var g:Int = lowValueColor.getGreen
    var b:Int = lowValueColor.getBlue()

    for ( i <- 0 until colorPosition ) {
      var rDistance:Int = r - highValueColor.getRed()
      var gDistance:Int = g - highValueColor.getGreen()
      var bDistance:Int = b - highValueColor.getBlue()
			
      if ( (math.abs(rDistance) >= math.abs(gDistance) ) 
        && (math.abs(rDistance) >= math.abs(bDistance) ) ) {
        r = changeColorValue(r, rDistance)
      } else if (math.abs(gDistance) >= math.abs(bDistance)) {
        g = changeColorValue(g, gDistance)
      } else {
        b = changeColorValue(b, bDistance)
      }
    }	
    new Color(r, g, b) 
  }


  private def getColorPosition(percentPosition:Double):Int = {
    math.round( colorValueDistance * math.pow(percentPosition, colorScale) ).toInt
  }

  private def updateColorDistance() = {
    var r1:Int = lowValueColor.getRed()
    var g1:Int = lowValueColor.getGreen()
    var b1:Int = lowValueColor.getBlue()
    var r2:Int = highValueColor.getRed()
    var g2:Int = highValueColor.getGreen()
    var b2:Int = highValueColor.getBlue()
		
    colorValueDistance = math.abs(r1 - r2)
    colorValueDistance += math.abs(g1 - g2)
    colorValueDistance += math.abs(b1 - b2)
  }

  private def changeColorValue(colorValue:Int, colorDistance:Int):Int = {
    if (colorDistance < 0) {
      colorValue+1
    } else if (colorDistance > 0) {
      colorValue-1
    } else {
      colorValue
    }
  }

  private def measureComponents() =  {
    var chartImage:BufferedImage = new BufferedImage(1, 1, BufferedImage.TYPE_INT_ARGB)
    var tempGraphics:Graphics2D = chartImage.createGraphics()
		

    if (xAxisLabel != null) {
      tempGraphics.setFont(axisLabelsFont)
      var metrics:FontMetrics = tempGraphics.getFontMetrics()
      xAxisLabelSize = new Dimension(metrics.stringWidth(xAxisLabel), metrics.getHeight())
      xAxisLabelDescent = metrics.getDescent()
    } else {
      xAxisLabelSize = new Dimension(0, 0)
    }

	
    if (yAxisLabel != null) {
      tempGraphics.setFont(axisLabelsFont)
      var metrics:FontMetrics = tempGraphics.getFontMetrics()
      yAxisLabelSize = new Dimension(metrics.stringWidth(yAxisLabel), metrics.getHeight())
      yAxisLabelAscent = metrics.getAscent()
    } else {
      yAxisLabelSize = new Dimension(0, 0)
    }
		

    if (showXAxisValues) {
      tempGraphics.setFont(axisValuesFont)
      var metrics:FontMetrics = tempGraphics.getFontMetrics()
      xAxisValuesHeight = metrics.getHeight()
      xAxisValuesWidthMax = 0
      
      for (i <- 0 until xValues.length) {
        var w:Int = metrics.stringWidth(xValues(i).toString())
          if (w > xAxisValuesWidthMax) {
            xAxisValuesWidthMax = w
          }
      }
    } else {
      xAxisValuesHeight = 0
    }

    if (showYAxisValues) {
      tempGraphics.setFont(axisValuesFont)
      var metrics:FontMetrics = tempGraphics.getFontMetrics()
      yAxisValuesHeight = metrics.getHeight()
      yAxisValuesAscent = metrics.getAscent()
      yAxisValuesWidthMax = 0
      
      for (i <-0 until yValues.length) {
        var w:Int = metrics.stringWidth(yValues(i).toString())
        if (w > yAxisValuesWidthMax) {
          yAxisValuesWidthMax = w
        }
      }
    } else {
      yAxisValuesHeight = 0
    }
      
    
    var heatMapWidth:Int = ( mat.ncols * cellSize.width)
    var heatMapHeight:Int = ( mat.nrows * cellSize.height)
    heatMapSize = new Dimension(heatMapWidth, heatMapHeight)

    var yValuesHorizontalSize:Int = 0
		
    if (yValuesHorizontal) {
      yValuesHorizontalSize = yAxisValuesWidthMax
    } else {
      yValuesHorizontalSize = yAxisValuesHeight
    }

    var xValuesVerticalSize:Int = 0
    if (xValuesHorizontal) {
      xValuesVerticalSize = xAxisValuesHeight
    } else {
      xValuesVerticalSize = xAxisValuesWidthMax
    }
		
    var chartWidth:Int = heatMapWidth + (2 * margin) + yAxisLabelSize.height + yValuesHorizontalSize + axisThickness
    var chartHeight:Int = heatMapHeight + (2 * margin) + xAxisLabelSize.height + xValuesVerticalSize + axisThickness
    chartSize = new Dimension(chartWidth, chartHeight)
  }

  private def updateCoordinates() {
    var x:Int = margin + axisThickness + yAxisLabelSize.height
    if (yValuesHorizontal) x+=yAxisValuesWidthMax else x+=yAxisValuesHeight
    var y:Int = margin
    heatMapTL = new Point(x, y)

    x = heatMapTL.x + heatMapSize.width
    y = heatMapTL.y + heatMapSize.height
    heatMapBR = new Point(x, y)
		
    x = heatMapTL.x + (heatMapSize.width / 2)
    y = heatMapTL.y + (heatMapSize.height / 2)
    heatMapC = new Point(x, y)
  }

  private def drawXLabel(chartGraphics:Graphics2D) = {
    if (xAxisLabel != null) {
      var yPosXAxisLabel:Int = chartSize.height - (margin / 2) - xAxisLabelDescent
      var xPosXAxisLabel:Int = heatMapC.x - (xAxisLabelSize.width / 2)

      chartGraphics.setFont(axisLabelsFont)
      chartGraphics.setColor(axisLabelColor)
      chartGraphics.drawString(xAxisLabel, xPosXAxisLabel, yPosXAxisLabel)
    }
  }
  
  private def drawYLabel(chartGraphics:Graphics2D) = {
    if (yAxisLabel != null) {
      var yPosYAxisLabel:Int = heatMapC.y + (yAxisLabelSize.width / 2)
      var xPosYAxisLabel:Int = (margin / 2) + yAxisLabelAscent
			
      chartGraphics.setFont(axisLabelsFont)
      chartGraphics.setColor(axisLabelColor)
			
      var transform:AffineTransform = chartGraphics.getTransform()
      var originalTransform:AffineTransform = transform.clone().asInstanceOf[AffineTransform]
      transform.rotate(math.toRadians(270), xPosYAxisLabel, yPosYAxisLabel)
      chartGraphics.setTransform(transform)
			
      chartGraphics.drawString(yAxisLabel, xPosYAxisLabel, yPosYAxisLabel)
			
      chartGraphics.setTransform(originalTransform)
    }
  }


  private def drawAxisBars(chartGraphics:Graphics2D) = {
    if (axisThickness > 0) {
      chartGraphics.setColor(axisColor)
			
      var x:Int = heatMapTL.x - axisThickness
      var y:Int = heatMapBR.y
      var width:Int = heatMapSize.width + axisThickness
      var height:Int = axisThickness
      chartGraphics.fillRect(x, y, width, height)
			
      x = heatMapTL.x - axisThickness
      y = heatMapTL.y
      width = axisThickness
      height = heatMapSize.height
      chartGraphics.fillRect(x, y, width, height)
    }
  }


  private def drawXValues(chartGraphics:Graphics2D) = {
    if (showXAxisValues) { 
		
      chartGraphics.setColor(axisValuesColor)
		
      for (i <- 0 until mat.ncols) {
        if (i % xAxisValuesFrequency == 0) {		
          var xValueStr:String = xValues(i).toString()		
          chartGraphics.setFont(axisValuesFont)
          var metrics:FontMetrics = chartGraphics.getFontMetrics()		
          var valueWidth:Int = metrics.stringWidth(xValueStr)
			
          if (xValuesHorizontal) {
            var valueXPos:Int = (i * cellSize.width) + ((cellSize.width / 2) - (valueWidth / 2))
            valueXPos += heatMapTL.x
            var valueYPos:Int = heatMapBR.y + metrics.getAscent() + 1
				
            chartGraphics.drawString(xValueStr, valueXPos, valueYPos)
          } else {
            var valueXPos:Int = heatMapTL.x + (i * cellSize.width) + ((cellSize.width / 2) + (xAxisValuesHeight / 2))
            var valueYPos:Int = heatMapBR.y + axisThickness + valueWidth
				
            var transform:AffineTransform = chartGraphics.getTransform()
            var originalTransform:AffineTransform = transform.clone().asInstanceOf[AffineTransform]
            transform.rotate(math.toRadians(270), valueXPos, valueYPos)
            chartGraphics.setTransform(transform)
				
            chartGraphics.drawString(xValueStr, valueXPos, valueYPos)
				
            chartGraphics.setTransform(originalTransform)
          }  
        }
      }
    }
  }


  private def drawYValues(chartGraphics:Graphics2D) = {
    if (showYAxisValues) {
      chartGraphics.setColor(axisValuesColor)	
      for (i <- 0 until mat.nrows ) {
        if (i % yAxisValuesFrequency == 0) {		
          var yValueStr:String = yValues(i).toString()		
          chartGraphics.setFont(axisValuesFont)
          var metrics:FontMetrics = chartGraphics.getFontMetrics()		
          var valueWidth:Int = metrics.stringWidth(yValueStr)		
          if (yValuesHorizontal) {
            var valueXPos:Int = margin + yAxisLabelSize.height + (yAxisValuesWidthMax - valueWidth)
            var valueYPos:Int = heatMapTL.y + (i * cellSize.height) + (cellSize.height/2) + (yAxisValuesAscent/2)
				
            chartGraphics.drawString(yValueStr, valueXPos, valueYPos)
          } else {
            var valueXPos:Int = margin + yAxisLabelSize.height + yAxisValuesAscent
            var valueYPos:Int = heatMapTL.y + (i * cellSize.height) + (cellSize.height/2) + (valueWidth/2)
				
            var transform:AffineTransform = chartGraphics.getTransform()
            var originalTransform:AffineTransform = transform.clone().asInstanceOf[AffineTransform]
            transform.rotate(math.toRadians(270), valueXPos, valueYPos)
            chartGraphics.setTransform(transform)
				
            chartGraphics.drawString(yValueStr, valueXPos, valueYPos)
				
            chartGraphics.setTransform(originalTransform)
          }
        }
      }
    }
  }


  def setXValues(xOffset:Double, xInterval:Double) = {		
    for (i <- 0 until mat.ncols) {
      xValues(i) = xOffset + (i * xInterval)
    }
  }

  def setYValues(yOffset:Double, yInterval:Double) = {		
    for (i <- 0 until mat.nrows) {
      yValues(i) = yOffset + (i * yInterval)
    }
  }


  def saveToFile(outputFile:File) = {
    var filename:String = outputFile.getName()
    var extPoint:Int = filename.lastIndexOf('.')

    if (extPoint < 0) {
      throw new IOException("Illegal filename: need a extension.")
    }

    var ext:String = filename.substring(extPoint + 1)
		
    if (ext.toLowerCase().equals("jpg") || ext.toLowerCase().equals("jpeg")) {
      var chart:BufferedImage = getChartImage(false)

      ImageIO.write(chart, ext, outputFile)
    } else {
      var chart:BufferedImage = getChartImage(true)
			
      ImageIO.write(chart, ext, outputFile)
    }
  }

}


