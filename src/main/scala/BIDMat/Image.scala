package BIDMat
import Mat._
import java.io._
import javax.imageio.ImageIO
import java.awt.image.BufferedImage 
import java.awt.Dimension
import java.awt.Graphics
import java.awt.Graphics2D
import java.awt.RenderingHints
import javax.swing._
import javax.swing.WindowConstants._

@SerialVersionUID(100L)
class Image(val img:BufferedImage) extends Serializable {
  
  var frame:JFrame = null;  
    
  final val width = img.getWidth
  
  final val height = img.getHeight
  
  def toIMat:IMat = {
    val mat = IMat(width,height)
    img.getRGB(0, 0, width, height, mat.data, 0, width)  // Actually returns packed ARGB values in the IMat
    mat
  }
  
  def toFMat:FMat = {
    val raster = img.getData                             // Seems to always hold RGB or ARGB data, independent of image encoding.
    val depth = raster.getNumDataElements 
    if (depth != 1) {
      throw new RuntimeException("Image is not grayscale, cannot push to an FMat")
    }
    val mat = FMat(width,height)
    raster.getPixels(0, 0, width, height, mat.data)
    mat
  }  
    
  def toFND:FND= {
    val height = img.getHeight
    val width = img.getWidth
    val raster = img.getData                            // Seems to always hold RGB or ARGB data, independent of image encoding.
    val depth = raster.getNumDataElements
    val mat = FND(depth, width,  height)
    raster.getPixels(0, 0, width, height, mat.data)
    mat
  }
  
  def redraw(mat:IMat):Image = {
    img.setRGB(0, 0, mat.nrows, mat.ncols, mat.data, 0, width);
    repaint
    this
  }
  
  def redraw(mat:FMat):Image = {
    img.getRaster.setPixels(0, 0, mat.nrows, mat.ncols, mat.data);
    repaint
    this
  } 
  
  def redraw(mat:FND):Image = {    
    val width = mat.dims(1);
    val height = mat.dim(2);
    img.getRaster.setPixels(0, 0, mat.dims(1), mat.dims(2), mat.data);
    repaint
    this;
  }
  
  def repaint = {
    if (frame != null) {
      frame.repaint();
    }
  }
  
  def resize(w0:Int, h0:Int):Image = {
    val w = if (w0 < 0) math.round(1f*width/height*h0).toInt else w0;
    val h = if (h0 < 0) math.round(1f*height/width*w0).toInt else h0;
    val im0 = img.getScaledInstance(w, h, java.awt.Image.SCALE_SMOOTH)
    val im = new BufferedImage(w, h, img.getType);
    val g:Graphics2D = im.createGraphics();
    g.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BILINEAR);
    g.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
    g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
    g.drawImage(im0, 0, 0, null);
    g.dispose; 
    new Image(im)    
  }
  
  def scale(factor:Double):Image = {
    resize(math.round(width*factor).toInt, math.round(height*factor).toInt)
  }
     
  def show(title0:String):BufferedImage = {
    if (Mat.inline) {
      img 
    } else {
    	val panel = new ImagePanel(img);
    	val title = if (title0 != null) title0 else "Image " + Image.imageCount;
    	frame = new JFrame(title);
    	Image.imageCount += 1;
    	frame.add(panel);
    	frame.pack;
    	frame.setDefaultCloseOperation(HIDE_ON_CLOSE);// or try DO_NOTHING_ON_CLOSE 
    	frame.setVisible(true);
    	Image.dummyImage.img
    }
  }
  
  def show():BufferedImage = show(null)
}

class ImagePanel(img:BufferedImage) extends JPanel                                                
{       
  setPreferredSize(new Dimension(img.getWidth, img.getHeight));
  
  override def paintComponent(g0:Graphics) =                                 
  {                                                                           
    if (null != img) {
      val scale = math.min(1.0*getWidth()/img.getWidth, 1.0*getHeight()/img.getHeight)
      val w = math.round(scale*img.getWidth()).toInt
      val h = math.round(scale*img.getHeight()).toInt
      val g = g0.asInstanceOf[Graphics2D]
      g.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BILINEAR);
      g.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
      g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
      g.drawImage(img, 0, 0, w, h, null)         
    }
  }                                                                           
} 


object Image {
  def loadImage(fname:String):Image = new Image(ImageIO.read(new File(fname)));
  
  private var imageCount = 0;
  
  var format = "jpg";
  
  val dummyImage = apply(MatFunctions.ones(1,1)*255f)
  
  def apply(mat:IMat):Image = {
    val width = mat.nrows
    val height = mat.ncols
    val img:BufferedImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB)
    img.setRGB(0, 0, width, height, mat.data, 0, width)
    new Image(img)
  }
  
  def apply(mat:FMat):Image = {
    val width = mat.nrows
    val height = mat.ncols
    val img:BufferedImage = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY)
    img.getRaster.setPixels(0, 0, width, height, mat.data)
    new Image(img)
  } 
  
  def apply(mat:FND):Image = { 	  
    val width = mat.dims(1)
    val height = mat.dim(2)
    val img:BufferedImage = mat.dims(0) match {
  	    case 1 => new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY)
  	    case 3 => new BufferedImage(width, height, BufferedImage.TYPE_3BYTE_BGR)
  	    case 4 => new BufferedImage(width, height, BufferedImage.TYPE_4BYTE_ABGR)
  	    case _ => throw new RuntimeException("Image from FND dimension not recognized")
  	  }
    img.getRaster.setPixels(0, 0, width, height, mat.data)
    new Image(img)
  }

}