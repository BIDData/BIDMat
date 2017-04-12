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
    
  def toFMat:FMat= {
    val height = img.getHeight;
    val width = img.getWidth;
    val mat = FMat.make(Array(4, width,  height));
    val ints = new Array[Int](height*width);
    img.getRGB(0, 0, width, height, ints, 0, width);      // Should be ARGB
    val mdata = mat.data;
    var i = 0;
    while (i < height*width) {                            // Produce an RGBA tensor
	  val ii =  ints(i);
	  mdata(i*4+2) = ii & 0xff;                           // B
	  mdata(i*4+1) = (ii >> 8) & 0xff;                    // G
	  mdata(i*4+0) = (ii >> 16) & 0xff;                   // R
	  mdata(i*4+3) = (ii >> 24) & 0xff;                   // A
	  i += 1;
    }
    mat
  }
  
  def redraw(mat:IMat):Image = {
    img.setRGB(0, 0, mat.nrows, mat.ncols, mat.data, 0, width);
    repaint
    this
  } 
  
  def redraw(mat:FMat):Image = {    
    val width = mat.dims(1);
    val height = mat.dims(2);
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
      val dd = mat.dims.length;
		  val width = if (dd > 2) mat.dims(1) else mat.dims(0);
		  val height = if (dd > 2) mat.dims(2) else mat.dims(1);
		  val ints = new Array[Int](width*height);
		  val mdata = mat.data;
		  val mult = 1+256+65536;
		  var i = 0;
		  val img:BufferedImage = 
		  		if (dd == 2 || mat.dims(0) == 1) {
		  			val im = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);
		  			while (i < height*width) {
		  				ints(i) = math.min(255,mdata(i).asInstanceOf[Int]) * mult;
		  				i += 1;
		  			}
		  			im;
		  		} else if (mat.dims(0) == 3) {
		  			val im = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		  			while (i < height*width) {
		  				ints(i) = ((((
		  				    (mdata(3*i).asInstanceOf[Int] & 0xff) << 8) +               // R
		  						(mdata(3*i+1).asInstanceOf[Int] & 0xff)) << 8) +             // G
		  						(mdata(3*i+2).asInstanceOf[Int] & 0xff));                      // B
		  				i += 1;
		  			}
		  			im; 
		  		} else if (mat.dims(0) == 4) {
		  			val im = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		  			while (i < height*width) {
		  				ints(i) = ((((((
		  				    (mdata(4*i+3).asInstanceOf[Int] & 0xff) << 8) +           // A
		  						(mdata(4*i+0).asInstanceOf[Int] & 0xff)) << 8) +           // R
		  						(mdata(4*i+1).asInstanceOf[Int] & 0xff)) << 8) +             // G
		  						(mdata(4*i+2).asInstanceOf[Int] & 0xff));                      // B
		  				i += 1;
		  			}
		  			im; 
		  		} else {
		  			throw new RuntimeException("Image from FMat dimension not recognized")
		  		}
		  img.setRGB(0, 0, width, height, ints, 0, width);
		  new Image(img);
  }

}
