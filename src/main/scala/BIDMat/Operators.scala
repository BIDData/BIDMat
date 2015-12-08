package BIDMat
import MatFunctions._

/*
 * Type coercions for operators
 */ 

trait Mop {
  def myname:String;
  /*
   * Operations on matching types. These need to be overridden by the actual operator.
   */
  def op(a:FMat, b:FMat, c:Mat):FMat = {notImplemented(myname, a, b); a}
  def op(a:DMat, b:DMat, c:Mat):DMat = {notImplemented(myname, a, b); a}
  def op(a:IMat, b:IMat, c:Mat):IMat = {notImplemented(myname, a, b); a}
  def op(a:LMat, b:LMat, c:Mat):LMat = {notImplemented(myname, a, b); a}
  def op(a:CMat, b:CMat, c:Mat):CMat = {notImplemented(myname, a, b); a}
  def op(a:GMat, b:GMat, c:Mat):GMat = {notImplemented(myname, a, b); a}
  def op(a:GMat, b:GSMat, c:Mat):GMat = {notImplemented(myname, a, b); a}
  def op(a:GDMat, b:GDMat, c:Mat):GDMat = {notImplemented(myname, a, b); a}
  def op(a:GDMat, b:GSDMat, c:Mat):GDMat = {notImplemented(myname, a, b); a}
  def op(a:GIMat, b:GIMat, c:Mat):GIMat = {notImplemented(myname, a, b); a}
  def op(a:GLMat, b:GLMat, c:Mat):GLMat = {notImplemented(myname, a, b); a}

  def op(a:SMat, b:SMat, c:Mat):SMat = {notImplemented(myname, a, b); a}
  def op(a:SDMat, b:SDMat, c:Mat):SDMat = {notImplemented(myname, a, b); a}
  def op(a:GSMat, b:GSMat, c:Mat):GSMat = {notImplemented(myname, a, b); a}
  def op(a:GSDMat, b:GSDMat, c:Mat):GSDMat = {notImplemented(myname, a, b); a}
  def op(a:GSMat, b:GMat, c:Mat):GMat = {notImplemented(myname, a, b); b}
  def op(a:GSDMat, b:GDMat, c:Mat):GDMat = {notImplemented(myname, a, b); b}

  def sop(a:SMat, b:FMat, c:Mat):SMat = {notImplemented(myname, a, b); a}
  def sop(a:SDMat, b:DMat, c:Mat):SDMat = {notImplemented(myname, a, b); a}
  def sop(a:GSMat, b:GMat, c:Mat):GSMat = {notImplemented(myname, a, b); a}
  def sop(a:GSDMat, b:GDMat, c:Mat):GSDMat = {notImplemented(myname, a, b); a}
  
  /*
   * Default coercions when first mat is an FMat. These can be overridden. 
   */
  def op(a:FMat, b:SMat, c:Mat):FMat = op(a, full(b), c)
  def op(a:FMat, b:DMat, c:Mat):DMat = op(DMat(a), b, c)
  def op(a:FMat, b:IMat, c:Mat):FMat = op(a, FMat(b), c)  
  def op(a:FMat, b:LMat, c:Mat):FMat = op(a, FMat(b), c)
  def op(a:FMat, b:CMat, c:Mat):CMat = op(CMat(a), b, c)  
  def op(a:FMat, b:GMat, c:Mat):GMat = op(GMat(a), b, c)

  def op(a:FMat, b:Mat, c:Mat):Mat = {
    b match {
      case bb:FMat => op(a, bb, c)
      case bb:SMat => op(a, bb, c)
      case bb:DMat => op(a, bb, c)
      case bb:IMat => op(a, bb, c)
      case bb:LMat => op(a, bb, c)
      case bb:CMat => op(a, bb, c)
      case bb:GMat => op(a, bb, c)
    }
  } 
  
  /*
   * Default coercions when first arg is a DMat. 
   */
  def op(a:DMat, b:FMat, c:Mat):DMat = op(a, DMat(b), c)
  def op(a:DMat, b:IMat, c:Mat):DMat = op(a, DMat(b), c)
  def op(a:DMat, b:LMat, c:Mat):DMat = op(a, DMat(b), c)
  def op(a:DMat, b:SDMat, c:Mat):DMat = op(a, DMat(b), c)
  def op(a:DMat, b:CMat, c:Mat):CMat = op(CMat(a), b, c)
  def op(a:DMat, b:GMat, c:Mat):GMat = op(GMat(a), b, c)
 
  def op(a:DMat, b:Mat, c:Mat):Mat = {
    b match {
      case bb:FMat => op(a, bb, c)
      case bb:DMat => op(a, bb, c)
      case bb:IMat => op(a, bb, c)
      case bb:LMat => op(a, bb, c)
      case bb:SDMat => op(a, bb, c)
      case bb:CMat => op(a, bb, c)
      case bb:GMat => op(a, bb, c)
    }
  }  
  
  /*
   * Default coercions when first arg is an IMat
   */
  def op(a:IMat, b:FMat, c:Mat):FMat = op(FMat(a), b, c)
  def op(a:IMat, b:DMat, c:Mat):DMat = op(DMat(a), b, c)
  def op(a:IMat, b:LMat, c:Mat):LMat = op(LMat(a), b, c)
  def op(a:IMat, b:SMat, c:Mat):FMat = op(FMat(a), b, c)
  def op(a:IMat, b:CMat, c:Mat):CMat = op(CMat(a), b, c)
  def op(a:IMat, b:GMat, c:Mat):GMat = op(GMat(a), b, c)
  
  def op(a:IMat, b:Mat, c:Mat):Mat = {
    b match {
      case bb:FMat => op(a, bb, c)
      case bb:DMat => op(a, bb, c)
      case bb:IMat => op(a, bb, c)
      case bb:LMat => op(a, bb, c)
      case bb:SMat => op(a, bb, c)
      case bb:CMat => op(a, bb, c)
      case bb:GMat => op(a, bb, c)
    }
  }
  
    /*
   * Default coercions when first arg is an LMat
   */
  def op(a:LMat, b:FMat, c:Mat):FMat = op(FMat(a), b, c)
  def op(a:LMat, b:DMat, c:Mat):DMat = op(DMat(a), b, c)
  def op(a:LMat, b:IMat, c:Mat):LMat = op(a, LMat(b), c)
  def op(a:LMat, b:SMat, c:Mat):FMat = op(FMat(a), b, c)
  def op(a:LMat, b:CMat, c:Mat):CMat = op(CMat(a), b, c)
//  def op(a:LMat, b:GMat, c:Mat):GMat = op(GMat(a), b, c)
  
  def op(a:LMat, b:Mat, c:Mat):Mat = {
    b match {
      case bb:FMat => op(a, bb, c)
      case bb:DMat => op(a, bb, c)
      case bb:IMat => op(a, bb, c)
      case bb:LMat => op(a, bb, c)
      case bb:SMat => op(a, bb, c)
      case bb:CMat => op(a, bb, c)
    }
  }
  
  /*
   * Default coercions when first arg is a CMat
   */
  def op(a:CMat, b:FMat, c:Mat):CMat = op(a, CMat(b), c)
  def op(a:CMat, b:DMat, c:Mat):CMat = op(a, CMat(b), c)
  def op(a:CMat, b:IMat, c:Mat):CMat = op(a, CMat(b), c)
  
  def op(a:CMat, b:Mat, c:Mat):Mat = {
    b match {
      case bb:FMat => op(a, bb, c)
      case bb:DMat => op(a, bb, c)
      case bb:IMat => op(a, bb, c)
      case bb:CMat => op(a, bb, c)
    }
  }
  
  /*
   * Default coercions when first arg is an SMat
   */
  def op(a:SMat, b:FMat, c:Mat):FMat = op(full(a), b, c)
  
  def op(a:SMat, b:Mat, c:Mat):Mat = {
    b match {
      case bb:FMat => sop(a, bb, c)
      case bb:SMat => op(a, bb, c)
    }
  }
  
  def sop(a:SMat, b:Mat, c:Mat):Mat = {
    b match {
      case bb:FMat => sop(a, bb, c)
      case bb:SMat => op(a, bb, c)
    }
  }
  
  /*
   * Default coercions when first arg is an SDMat
   */
  def op(a:SDMat, b:DMat, c:Mat):DMat = op(full(a), b, c)

  
  def op(a:SDMat, b:Mat, c:Mat):Mat = {
    b match {
      case bb:DMat => op(a, bb, c)
      case bb:SDMat => sop(a, bb, c)
    }
  }
  
  def sop(a:SDMat, b:Mat, c:Mat):Mat = {
    b match {
      case bb:DMat => sop(a, bb, c)
      case bb:SDMat => op(a, bb, c)
    }
  }
  
  /*
   * Default coercions when first arg is a GMat
   */
  def op(a:GMat, b:FMat, c:Mat):GMat = op(a, GMat(b), c)
  def op(a:GMat, b:IMat, c:Mat):GMat = op(a, GMat(b), c)
  def op(a:GMat, b:DMat, c:Mat):GMat = op(a, GMat(b), c)
//  def op(a:GMat, b:GSMat, c:Mat):GMat = op(a, b, c)
  
  
  def op(a:GMat, b:Mat, c:Mat):Mat = {
    b match {
      case bb:FMat => op(a, bb, c)
      case bb:IMat => op(a, bb, c)
      case bb:DMat => op(a, bb, c)
      case bb:GMat => op(a, bb, c)
      case bb:GSMat => op(a, bb, c)
    }
  }      
  
  /*
   * Default coercions when first arg is a GIMat
   */
  def op(a:GIMat, b:FMat, c:Mat):GMat = op(GMat(a), GMat(b), c)
  def op(a:GIMat, b:IMat, c:Mat):GIMat = op(a, GIMat(b), c)
  def op(a:GIMat, b:DMat, c:Mat):GDMat = op(GDMat(a), GDMat(b), c)
  def op(a:GIMat, b:GMat, c:Mat):GMat = op(GMat(a), b, c)
  
  
  def op(a:GIMat, b:Mat, c:Mat):Mat = {
    b match {
      case bb:FMat => op(a, bb, c)
      case bb:IMat => op(a, bb, c)
      case bb:DMat => op(a, bb, c)
      case bb:GMat => op(a, bb, c)
      case bb:GIMat => op(a, bb, c)
    }
  }      
  
    /*
   * Default coercions when first arg is a GLMat
   */
  def op(a:GLMat, b:FMat, c:Mat):GDMat = op(GDMat(a), GDMat(b), c)
  def op(a:GLMat, b:IMat, c:Mat):GLMat = op(a, GLMat(LMat(b)), c)
  def op(a:GLMat, b:LMat, c:Mat):GDMat = op(GDMat(a), GDMat(b), c)
  def op(a:GLMat, b:GMat, c:Mat):GDMat = op(GDMat(a), GDMat(b), c)
  def op(a:GLMat, b:GIMat, c:Mat):GLMat = op(a, GLMat(b), c)
  
  
  def op(a:GLMat, b:Mat, c:Mat):Mat = {
    b match {
      case bb:FMat => op(a, bb, c)
      case bb:IMat => op(a, bb, c)
      case bb:LMat => op(a, bb, c)
      case bb:GMat => op(a, bb, c)
      case bb:GIMat => op(a, bb, c)
    }
  }
  /*
   * Default coercions when first arg is a GDMat
   */
  def op(a:GDMat, b:FMat, c:Mat):GDMat = op(a, GDMat(b), c)
  def op(a:GDMat, b:IMat, c:Mat):GDMat = op(a, GDMat(b), c)
  def op(a:GDMat, b:DMat, c:Mat):GDMat = op(a, GDMat(b), c)
//  def op(a:GDMat, b:GSDMat, c:Mat):GDMat = op(a, b, c)
  
  def op(a:GDMat, b:Mat, c:Mat):Mat = {
    b match {
      case bb:FMat => op(a, bb, c)
      case bb:IMat => op(a, bb, c)
      case bb:DMat => op(a, bb, c)
      case bb:GMat => op(a, bb, c)
      case bb:GDMat => op(a, bb, c)
      case bb:GSMat => op(a, bb, c)
      case bb:GSDMat => op(a, bb, c)
    }
  }
  /*
   * Default coercions when first arg is a GSMat
   */
//  def op(a:GSMat, b:GMat, c:Mat):GMat = op(a, b, c)
  def op(a:GSMat, b:FMat, c:Mat):GMat = op(a, GMat(b), c)
  def op(a:GSMat, b:IMat, c:Mat):GMat = op(a, GMat(b), c)
  def op(a:GSMat, b:DMat, c:Mat):GMat = op(a, GMat(b), c)
  
  
  def op(a:GSMat, b:Mat, c:Mat):Mat = {
    b match {
      case bb:FMat => op(a, bb, c)
      case bb:IMat => op(a, bb, c)
      case bb:DMat => op(a, bb, c)
      case bb:GMat => op(a, bb, c)
      case bb:GIMat => op(a, bb, c)
      case bb:GSMat => op(a, bb, c)
    }
  }
  
  def sop(a:GSMat, b:Mat, c:Mat):Mat = {
    b match {
      case bb:GMat => sop(a, bb, c)
      case bb:GSMat => sop(a, bb, c)
    }
  }
  
  /*
   * Default coercions when first arg is a GSDMat
   */
//  def op(a:GSDMat, b:GDMat, c:Mat):GDMat = op(a, b, c)
  def op(a:GSDMat, b:FMat, c:Mat):GDMat = op(a, GDMat(b), c)
  def op(a:GSDMat, b:IMat, c:Mat):GDMat = op(a, GDMat(b), c)
  def op(a:GSDMat, b:DMat, c:Mat):GDMat = op(a, GDMat(b), c)
  
  def op(a:GSDMat, b:Mat, c:Mat):Mat = {
    b match {
      case bb:FMat => op(a, bb, c)
      case bb:IMat => op(a, bb, c)
      case bb:DMat => op(a, bb, c)
      case bb:GMat => op(a, bb, c)
      case bb:GDMat => op(a, bb, c)
      case bb:GIMat => op(a, bb, c)
      case bb:GSMat => op(a, bb, c)
    }
  }
  
  def sop(a:GSDMat, b:Mat, c:Mat):Mat = {
    b match {
      case bb:GDMat => sop(a, bb, c)
      case bb:GSDMat => sop(a, bb, c)
    }
  }
  
  /*
   * No clue
   */
  def op(a:Mat, b:Mat, c:Mat):Mat = {
    a match {
      case aa:FMat => op(aa, b, c)
      case aa:DMat => op(aa, b, c)
      case aa:IMat => op(aa, b, c)
      case aa:CMat => op(aa, b, c)
      case aa:SMat => op(aa, b, c)
      case aa:SDMat => op(aa, b, c)
      case aa:GMat => op(aa, b, c)
      case aa:GIMat => op(aa, b, c)
      case aa:GSMat => op(aa, b, c)
//      case aa:GIMat => op(aa, b, c)
    }
  }
 
  def notImplemented(s:String, m:Mat):Mat = { 
    throw new RuntimeException("operator "+s+" not implemented for "+m.mytype)
  }
  def notImplemented(s:String, m:Mat, n:Mat):Mat = { 
    throw new RuntimeException("operator "+s+" not implemented for "+m.mytype+" and "+n.mytype)
  }
    
  /*
   * Right now these just create a new pair. 
   */
  def getFPair(c:Mat, a:FMat):FPair = new FPair(c, a) 
  def getDPair(c:Mat, a:DMat):DPair = new DPair(c, a)    
  def getIPair(c:Mat, a:IMat):IPair = new IPair(c, a)  
  def getLPair(c:Mat, a:LMat):LPair = new LPair(c, a)
  def getCPair(c:Mat, a:CMat):CPair = new CPair(c, a)  
  def getSPair(c:Mat, a:SMat):SPair = new SPair(c, a) 
  def getSDPair(c:Mat, a:SDMat):SDPair = new SDPair(c, a)
  def getGPair(c:Mat, a:GMat):GPair = new GPair(c, a) 
  def getGDPair(c:Mat, a:GDMat):GDPair = new GDPair(c, a)
  def getGIPair(c:Mat, a:GIMat):GIPair = new GIPair(c, a)
  def getGLPair(c:Mat, a:GLMat):GLPair = new GLPair(c, a)
  def getGSPair(c:Mat, a:GSMat):GSPair = new GSPair(c, a)
  def getGSDPair(c:Mat, a:GSDMat):GSDPair = new GSDPair(c, a)
}

object Mop_Plus extends Mop { // Plus includes only operations on matching types
  override def myname = "+"
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) + b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) + b
  override def op(a:IMat, b:IMat, c:Mat):IMat = getIPair(c, a) + b
  override def op(a:LMat, b:LMat, c:Mat):LMat = getLPair(c, a) + b
  override def op(a:CMat, b:CMat, c:Mat):CMat = getCPair(c, a) + b
  override def op(a:GMat, b:GMat, c:Mat):GMat = getGPair(c, a) + b
  override def op(a:GDMat, b:GDMat, c:Mat):GDMat = getGDPair(c, a) + b
  override def op(a:GIMat, b:GIMat, c:Mat):GIMat = getGIPair(c, a) + b
  override def op(a:GLMat, b:GLMat, c:Mat):GLMat = getGLPair(c, a) + b
  override def op(a:SMat, b:SMat, c:Mat):SMat = getSPair(c, a) + b
  override def sop(a:SMat, b:FMat, c:Mat):SMat = getSPair(c, a) + b
  override def sop(a:GSMat, b:GMat, c:Mat):GSMat = getGSPair(c, a) + b
  override def sop(a:GSDMat, b:GDMat, c:Mat):GSDMat = getGSDPair(c, a) + b
}

object Mop_Minus extends Mop { 
  override def myname = "-"
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) - b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) - b
  override def op(a:IMat, b:IMat, c:Mat):IMat = getIPair(c, a) - b
  override def op(a:LMat, b:LMat, c:Mat):LMat = getLPair(c, a) - b
  override def op(a:CMat, b:CMat, c:Mat):CMat = getCPair(c, a) - b
  override def op(a:GMat, b:GMat, c:Mat):GMat = getGPair(c, a) - b
  override def op(a:GIMat, b:GIMat, c:Mat):GIMat = getGIPair(c, a) - b
  override def op(a:GDMat, b:GDMat, c:Mat):GDMat = getGDPair(c, a) - b
  override def op(a:GLMat, b:GLMat, c:Mat):GLMat = getGLPair(c, a) - b
  override def op(a:SMat, b:SMat, c:Mat):SMat = getSPair(c, a) - b
  override def sop(a:SMat, b:FMat, c:Mat):SMat = getSPair(c, a) - b
  override def sop(a:GSMat, b:GMat, c:Mat):GSMat = getGSPair(c, a) - b
  override def sop(a:GSDMat, b:GDMat, c:Mat):GSDMat = getGSDPair(c, a) - b
}

object Mop_Times extends Mop { // Times includes dense-sparse combinations
  override def myname = "*"
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) * b
  override def op(a:FMat, b:SMat, c:Mat):FMat = getFPair(c, a) * b
  override def op(a:SMat, b:FMat, c:Mat):FMat = getSPair(c, a) * b
  override def op(a:SMat, b:SMat, c:Mat):SMat = getSPair(c, a) * b 
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) * b
  override def op(a:DMat, b:SDMat, c:Mat):DMat = getDPair(c, a) * b
  override def op(a:SDMat, b:DMat, c:Mat):DMat = getSDPair(c, a) * b
  override def op(a:SDMat, b:SDMat, c:Mat):SDMat = getSDPair(c, a) * b
  override def op(a:IMat, b:IMat, c:Mat):IMat = getIPair(c, a) * b
  override def op(a:LMat, b:LMat, c:Mat):LMat = getLPair(c, a) * b
  override def op(a:CMat, b:CMat, c:Mat):CMat = getCPair(c, a) * b
  override def op(a:GMat, b:GMat, c:Mat):GMat = getGPair(c, a) * b
  override def op(a:GDMat, b:GDMat, c:Mat):GDMat = getGDPair(c, a) * b
  override def op(a:GMat, b:GSMat, c:Mat):GMat = getGPair(c, a) * b
  override def op(a:GDMat, b:GSDMat, c:Mat):GDMat = getGDPair(c, a) * b
  override def op(a:GSMat, b:GMat, c:Mat):GMat = getGSPair(c, a) * b
  override def op(a:GSDMat, b:GDMat, c:Mat):GDMat = getGSDPair(c, a) * b
}

object Mop_TimesT extends Mop { 
  override def myname = "*^"
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) *^ b
  override def op(a:FMat, b:SMat, c:Mat):FMat = getFPair(c, a) *^ b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) *^ b
  override def op(a:DMat, b:SDMat, c:Mat):DMat = getDPair(c, a) *^ b
  override def op(a:GMat, b:GMat, c:Mat):GMat = getGPair(c, a) *^ b
  override def op(a:GDMat, b:GDMat, c:Mat):GDMat = getGDPair(c, a) *^ b
  override def op(a:GMat, b:GSMat, c:Mat):GMat = getGPair(c, a) *^ b
  override def op(a:GDMat, b:GSDMat, c:Mat):GDMat = getGDPair(c, a) *^ b
}

object Mop_TTimes extends Mop { 
  override def myname = "^*"
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) Tx b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) Tx b
  override def op(a:SMat, b:FMat, c:Mat):FMat = getSPair(c, a) Tx b
  override def op(a:SDMat, b:DMat, c:Mat):DMat = getSDPair(c, a) Tx b
  override def op(a:GMat, b:GMat, c:Mat):GMat = getGPair(c, a) Tx b
  override def op(a:GDMat, b:GDMat, c:Mat):GDMat = getGDPair(c, a) Tx b
  override def op(a:GSMat, b:GMat, c:Mat):GMat = getGSPair(c, a) Tx b
  override def op(a:GSDMat, b:GDMat, c:Mat):GDMat = getGSDPair(c, a) Tx b
}

object Mop_Div extends Mop { 
	override def myname = "/<"
  override def op(a:FMat, b:FMat, c:Mat):FMat = a /< b
  override def op(a:DMat, b:DMat, c:Mat):DMat = a /< b
  override def op(a:CMat, b:CMat, c:Mat):CMat = a /< b
}

object Mop_RSolve extends Mop { 
  override def myname = "\\"
  override def op(a:FMat, b:FMat, c:Mat):FMat = a \\ b
  override def op(a:DMat, b:DMat, c:Mat):DMat = a \\ b
  override def op(a:CMat, b:CMat, c:Mat):CMat = a \\ b
}

object Mop_ETimes extends Mop { 
  override def myname = "*@"
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) *@ b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) *@ b
  override def op(a:IMat, b:IMat, c:Mat):IMat = getIPair(c, a) *@ b
  override def op(a:LMat, b:LMat, c:Mat):LMat = getLPair(c, a) *@ b
  override def op(a:CMat, b:CMat, c:Mat):CMat = getCPair(c, a) *@ b
  override def op(a:GMat, b:GMat, c:Mat):GMat = getGPair(c, a) *@ b
  override def op(a:GDMat, b:GDMat, c:Mat):GDMat = getGDPair(c, a) *@ b
  override def op(a:GIMat, b:GIMat, c:Mat):GIMat = getGIPair(c, a) *@ b
  override def op(a:GLMat, b:GLMat, c:Mat):GLMat = getGLPair(c, a) *@ b
  override def op(a:SMat, b:SMat, c:Mat):SMat = getSPair(c, a) *@ b
  override def sop(a:SMat, b:FMat, c:Mat):SMat = getSPair(c, a) *@ b
  override def sop(a:GSMat, b:GMat, c:Mat):GSMat = getGSPair(c, a) *@ b
  override def sop(a:GSDMat, b:GDMat, c:Mat):GSDMat = getGSDPair(c, a) *@ b
}

object Mop_Kron extends Mop { 
  override def myname = "**"
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) ** b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) ** b
  override def op(a:IMat, b:IMat, c:Mat):IMat = getIPair(c, a) ** b
  override def op(a:LMat, b:LMat, c:Mat):LMat = getLPair(c, a) ** b
}

object Mop_EDiv extends Mop {
  override def myname = "/"
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) / b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) / b
  override def op(a:IMat, b:IMat, c:Mat):IMat = getIPair(c, a) / b
  override def op(a:LMat, b:LMat, c:Mat):LMat = getLPair(c, a) / b
  override def op(a:CMat, b:CMat, c:Mat):CMat = getCPair(c, a) / b
  override def op(a:GMat, b:GMat, c:Mat):GMat = getGPair(c, a) / b
  override def op(a:GDMat, b:GDMat, c:Mat):GDMat = getGDPair(c, a) / b
  override def op(a:GIMat, b:GIMat, c:Mat):GIMat = getGIPair(c, a) / b
  override def op(a:GLMat, b:GLMat, c:Mat):GLMat = getGLPair(c, a) / b
  override def op(a:SMat, b:SMat, c:Mat):SMat = getSPair(c, a) / b
  override def sop(a:SMat, b:FMat, c:Mat):SMat = getSPair(c, a) / b
  override def sop(a:GSMat, b:GMat, c:Mat):GSMat = getGSPair(c, a) / b
  override def sop(a:GSDMat, b:GDMat, c:Mat):GSDMat = getGSDPair(c, a) / b
}

object Mop_Dot extends Mop { 
  override def myname = "dot"
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) dot b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) dot b
  override def op(a:IMat, b:IMat, c:Mat):IMat = getIPair(c, a) dot b
  override def op(a:LMat, b:LMat, c:Mat):LMat = getLPair(c, a) dot b
  override def op(a:CMat, b:CMat, c:Mat):CMat = getCPair(c, a) dot b
  override def op(a:GMat, b:GMat, c:Mat):GMat = getGPair(c, a) dot b
  override def op(a:GDMat, b:GDMat, c:Mat):GDMat = getGDPair(c, a) dot b
}

object Mop_Dotr extends Mop { 
  override def myname = "dotr"
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) dotr b
  override def op(a:IMat, b:IMat, c:Mat):IMat = getIPair(c, a) dotr b
  override def op(a:LMat, b:LMat, c:Mat):LMat = getLPair(c, a) dotr b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) dotr b
  override def op(a:GMat, b:GMat, c:Mat):GMat = getGPair(c, a) dotr b
  override def op(a:GDMat, b:GDMat, c:Mat):GDMat = getGDPair(c, a) dotr b
}

object Mop_Pow extends Mop {
  override def myname = "^"
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) ^ b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) ^ b
  override def op(a:GMat, b:GMat, c:Mat):GMat = getGPair(c, a) ^ b
  override def op(a:GDMat, b:GDMat, c:Mat):GDMat = getGDPair(c, a) ^ b
}

object Mop_HCat extends Mop { 
	override def myname = "\\"
  override def op(a:FMat, b:FMat, c:Mat):FMat = a \ b
  override def op(a:DMat, b:DMat, c:Mat):DMat = a \ b
  override def op(a:IMat, b:IMat, c:Mat):IMat = a \ b
  override def op(a:LMat, b:LMat, c:Mat):LMat = a \ b
  override def op(a:SMat, b:SMat, c:Mat):SMat = a \ b
  override def op(a:SDMat, b:SDMat, c:Mat):SDMat = a \ b
  override def op(a:CMat, b:CMat, c:Mat):CMat = a \ b
  override def op(a:GMat, b:GMat, c:Mat):GMat = getGPair(c, a) \ b
  override def op(a:GDMat, b:GDMat, c:Mat):GDMat = getGDPair(c, a) \ b
  override def op(a:GIMat, b:GIMat, c:Mat):GIMat = getGIPair(c, a) \ b
  override def op(a:GLMat, b:GLMat, c:Mat):GLMat = getGLPair(c, a) \ b
}

object Mop_VCat extends Mop {
	override def myname = "on"
  override def op(a:FMat, b:FMat, c:Mat):FMat = a on b
  override def op(a:DMat, b:DMat, c:Mat):DMat = a on b
  override def op(a:IMat, b:IMat, c:Mat):IMat = a on b
  override def op(a:LMat, b:LMat, c:Mat):LMat = a on b
  override def op(a:SMat, b:SMat, c:Mat):SMat = a on b
  override def op(a:SDMat, b:SDMat, c:Mat):SDMat = a on b
  override def op(a:CMat, b:CMat, c:Mat):CMat = a on b
  override def op(a:GMat, b:GMat, c:Mat):GMat = getGPair(c, a) on b
  override def op(a:GDMat, b:GDMat, c:Mat):GDMat = getGDPair(c, a) on b
  override def op(a:GIMat, b:GIMat, c:Mat):GIMat = getGIPair(c, a) on b
  override def op(a:GLMat, b:GLMat, c:Mat):GLMat = getGLPair(c, a) on b
}

object Mop_LT extends Mop { 
	override def myname = "<"
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) < b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) < b
  override def op(a:IMat, b:IMat, c:Mat):IMat = getIPair(c, a) < b
  override def op(a:LMat, b:LMat, c:Mat):LMat = getLPair(c, a) < b
  override def op(a:GMat, b:GMat, c:Mat):GMat = getGPair(c, a) < b
  override def op(a:GDMat, b:GDMat, c:Mat):GDMat = getGDPair(c, a) < b
  override def op(a:GIMat, b:GIMat, c:Mat):GIMat = getGIPair(c, a) < b
  override def op(a:GLMat, b:GLMat, c:Mat):GLMat = getGLPair(c, a) < b
  override def sop(a:SMat, b:FMat, c:Mat):SMat = getSPair(c, a) < b
  override def sop(a:GSMat, b:GMat, c:Mat):GSMat = getGSPair(c, a) < b
  override def sop(a:GSDMat, b:GDMat, c:Mat):GSDMat = getGSDPair(c, a) < b
}

object Mop_GT extends Mop { 
	override def myname = ">"
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) > b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) > b
  override def op(a:IMat, b:IMat, c:Mat):IMat = getIPair(c, a) > b
  override def op(a:LMat, b:LMat, c:Mat):LMat = getLPair(c, a) > b
  override def op(a:GMat, b:GMat, c:Mat):GMat = getGPair(c, a) > b
  override def op(a:GDMat, b:GDMat, c:Mat):GDMat = getGDPair(c, a) > b
  override def op(a:GIMat, b:GIMat, c:Mat):GIMat = getGIPair(c, a) > b
  override def op(a:GLMat, b:GLMat, c:Mat):GLMat = getGLPair(c, a) > b
  override def sop(a:SMat, b:FMat, c:Mat):SMat = getSPair(c, a) > b
  override def sop(a:GSMat, b:GMat, c:Mat):GSMat = getGSPair(c, a) > b
  override def sop(a:GSDMat, b:GDMat, c:Mat):GSDMat = getGSDPair(c, a) > b
  
}

object Mop_LE extends Mop {
	override def myname = "<="
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) <= b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) <= b
  override def op(a:IMat, b:IMat, c:Mat):IMat = getIPair(c, a) <= b
  override def op(a:LMat, b:LMat, c:Mat):LMat = getLPair(c, a) <= b
  override def op(a:GMat, b:GMat, c:Mat):GMat = getGPair(c, a) <= b
  override def op(a:GDMat, b:GDMat, c:Mat):GDMat = getGDPair(c, a) <= b
  override def op(a:GIMat, b:GIMat, c:Mat):GIMat = getGIPair(c, a) <= b
  override def op(a:GLMat, b:GLMat, c:Mat):GLMat = getGLPair(c, a) <= b
  override def sop(a:SMat, b:FMat, c:Mat):SMat = getSPair(c, a) <= b
  override def sop(a:GSMat, b:GMat, c:Mat):GSMat = getGSPair(c, a) <= b
  override def sop(a:GSDMat, b:GDMat, c:Mat):GSDMat = getGSDPair(c, a) <= b
}

object Mop_GE extends Mop { 
	override def myname = ">="
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) >= b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) >= b
  override def op(a:IMat, b:IMat, c:Mat):IMat = getIPair(c, a) >= b
  override def op(a:LMat, b:LMat, c:Mat):LMat = getLPair(c, a) >= b
  override def op(a:GMat, b:GMat, c:Mat):GMat = getGPair(c, a) >= b
  override def op(a:GDMat, b:GDMat, c:Mat):GDMat = getGDPair(c, a) >= b
  override def op(a:GIMat, b:GIMat, c:Mat):GIMat = getGIPair(c, a) >= b
  override def op(a:GLMat, b:GLMat, c:Mat):GLMat = getGLPair(c, a) >= b
  override def sop(a:SMat, b:FMat, c:Mat):SMat = getSPair(c, a) >= b
  override def sop(a:GSMat, b:GMat, c:Mat):GSMat = getGSPair(c, a) >= b
  override def sop(a:GSDMat, b:GDMat, c:Mat):GSDMat = getGSDPair(c, a) >= b
}

object Mop_EQ extends Mop { 
	override def myname = "=="
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) == b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) == b
  override def op(a:IMat, b:IMat, c:Mat):IMat = getIPair(c, a) == b
  override def op(a:LMat, b:LMat, c:Mat):LMat = getLPair(c, a) == b
  override def op(a:CMat, b:CMat, c:Mat):CMat = getCPair(c, a) == b
  override def op(a:GMat, b:GMat, c:Mat):GMat = getGPair(c, a) == b
  override def op(a:GDMat, b:GDMat, c:Mat):GDMat = getGDPair(c, a) == b
  override def op(a:GIMat, b:GIMat, c:Mat):GIMat = getGIPair(c, a) == b
  override def op(a:GLMat, b:GLMat, c:Mat):GLMat = getGLPair(c, a) == b
  override def sop(a:SMat, b:FMat, c:Mat):SMat = getSPair(c, a) == b
  override def sop(a:GSMat, b:GMat, c:Mat):GSMat = getGSPair(c, a) == b
  override def sop(a:GSDMat, b:GDMat, c:Mat):GSDMat = getGSDPair(c, a) == b
}

object Mop_NE extends Mop { 
	override def myname = "!="
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) != b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) != b
  override def op(a:IMat, b:IMat, c:Mat):IMat = getIPair(c, a) != b
  override def op(a:LMat, b:LMat, c:Mat):LMat = getLPair(c, a) != b
  override def op(a:CMat, b:CMat, c:Mat):CMat = getCPair(c, a) != b
  override def op(a:GMat, b:GMat, c:Mat):GMat = getGPair(c, a) != b
  override def op(a:GDMat, b:GDMat, c:Mat):GDMat = getGDPair(c, a) != b
  override def op(a:GIMat, b:GIMat, c:Mat):GIMat = getGIPair(c, a) != b
  override def op(a:GLMat, b:GLMat, c:Mat):GLMat = getGLPair(c, a) != b
  override def sop(a:SMat, b:FMat, c:Mat):SMat = getSPair(c, a) != b
  override def sop(a:GSMat, b:GMat, c:Mat):GSMat = getGSPair(c, a) != b
  override def sop(a:GSDMat, b:GDMat, c:Mat):GSDMat = getGSDPair(c, a) != b
}
