package BIDMat
import MatFunctions._

/*
 * Type coercions for operators
 */ 

/* op is the default binary operator and for two same-typed matrices returns a matrix of matching type.
 * sop is a sparse operation that is only relevant between matrices of different types to resolve the ambuiguity of the result.
 * A sparse operation between sparse A and dense B could return a sparse or dense result depending on the operation.
 * We use the "op" operator to return a dense result, while "sop" returns a sparse result. 
 */

trait Mop {
  def myname:String;
  /*
   * Operations on matching types. These need to be overridden by the actual operator.
   */
  def op(a:BMat, b:BMat, c:Mat):BMat = {notImplemented(myname, a, b); a}
  def op(a:FMat, b:FMat, c:Mat):FMat = {notImplemented(myname, a, b); a}
  def op(a:DMat, b:DMat, c:Mat):DMat = {notImplemented(myname, a, b); a}
  def op(a:IMat, b:IMat, c:Mat):IMat = {notImplemented(myname, a, b); a}
  def op(a:LMat, b:LMat, c:Mat):LMat = {notImplemented(myname, a, b); a}
  def op(a:CMat, b:CMat, c:Mat):CMat = {notImplemented(myname, a, b); a}

  def op(a:SMat, b:SMat, c:Mat):SMat = {notImplemented(myname, a, b); a}
  def op(a:SDMat, b:SDMat, c:Mat):SDMat = {notImplemented(myname, a, b); a}
  
  def sop(a:SMat, b:SMat, c:Mat):SMat = {notImplemented(myname, a, b); a}
  def sop(a:SDMat, b:SDMat, c:Mat):SDMat = {notImplemented(myname, a, b); a}

  def sop(a:SMat, b:FMat, c:Mat):SMat = {notImplemented(myname, a, b); a}
  def sop(a:SDMat, b:DMat, c:Mat):SDMat = {notImplemented(myname, a, b); a}

  def top(a:TMat, b:FMat, c:Mat):Mat = {notImplemented(myname, a, b); a}
  def top(a:TMat, b:SMat, c:Mat):Mat = {notImplemented(myname, a, b); a}
  def top(a:TMat, b:GMat, c:Mat):Mat = {notImplemented(myname, a, b); a}
  def top(a:TMat, b:GSMat, c:Mat):Mat = {notImplemented(myname, a, b); a}
  
  /*
   * Default coercions when first mat is an FMat. These can be overridden. 
   */
  def op(a:FMat, b:SMat, c:Mat):FMat = op(a, full(b), c)
  def op(a:FMat, b:DMat, c:Mat):DMat = op(DMat(a), b, c)
  def op(a:FMat, b:IMat, c:Mat):FMat = op(a, FMat(b), c)  
  def op(a:FMat, b:LMat, c:Mat):FMat = op(a, FMat(b), c)
  def op(a:FMat, b:CMat, c:Mat):CMat = op(CMat(a), b, c)  

  def op(a:FMat, b:Mat, c:Mat):Mat = {
    b match {
      case bb:FMat => op(a, bb, c)
      case bb:SMat => op(a, bb, c)
      case bb:DMat => op(a, bb, c)
      case bb:IMat => op(a, bb, c)
      case bb:LMat => op(a, bb, c)
      case bb:CMat => op(a, bb, c)
    }
  } 
  
  /*
   * Default coercions when first arg is a DMat. 
   */
  def op(a:DMat, b:FMat, c:Mat):DMat = op(a, DMat(b), c)
  def op(a:DMat, b:IMat, c:Mat):DMat = op(a, DMat(b), c)
  def op(a:DMat, b:LMat, c:Mat):DMat = op(a, DMat(b), c)
  def op(a:DMat, b:SDMat, c:Mat):DMat = op(a, full(b), c)
  def op(a:DMat, b:CMat, c:Mat):CMat = op(CMat(a), b, c)
 
  def op(a:DMat, b:Mat, c:Mat):Mat = {
    b match {
      case bb:FMat => op(a, bb, c)
      case bb:DMat => op(a, bb, c)
      case bb:IMat => op(a, bb, c)
      case bb:LMat => op(a, bb, c)
      case bb:SDMat => op(a, bb, c)
      case bb:CMat => op(a, bb, c)
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
  
  def op(a:IMat, b:Mat, c:Mat):Mat = {
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
      case bb:SDMat => op(a, bb, c)
    }
  }
  
  def sop(a:SDMat, b:Mat, c:Mat):Mat = {
    b match {
      case bb:DMat => op(a, bb, c)
      case bb:SDMat => op(a, bb, c)
    }
  }

  /*
   * TMat coercion
   */ 

  def op(a:TMat, b:Mat, c:Mat):Mat = {
    b match {
      case bb:GMat => top(a, bb, c);
      case bb:GSMat => top(a, bb, c);
      case bb:FMat => top(a, bb, c);
      case bb:SMat => top(a, bb, c);
    }
  }

  
  /*
   * No clue
   */
  def op(a:Mat, b:Mat, c:Mat):Mat = {
    a match {
      case aa:BMat => op(aa, b, c)
      case aa:FMat => op(aa, b, c)
      case aa:DMat => op(aa, b, c)
      case aa:IMat => op(aa, b, c)
      case aa:CMat => op(aa, b, c)
      case aa:SMat => op(aa, b, c)
      case aa:SDMat => op(aa, b, c)
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
  def getTPair(c:Mat, a:TMat):TPair = new TPair(c, a)
}

object Mop_Plus extends Mop { // Plus includes only operations on matching types
  override def myname = "+"
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) + b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) + b
  override def op(a:IMat, b:IMat, c:Mat):IMat = getIPair(c, a) + b
  override def op(a:LMat, b:LMat, c:Mat):LMat = getLPair(c, a) + b
  override def op(a:CMat, b:CMat, c:Mat):CMat = getCPair(c, a) + b
  override def op(a:SMat, b:SMat, c:Mat):SMat = getSPair(c, a) + b
  override def sop(a:SMat, b:FMat, c:Mat):SMat = getSPair(c, a) + b

}

object Mop_Minus extends Mop { 
  override def myname = "-"
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) - b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) - b
  override def op(a:IMat, b:IMat, c:Mat):IMat = getIPair(c, a) - b
  override def op(a:LMat, b:LMat, c:Mat):LMat = getLPair(c, a) - b
  override def op(a:CMat, b:CMat, c:Mat):CMat = getCPair(c, a) - b
  override def op(a:SMat, b:SMat, c:Mat):SMat = getSPair(c, a) - b
  override def sop(a:SMat, b:FMat, c:Mat):SMat = getSPair(c, a) - b
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
  override def op(a:TMat, b:Mat, c:Mat):Mat = getTPair(c, a) * b
}

object Mop_TimesT extends Mop { 
  override def myname = "*^"
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) *^ b
  override def op(a:FMat, b:SMat, c:Mat):FMat = getFPair(c, a) *^ b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) *^ b
  override def op(a:DMat, b:SDMat, c:Mat):DMat = getDPair(c, a) *^ b
}

object Mop_TTimes extends Mop { 
  override def myname = "^*"
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) Tx b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) Tx b
  override def op(a:SMat, b:FMat, c:Mat):FMat = getSPair(c, a) Tx b
  override def op(a:SDMat, b:DMat, c:Mat):DMat = getSDPair(c, a) Tx b
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
  override def op(a:SMat, b:SMat, c:Mat):SMat = getSPair(c, a) *@ b
  override def sop(a:SMat, b:FMat, c:Mat):SMat = getSPair(c, a) *@ b
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
  override def op(a:SMat, b:SMat, c:Mat):SMat = getSPair(c, a) / b
  override def sop(a:SMat, b:FMat, c:Mat):SMat = getSPair(c, a) / b
}

object Mop_Dot extends Mop { 
  override def myname = "dot"
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) dot b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) dot b
  override def op(a:IMat, b:IMat, c:Mat):IMat = getIPair(c, a) dot b
  override def op(a:LMat, b:LMat, c:Mat):LMat = getLPair(c, a) dot b
  override def op(a:CMat, b:CMat, c:Mat):CMat = getCPair(c, a) dot b
}

object Mop_Dotr extends Mop { 
  override def myname = "dotr"
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) dotr b
  override def op(a:IMat, b:IMat, c:Mat):IMat = getIPair(c, a) dotr b
  override def op(a:LMat, b:LMat, c:Mat):LMat = getLPair(c, a) dotr b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) dotr b
}

object Mop_Pow extends Mop {
  override def myname = "^"
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) ^ b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) ^ b
//  override def op(a:TMat, b:TMat, c:Mat):TMat = getTPair(c, a) ^ b
}

object Mop_HCat extends Mop { 
	override def myname = "\\"
  override def op(a:BMat, b:BMat, c:Mat):BMat = a \ b
  override def op(a:FMat, b:FMat, c:Mat):FMat = a \ b
  override def op(a:DMat, b:DMat, c:Mat):DMat = a \ b
  override def op(a:IMat, b:IMat, c:Mat):IMat = a \ b
  override def op(a:LMat, b:LMat, c:Mat):LMat = a \ b
  override def op(a:SMat, b:SMat, c:Mat):SMat = a \ b
  override def op(a:SDMat, b:SDMat, c:Mat):SDMat = a \ b
  override def op(a:CMat, b:CMat, c:Mat):CMat = a \ b
}

object Mop_VCat extends Mop {
  override def myname = "on"
  override def op(a:BMat, b:BMat, c:Mat):BMat = a on b
  override def op(a:FMat, b:FMat, c:Mat):FMat = a on b
  override def op(a:DMat, b:DMat, c:Mat):DMat = a on b
  override def op(a:IMat, b:IMat, c:Mat):IMat = a on b
  override def op(a:LMat, b:LMat, c:Mat):LMat = a on b
  override def op(a:SMat, b:SMat, c:Mat):SMat = a on b
  override def op(a:SDMat, b:SDMat, c:Mat):SDMat = a on b
  override def op(a:CMat, b:CMat, c:Mat):CMat = a on b
}

object Mop_LT extends Mop { 
	override def myname = "<"
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) < b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) < b
  override def op(a:IMat, b:IMat, c:Mat):IMat = getIPair(c, a) < b
  override def op(a:LMat, b:LMat, c:Mat):LMat = getLPair(c, a) < b
  override def sop(a:SMat, b:FMat, c:Mat):SMat = getSPair(c, a) < b
  override def sop(a:SDMat, b:DMat, c:Mat):SDMat = getSDPair(c, a) < b
}

object Mop_GT extends Mop { 
	override def myname = ">"
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) > b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) > b
  override def op(a:IMat, b:IMat, c:Mat):IMat = getIPair(c, a) > b
  override def op(a:LMat, b:LMat, c:Mat):LMat = getLPair(c, a) > b
  override def sop(a:SMat, b:FMat, c:Mat):SMat = getSPair(c, a) > b
  override def sop(a:SDMat, b:DMat, c:Mat):SDMat = getSDPair(c, a) > b
  
}

object Mop_LE extends Mop {
	override def myname = "<="
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) <= b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) <= b
  override def op(a:IMat, b:IMat, c:Mat):IMat = getIPair(c, a) <= b
  override def op(a:LMat, b:LMat, c:Mat):LMat = getLPair(c, a) <= b
  override def sop(a:SMat, b:FMat, c:Mat):SMat = getSPair(c, a) <= b
  override def sop(a:SDMat, b:DMat, c:Mat):SDMat = getSDPair(c, a) <= b
}

object Mop_GE extends Mop { 
	override def myname = ">="
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) >= b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) >= b
  override def op(a:IMat, b:IMat, c:Mat):IMat = getIPair(c, a) >= b
  override def op(a:LMat, b:LMat, c:Mat):LMat = getLPair(c, a) >= b
  override def sop(a:SMat, b:FMat, c:Mat):SMat = getSPair(c, a) >= b
  override def sop(a:SDMat, b:DMat, c:Mat):SDMat = getSDPair(c, a) >= b
}

object Mop_EQ extends Mop { 
	override def myname = "=="
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) == b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) == b
  override def op(a:IMat, b:IMat, c:Mat):IMat = getIPair(c, a) == b
  override def op(a:LMat, b:LMat, c:Mat):LMat = getLPair(c, a) == b
  override def op(a:CMat, b:CMat, c:Mat):CMat = getCPair(c, a) == b
  override def sop(a:SMat, b:FMat, c:Mat):SMat = getSPair(c, a) == b
  override def sop(a:SDMat, b:DMat, c:Mat):SDMat = getSDPair(c, a) == b
}

object Mop_NE extends Mop { 
	override def myname = "!="
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) != b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) != b
  override def op(a:IMat, b:IMat, c:Mat):IMat = getIPair(c, a) != b
  override def op(a:LMat, b:LMat, c:Mat):LMat = getLPair(c, a) != b
  override def op(a:CMat, b:CMat, c:Mat):CMat = getCPair(c, a) != b
  override def sop(a:SMat, b:FMat, c:Mat):SMat = getSPair(c, a) != b
  override def sop(a:SDMat, b:DMat, c:Mat):SDMat = getSDPair(c, a) != b
}

object Mop_Max extends Mop { 
  override def myname = "max"
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) max b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) max b
  override def op(a:IMat, b:IMat, c:Mat):IMat = getIPair(c, a) max b
  override def op(a:LMat, b:LMat, c:Mat):LMat = getLPair(c, a) max b
  override def sop(a:SMat, b:FMat, c:Mat):SMat = getSPair(c, a) max b
  override def sop(a:SDMat, b:DMat, c:Mat):SDMat = getSDPair(c, a) max b
}

object Mop_Min extends Mop { 
  override def myname = "min"
  override def op(a:FMat, b:FMat, c:Mat):FMat = getFPair(c, a) min b
  override def op(a:DMat, b:DMat, c:Mat):DMat = getDPair(c, a) min b
  override def op(a:IMat, b:IMat, c:Mat):IMat = getIPair(c, a) min b
  override def op(a:LMat, b:LMat, c:Mat):LMat = getLPair(c, a) min b
  override def sop(a:SMat, b:FMat, c:Mat):SMat = getSPair(c, a) min b
  override def sop(a:SDMat, b:DMat, c:Mat):SDMat = getSDPair(c, a) min b
}
