package BIDMat
import MatFunctions._

object Operator {
  def applyMat(a:FMat, b:Mat, c:Mat, op:Mop):Mat = {
    b match {
      case fb:FMat => op.fop(a, fb, c)
      case sb:SMat => op.fop(a, sb, c)
      case db:DMat => op.dop(DMat(a), db, c)
      case ib:IMat => op.fop(a, FMat(ib), c)
      case cb:CMat => op.cop(CMat(a), cb, c)
      case gb:GMat => op.gop(GMat(a), gb, c)
    }
  }
  
  def applyMat(a:DMat, b:Mat, c:Mat, op:Mop):Mat = {
    b match {
      case fb:FMat => op.dop(a, DMat(fb), c)
      case db:DMat => op.dop(a, db, c)
      case ib:IMat => op.dop(a, DMat(ib), c)
      case cb:CMat => op.cop(CMat(a), cb, c)
    }
  }
  
  def applyMat(a:IMat, b:Mat, c:Mat, op:Mop):Mat = {
    b match {
      case fb:FMat => op.fop(FMat(a), fb, c)
      case db:DMat => op.dop(DMat(a), db, c)
      case ib:IMat => op.iop(a, ib, c)
      case cb:CMat => op.cop(CMat(a), cb, c)
      case gb:GMat => op.gop(GMat(FMat(a)), gb, c)
    }
  }
  
  def applyMat(a:CMat, b:Mat, c:Mat, op:Mop):Mat = {
    b match {
      case fb:FMat => op.cop(a, CMat(fb), c)
      case db:DMat => op.cop(a, CMat(db), c)
      case ib:IMat => op.cop(a, CMat(ib), c)
      case cb:CMat => op.cop(CMat(a), cb, c)
    }
  }
  
  def applyMat(a:GMat, b:Mat, c:Mat, op:Mop):Mat = {
    b match {
      case gb:GMat => op.gop(a, gb, c)
    }
  }
  
  def applyMat(a:SMat, b:Mat, c:Mat, op:Mop):Mat = {
    b match {
      case sb:SMat => op.sop(a, sb, c)
    }
  }
  
  def multDim1(a:Mat, b:Mat):Int = {
    if (a.nrows == 1 && a.ncols == 1) {
      b.nrows
    } else {
      a.nrows
    }
  } 
  
  def multDim2(a:Mat, b:Mat):Int = {
    if (b.nrows == 1 && b.ncols == 1) {
      a.ncols
    } else {
      b.ncols
    }
  } 
  
  def getFPair(c:Mat, a:FMat):FPair = {
    if (c.asInstanceOf[AnyRef] != null) {
      new FPair(c, a)
    } else {
      new FPair(FMat(a.nrows, a.ncols), a)
    }
  }
  
  def getFPair(c:Mat, a:FMat, b:FMat):FPair = {
    if (c.asInstanceOf[AnyRef] != null) {
      new FPair(c, a)
    } else {
      new FPair(FMat(multDim1(a,b), multDim2(a,b)), a)
    }
  }
  
  def getDPair(c:Mat, a:DMat):DPair = {
    if (c.asInstanceOf[AnyRef] != null) {
      new DPair(c, a)
    } else {
      new DPair(DMat(a.nrows, a.ncols), a)
    }
  }
  
  def getDPair(c:Mat, a:DMat, b:DMat):DPair = {
    if (c.asInstanceOf[AnyRef] != null) {
      new DPair(c, a)
    } else {
      new DPair(DMat(multDim1(a,b), multDim2(a,b)), a)
    }
  }
  
  def getIPair(c:Mat, a:IMat):IPair = {
    if (c.asInstanceOf[AnyRef] != null) {
      new IPair(c, a)
    } else {
      new IPair(IMat(a.nrows, a.ncols), a)
    }
  }
  
  def getIPair(c:Mat, a:IMat, b:IMat):IPair = {
    if (c.asInstanceOf[AnyRef] != null) {
      new IPair(c, a)
    } else {
      new IPair(IMat(multDim1(a,b), multDim2(a,b)), a)
    }
  }
  
  def getCPair(c:Mat, a:CMat):CPair = {
    if (c.asInstanceOf[AnyRef] != null) {
      new CPair(c, a)
    } else {
      new CPair(CMat(a.nrows, a.ncols), a)
    }
  }
  
  def getCPair(c:Mat, a:CMat, b:CMat):CPair = {
    if (c.asInstanceOf[AnyRef] != null) {
      new CPair(c, a)
    } else {
      new CPair(CMat(multDim1(a,b), multDim2(a,b)), a)
    }
  }
  
  def getGPair(c:Mat, a:GMat):GPair = {
    if (c.asInstanceOf[AnyRef] != null) {
      new GPair(c, a)
    } else {
      new GPair(GMat(a.nrows, a.ncols), a)
    }
  }
  
  def getGPair(c:Mat, a:GMat, b:GMat):GPair = {
    if (c.asInstanceOf[AnyRef] != null) {
      new GPair(c, a)
    } else {
      new GPair(GMat(multDim1(a,b), multDim2(a,b)), a)
    }
  }
  
   def getSPair(c:Mat, a:SMat):SPair = {
    if (c.asInstanceOf[AnyRef] != null) {
      new SPair(c, a)
    } else {
      new SPair(SMat(a.nrows, a.ncols, a.nnz), a)
    }
  }
}

trait Mop {
  def fop(a:FMat, b:FMat, c:Mat):FMat
  def fop(a:FMat, b:SMat, c:Mat):FMat
  def dop(a:DMat, b:DMat, c:Mat):DMat
  def iop(a:IMat, b:IMat, c:Mat):IMat 
  def cop(a:CMat, b:CMat, c:Mat):CMat
  def gop(a:GMat, b:GMat, c:Mat):GMat
  def sop(a:SMat, b:SMat, c:Mat):SMat
  def notImplemented0(s:String, m:Mat):Mat = { 
    throw new RuntimeException("operator "+s+" not implemented for "+m.mytype)
  }
}

object Mop_Plus extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = Operator.getFPair(c, a) + b
  override def fop(a:FMat, b:SMat, c:Mat):FMat = Operator.getFPair(c, a) + full(b)
  override def dop(a:DMat, b:DMat, c:Mat):DMat = Operator.getDPair(c, a) + b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = Operator.getIPair(c, a) + b
  override def cop(a:CMat, b:CMat, c:Mat):CMat = Operator.getCPair(c, a) + b
  override def gop(a:GMat, b:GMat, c:Mat):GMat = Operator.getGPair(c, a) + b
  override def sop(a:SMat, b:SMat, c:Mat):SMat = Operator.getSPair(c, a) + b
}

object Mop_Minus extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = Operator.getFPair(c, a) - b
  override def fop(a:FMat, b:SMat, c:Mat):FMat = Operator.getFPair(c, a) - full(b)
  override def dop(a:DMat, b:DMat, c:Mat):DMat = Operator.getDPair(c, a) - b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = Operator.getIPair(c, a) - b
  override def cop(a:CMat, b:CMat, c:Mat):CMat = Operator.getCPair(c, a) - b
  override def gop(a:GMat, b:GMat, c:Mat):GMat = Operator.getGPair(c, a) - b
  override def sop(a:SMat, b:SMat, c:Mat):SMat = Operator.getSPair(c, a) - b
}

object Mop_Times extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = Operator.getFPair(c, a, b) * b
  override def fop(a:FMat, b:SMat, c:Mat):FMat = Operator.getFPair(c, a) * b
  override def dop(a:DMat, b:DMat, c:Mat):DMat = Operator.getDPair(c, a, b) * b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = Operator.getIPair(c, a, b) * b
  override def cop(a:CMat, b:CMat, c:Mat):CMat = Operator.getCPair(c, a, b) * b
  override def gop(a:GMat, b:GMat, c:Mat):GMat = Operator.getGPair(c, a, b) * b
  override def sop(a:SMat, b:SMat, c:Mat):SMat = {notImplemented0("*", a); a}
}

object Mop_Div extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = a /< b
  override def fop(a:FMat, b:SMat, c:Mat):FMat = {notImplemented0("/<", a); a}
  override def dop(a:DMat, b:DMat, c:Mat):DMat = a /< b
  override def cop(a:CMat, b:CMat, c:Mat):CMat = a /< b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = {notImplemented0("/<", a); a}
  override def gop(a:GMat, b:GMat, c:Mat):GMat = {notImplemented0("/<", a); a}
  override def sop(a:SMat, b:SMat, c:Mat):SMat = {notImplemented0("/<", a); a}
}

object Mop_RSolve extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = a \\ b
  override def fop(a:FMat, b:SMat, c:Mat):FMat = {notImplemented0("\\\\", a); a}
  override def dop(a:DMat, b:DMat, c:Mat):DMat = a \\ b
  override def cop(a:CMat, b:CMat, c:Mat):CMat = a \\ b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = {notImplemented0("\\\\", a); a}
  override def gop(a:GMat, b:GMat, c:Mat):GMat = {notImplemented0("\\\\", a); a}
  override def sop(a:SMat, b:SMat, c:Mat):SMat = {notImplemented0("\\\\", a); a}
}

object Mop_ETimes extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = Operator.getFPair(c, a) *@ b
  override def fop(a:FMat, b:SMat, c:Mat):FMat = Operator.getFPair(c, a) *@ full(b)
  override def dop(a:DMat, b:DMat, c:Mat):DMat = Operator.getDPair(c, a) *@ b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = Operator.getIPair(c, a) *@ b
  override def cop(a:CMat, b:CMat, c:Mat):CMat = Operator.getCPair(c, a) *@ b
  override def gop(a:GMat, b:GMat, c:Mat):GMat = Operator.getGPair(c, a) *@ b
  override def sop(a:SMat, b:SMat, c:Mat):SMat = Operator.getSPair(c, a) *@ b
}

object Mop_EDiv extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = Operator.getFPair(c, a) / b
  override def fop(a:FMat, b:SMat, c:Mat):FMat = Operator.getFPair(c, a) / full(b)
  override def dop(a:DMat, b:DMat, c:Mat):DMat = Operator.getDPair(c, a) / b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = {notImplemented0("/", a); a}
  override def cop(a:CMat, b:CMat, c:Mat):CMat = Operator.getCPair(c, a) / b
  override def gop(a:GMat, b:GMat, c:Mat):GMat = Operator.getGPair(c, a) / b
  override def sop(a:SMat, b:SMat, c:Mat):SMat = Operator.getSPair(c, a) / b
}

object Mop_HCat extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = a \ b
  override def fop(a:FMat, b:SMat, c:Mat):FMat = {notImplemented0("\\", a); a}
  override def dop(a:DMat, b:DMat, c:Mat):DMat = a \ b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = a \ b
  override def cop(a:CMat, b:CMat, c:Mat):CMat = a \ b
  override def gop(a:GMat, b:GMat, c:Mat):GMat = {notImplemented0("\\", a); a}
  override def sop(a:SMat, b:SMat, c:Mat):SMat = {notImplemented0("\\", a); a}
}

object Mop_VCat extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = a on b
  override def fop(a:FMat, b:SMat, c:Mat):FMat = {notImplemented0("on", a); a}
  override def dop(a:DMat, b:DMat, c:Mat):DMat = a on b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = a on b
  override def cop(a:CMat, b:CMat, c:Mat):CMat = a on b
  override def gop(a:GMat, b:GMat, c:Mat):GMat = {notImplemented0("on", a); a}
  override def sop(a:SMat, b:SMat, c:Mat):SMat = {notImplemented0("on", a); a}
}

object Mop_LT extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = Operator.getFPair(c, a) < b
  override def fop(a:FMat, b:SMat, c:Mat):FMat = Operator.getFPair(c, a) < full(b)
  override def dop(a:DMat, b:DMat, c:Mat):DMat = Operator.getDPair(c, a) < b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = Operator.getIPair(c, a) < b
  override def cop(a:CMat, b:CMat, c:Mat):CMat = {notImplemented0("<", a); a}
  override def gop(a:GMat, b:GMat, c:Mat):GMat = Operator.getGPair(c, a) < b
  override def sop(a:SMat, b:SMat, c:Mat):SMat = {notImplemented0("<", a); a}
}

object Mop_GT extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = Operator.getFPair(c, a) > b
  override def fop(a:FMat, b:SMat, c:Mat):FMat = Operator.getFPair(c, a) > full(b)
  override def dop(a:DMat, b:DMat, c:Mat):DMat = Operator.getDPair(c, a) > b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = Operator.getIPair(c, a) > b
  override def cop(a:CMat, b:CMat, c:Mat):CMat = {notImplemented0(">", a); a}
  override def gop(a:GMat, b:GMat, c:Mat):GMat = Operator.getGPair(c, a) > b
  override def sop(a:SMat, b:SMat, c:Mat):SMat = {notImplemented0(">", a); a}
}

object Mop_LE extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = Operator.getFPair(c, a) <= b
  override def fop(a:FMat, b:SMat, c:Mat):FMat = Operator.getFPair(c, a) <= full(b)
  override def dop(a:DMat, b:DMat, c:Mat):DMat = Operator.getDPair(c, a) <= b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = Operator.getIPair(c, a) <= b
  override def cop(a:CMat, b:CMat, c:Mat):CMat = {notImplemented0("<=", a); a}
  override def gop(a:GMat, b:GMat, c:Mat):GMat = Operator.getGPair(c, a) <= b
  override def sop(a:SMat, b:SMat, c:Mat):SMat = {notImplemented0("<=", a); a}
}

object Mop_GE extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = Operator.getFPair(c, a) >= b
  override def fop(a:FMat, b:SMat, c:Mat):FMat = Operator.getFPair(c, a) >= full(b)
  override def dop(a:DMat, b:DMat, c:Mat):DMat = Operator.getDPair(c, a) >= b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = Operator.getIPair(c, a) >= b
  override def cop(a:CMat, b:CMat, c:Mat):CMat = {notImplemented0(">=", a); a}
  override def gop(a:GMat, b:GMat, c:Mat):GMat = Operator.getGPair(c, a) >= b
  override def sop(a:SMat, b:SMat, c:Mat):SMat = {notImplemented0(">=", a); a}
}

object Mop_EQ extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = Operator.getFPair(c, a) == b
  override def fop(a:FMat, b:SMat, c:Mat):FMat = Operator.getFPair(c, a) == full(b)
  override def dop(a:DMat, b:DMat, c:Mat):DMat = Operator.getDPair(c, a) == b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = Operator.getIPair(c, a) == b
  override def cop(a:CMat, b:CMat, c:Mat):CMat = Operator.getCPair(c, a) == b
  override def gop(a:GMat, b:GMat, c:Mat):GMat = Operator.getGPair(c, a) == b
  override def sop(a:SMat, b:SMat, c:Mat):SMat = {notImplemented0("==", a); a}
}

object Mop_NE extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = Operator.getFPair(c, a) != b
  override def fop(a:FMat, b:SMat, c:Mat):FMat = Operator.getFPair(c, a) != full(b)
  override def dop(a:DMat, b:DMat, c:Mat):DMat = Operator.getDPair(c, a) != b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = Operator.getIPair(c, a) != b
  override def cop(a:CMat, b:CMat, c:Mat):CMat = Operator.getCPair(c, a) != b
  override def gop(a:GMat, b:GMat, c:Mat):GMat = Operator.getGPair(c, a) != b
  override def sop(a:SMat, b:SMat, c:Mat):SMat = {notImplemented0("!=", a); a}
}
