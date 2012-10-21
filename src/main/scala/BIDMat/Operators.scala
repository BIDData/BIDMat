package BIDMat

object Operator {
  def applyMat(a:FMat, b:Mat, c:Mop):Mat = {
    b match {
      case fb:FMat => c.fop(a, fb)
      case db:DMat => c.dop(DMat(a), db)
      case ib:IMat => c.fop(a, FMat(ib))
      case cb:CMat => c.cop(CMat(a), cb)
    }
  }
  
  def applyMat(a:DMat, b:Mat, c:Mop):Mat = {
    b match {
      case fb:FMat => c.dop(a, DMat(fb))
      case db:DMat => c.dop(a, db)
      case ib:IMat => c.dop(a, DMat(ib))
      case cb:CMat => c.cop(CMat(a), cb)
    }
  }
  
  def applyMat(a:IMat, b:Mat, c:Mop):Mat = {
    b match {
      case fb:FMat => c.fop(FMat(a), fb)
      case db:DMat => c.dop(DMat(a), db)
      case ib:IMat => c.iop(a, ib)
      case cb:CMat => c.cop(CMat(a), cb)
    }
  }
  
  def applyMat(a:CMat, b:Mat, c:Mop):Mat = {
    b match {
      case fb:FMat => c.cop(a, CMat(fb))
      case db:DMat => c.cop(a, CMat(db))
      case ib:IMat => c.cop(a, CMat(ib))
      case cb:CMat => c.cop(CMat(a), cb)
    }
  }
}

trait Mop {
  def fop(a:FMat, b:FMat):FMat
  def dop(a:DMat, b:DMat):DMat
  def iop(a:IMat, b:IMat):IMat 
  def cop(a:CMat, b:CMat):CMat
  def notImplemented0(s:String, m:Mat):Mat = { 
    throw new RuntimeException("operator "+s+" not implemented for "+m)
  }
}

object Mop_Plus extends Mop { 
  override def fop(a:FMat, b:FMat):FMat = a + b
  override def dop(a:DMat, b:DMat):DMat = a + b
  override def iop(a:IMat, b:IMat):IMat = a + b
  override def cop(a:CMat, b:CMat):CMat = a + b
}

object Mop_Minus extends Mop { 
  override def fop(a:FMat, b:FMat):FMat = a - b
  override def dop(a:DMat, b:DMat):DMat = a - b
  override def iop(a:IMat, b:IMat):IMat = a - b
  override def cop(a:CMat, b:CMat):CMat = a - b
}

object Mop_Times extends Mop { 
  override def fop(a:FMat, b:FMat):FMat = a * b
  override def dop(a:DMat, b:DMat):DMat = a * b
  override def iop(a:IMat, b:IMat):IMat = a * b
  override def cop(a:CMat, b:CMat):CMat = a * b
}

object Mop_Div extends Mop { 
  override def fop(a:FMat, b:FMat):FMat = a / b
  override def dop(a:DMat, b:DMat):DMat = a / b
  override def cop(a:CMat, b:CMat):CMat = a / b
  override def iop(a:IMat, b:IMat):IMat = {notImplemented0("/", a); a}
}

object Mop_RSolve extends Mop { 
  override def fop(a:FMat, b:FMat):FMat = a \\ b
  override def dop(a:DMat, b:DMat):DMat = a \\ b
  override def cop(a:CMat, b:CMat):CMat = a \\ b
  override def iop(a:IMat, b:IMat):IMat = {notImplemented0("/", a); a}
}

object Mop_ETimes extends Mop { 
  override def fop(a:FMat, b:FMat):FMat = a *@ b
  override def dop(a:DMat, b:DMat):DMat = a *@ b
  override def iop(a:IMat, b:IMat):IMat = a *@ b
  override def cop(a:CMat, b:CMat):CMat = a *@ b
}

object Mop_EDiv extends Mop { 
  override def fop(a:FMat, b:FMat):FMat = a /@ b
  override def dop(a:DMat, b:DMat):DMat = a /@ b
  override def iop(a:IMat, b:IMat):IMat = a /@ b
  override def cop(a:CMat, b:CMat):CMat = a /@ b
}

object Mop_HCat extends Mop { 
  override def fop(a:FMat, b:FMat):FMat = a \ b
  override def dop(a:DMat, b:DMat):DMat = a \ b
  override def iop(a:IMat, b:IMat):IMat = a \ b
  override def cop(a:CMat, b:CMat):CMat = {notImplemented0("\\", a); a}
}

object Mop_VCat extends Mop { 
  override def fop(a:FMat, b:FMat):FMat = a on b
  override def dop(a:DMat, b:DMat):DMat = a on b
  override def iop(a:IMat, b:IMat):IMat = a on b
  override def cop(a:CMat, b:CMat):CMat = {notImplemented0("on", a); a}
}

object Mop_LT extends Mop { 
  override def fop(a:FMat, b:FMat):FMat = a < b
  override def dop(a:DMat, b:DMat):DMat = a < b
  override def iop(a:IMat, b:IMat):IMat = a < b
  override def cop(a:CMat, b:CMat):CMat = {notImplemented0("<", a); a}
}

object Mop_GT extends Mop { 
  override def fop(a:FMat, b:FMat):FMat = a > b
  override def dop(a:DMat, b:DMat):DMat = a > b
  override def iop(a:IMat, b:IMat):IMat = a > b
  override def cop(a:CMat, b:CMat):CMat = {notImplemented0(">", a); a}
}

object Mop_LE extends Mop { 
  override def fop(a:FMat, b:FMat):FMat = a <= b
  override def dop(a:DMat, b:DMat):DMat = a <= b
  override def iop(a:IMat, b:IMat):IMat = a <= b
  override def cop(a:CMat, b:CMat):CMat = {notImplemented0("<=", a); a}
}

object Mop_GE extends Mop { 
  override def fop(a:FMat, b:FMat):FMat = a >= b
  override def dop(a:DMat, b:DMat):DMat = a >= b
  override def iop(a:IMat, b:IMat):IMat = a >= b
  override def cop(a:CMat, b:CMat):CMat = {notImplemented0(">=", a); a}
}

object Mop_EQ extends Mop { 
  override def fop(a:FMat, b:FMat):FMat = a == b
  override def dop(a:DMat, b:DMat):DMat = a == b
  override def iop(a:IMat, b:IMat):IMat = a == b
  override def cop(a:CMat, b:CMat):CMat = {notImplemented0("==", a); a}
}

object Mop_NE extends Mop { 
  override def fop(a:FMat, b:FMat):FMat = a != b
  override def dop(a:DMat, b:DMat):DMat = a != b
  override def iop(a:IMat, b:IMat):IMat = a != b
  override def cop(a:CMat, b:CMat):CMat = {notImplemented0("!=", a); a}
}
