package BIDMat

object Operator {
  def applyMat(a:FMat, b:Mat, c:Mat, op:Mop):Mat = {
    b match {
      case fb:FMat => op.fop(a, fb, c)
      case db:DMat => op.dop(DMat(a), db, c)
      case ib:IMat => op.fop(a, FMat(ib), c)
      case cb:CMat => op.cop(CMat(a), cb, c)
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
}

trait Mop {
  def fop(a:FMat, b:FMat, c:Mat):FMat
  def dop(a:DMat, b:DMat, c:Mat):DMat
  def iop(a:IMat, b:IMat, c:Mat):IMat 
  def cop(a:CMat, b:CMat, c:Mat):CMat
  def gop(a:GMat, b:GMat, c:Mat):GMat
  def notImplemented0(s:String, m:Mat):Mat = { 
    throw new RuntimeException("operator "+s+" not implemented for "+m)
  }
}

object Mop_Plus extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = c.asInstanceOf[FMat] ~ a + b
  override def dop(a:DMat, b:DMat, c:Mat):DMat = c.asInstanceOf[DMat] ~ a + b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = c.asInstanceOf[IMat] ~ a + b
  override def cop(a:CMat, b:CMat, c:Mat):CMat = c.asInstanceOf[CMat] ~ a + b
  override def gop(a:GMat, b:GMat, c:Mat):GMat = c.asInstanceOf[GMat] ~ a + b
}

object Mop_Minus extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = c.asInstanceOf[FMat] ~ a - b
  override def dop(a:DMat, b:DMat, c:Mat):DMat = c.asInstanceOf[DMat] ~ a - b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = c.asInstanceOf[IMat] ~ a - b
  override def cop(a:CMat, b:CMat, c:Mat):CMat = a - b
  override def gop(a:GMat, b:GMat, c:Mat):GMat = c.asInstanceOf[GMat] ~ a - b
}

object Mop_Times extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = c.asInstanceOf[FMat] ~ a * b
  override def dop(a:DMat, b:DMat, c:Mat):DMat = c.asInstanceOf[DMat] ~ a * b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = c.asInstanceOf[IMat] ~ a * b
  override def cop(a:CMat, b:CMat, c:Mat):CMat = a * b
  override def gop(a:GMat, b:GMat, c:Mat):GMat = c.asInstanceOf[GMat] ~ a * b
}

object Mop_Div extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = a / b
  override def dop(a:DMat, b:DMat, c:Mat):DMat = a / b
  override def cop(a:CMat, b:CMat, c:Mat):CMat = a / b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = {notImplemented0("/", a); a}
  override def gop(a:GMat, b:GMat, c:Mat):GMat = {notImplemented0("/", a); a}
}

object Mop_RSolve extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = a \\ b
  override def dop(a:DMat, b:DMat, c:Mat):DMat = a \\ b
  override def cop(a:CMat, b:CMat, c:Mat):CMat = a \\ b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = {notImplemented0("/", a); a}
  override def gop(a:GMat, b:GMat, c:Mat):GMat = {notImplemented0("/", a); a}
}

object Mop_ETimes extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = c.asInstanceOf[FMat] ~ a *@ b
  override def dop(a:DMat, b:DMat, c:Mat):DMat = c.asInstanceOf[DMat] ~ a *@ b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = c.asInstanceOf[IMat] ~ a *@ b
  override def cop(a:CMat, b:CMat, c:Mat):CMat = a *@ b
  override def gop(a:GMat, b:GMat, c:Mat):GMat = c.asInstanceOf[GMat] ~ a *@ b
}

object Mop_EDiv extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = c.asInstanceOf[FMat] ~ a /@ b
  override def dop(a:DMat, b:DMat, c:Mat):DMat = c.asInstanceOf[DMat] ~ a /@ b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = {notImplemented0("/", a); a}
  override def cop(a:CMat, b:CMat, c:Mat):CMat = a /@ b
  override def gop(a:GMat, b:GMat, c:Mat):GMat = c.asInstanceOf[GMat] ~ a /@ b
}

object Mop_HCat extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = a \ b
  override def dop(a:DMat, b:DMat, c:Mat):DMat = a \ b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = a \ b
  override def cop(a:CMat, b:CMat, c:Mat):CMat = a \ b
  override def gop(a:GMat, b:GMat, c:Mat):GMat = {notImplemented0("/", a); a}
}

object Mop_VCat extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = a on b
  override def dop(a:DMat, b:DMat, c:Mat):DMat = a on b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = a on b
  override def cop(a:CMat, b:CMat, c:Mat):CMat = a on b
  override def gop(a:GMat, b:GMat, c:Mat):GMat = {notImplemented0("/", a); a}
}

object Mop_LT extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = c.asInstanceOf[FMat] ~ a < b
  override def dop(a:DMat, b:DMat, c:Mat):DMat = c.asInstanceOf[DMat] ~ a < b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = c.asInstanceOf[IMat] ~ a < b
  override def cop(a:CMat, b:CMat, c:Mat):CMat = {notImplemented0("<", a); a}
  override def gop(a:GMat, b:GMat, c:Mat):GMat = c.asInstanceOf[GMat] ~ a < b
}

object Mop_GT extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = c.asInstanceOf[FMat] ~ a > b
  override def dop(a:DMat, b:DMat, c:Mat):DMat = c.asInstanceOf[DMat] ~ a > b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = c.asInstanceOf[IMat] ~ a > b
  override def cop(a:CMat, b:CMat, c:Mat):CMat = {notImplemented0(">", a); a}
  override def gop(a:GMat, b:GMat, c:Mat):GMat = c.asInstanceOf[GMat] ~ a > b
}

object Mop_LE extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = c.asInstanceOf[FMat] ~ a <= b
  override def dop(a:DMat, b:DMat, c:Mat):DMat = c.asInstanceOf[DMat] ~ a <= b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = c.asInstanceOf[IMat] ~ a <= b
  override def cop(a:CMat, b:CMat, c:Mat):CMat = {notImplemented0("<=", a); a}
  override def gop(a:GMat, b:GMat, c:Mat):GMat = c.asInstanceOf[GMat] ~ a <= b
}

object Mop_GE extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = c.asInstanceOf[FMat] ~ a >= b
  override def dop(a:DMat, b:DMat, c:Mat):DMat = c.asInstanceOf[DMat] ~ a >= b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = c.asInstanceOf[IMat] ~ a >= b
  override def cop(a:CMat, b:CMat, c:Mat):CMat = {notImplemented0(">=", a); a}
  override def gop(a:GMat, b:GMat, c:Mat):GMat = c.asInstanceOf[GMat] ~ a >= b
}

object Mop_EQ extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = c.asInstanceOf[FMat] ~ a == b
  override def dop(a:DMat, b:DMat, c:Mat):DMat = c.asInstanceOf[DMat] ~ a == b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = c.asInstanceOf[IMat] ~ a == b
  override def cop(a:CMat, b:CMat, c:Mat):CMat = {notImplemented0("==", a); a}
  override def gop(a:GMat, b:GMat, c:Mat):GMat = c.asInstanceOf[GMat] ~ a == b
}

object Mop_NE extends Mop { 
  override def fop(a:FMat, b:FMat, c:Mat):FMat = c.asInstanceOf[FMat] ~ a != b
  override def dop(a:DMat, b:DMat, c:Mat):DMat = c.asInstanceOf[DMat] ~ a != b
  override def iop(a:IMat, b:IMat, c:Mat):IMat = c.asInstanceOf[IMat] ~ a != b
  override def cop(a:CMat, b:CMat, c:Mat):CMat = {notImplemented0("!=", a); a}
  override def gop(a:GMat, b:GMat, c:Mat):GMat = c.asInstanceOf[GMat] ~ a != b
}
