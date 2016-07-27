package BIDMat


trait HDFSIOtrait extends Serializable {
  def writeMat(name:String, Mat:Mat, comp:Int)
  def writeND(name:String, nd:ND, comp:Int)
  def writeMats(name:String, Mats:Array[ND], comp:Int)
  def writeNDs(name:String, nds:Array[ND], comp:Int)
  def readMat(name:String, omat:Mat):Mat
  def readND(name:String, omat:ND):ND
  def readMats(name:String, omats:Array[ND]):Array[ND]
  def readNDs(name:String, omat:Array[ND]):Array[ND]
}

trait MatIOtrait extends Serializable {
  def get:Array[ND]
}