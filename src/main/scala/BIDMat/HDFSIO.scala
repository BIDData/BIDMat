package BIDMat


trait HDFSIOtrait {
  def writeMat(name:String, Mat:Mat, comp:Int)
  def writeND(name:String, nd:ND, comp:Int)
  def writeMats(name:String, Mats:Array[Mat], comp:Int)
  def writeNDs(name:String, nds:Array[ND], comp:Int)
  def readMat(name:String, omat:Mat):Mat
  def readND(name:String, omat:ND):ND
  def readMats(name:String, omats:Array[Mat]):Array[Mat]
  def readNDs(name:String, omat:Array[ND]):Array[ND]
}

trait MatIOtrait {
  def get:Array[Mat]
}