package BIDMat


trait HDFSIOtrait extends Serializable {
  def writeMat(name:String, Mat:Mat, comp:Int)
  def writeMats(name:String, Mats:Array[Mat], comp:Int)
  def readMat(name:String, omat:Mat):Mat
  def readMats(name:String, omats:Array[Mat]):Array[Mat]
}

trait MatIOtrait extends Serializable {
  def get:Array[Mat]
}