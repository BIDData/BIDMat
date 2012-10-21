package BIDMat

import MatFunctions._
import SciFunctions._
import Learner._

object TestLearner {
  def main(args: Array[String]): Unit = {
  	val dirname = "d:\\sentiment\\sorted_data\\books\\parts\\"
  	val revtrain:SDMat = load(dirname+"part1.mat", "revtrain")
  	val rt = SMat(revtrain) 
  	val scrtrain:IMat = load(dirname+"part1.mat", "scrtrain")
  	val st = FMat(scrtrain).t > 4
  	val learner = Learner(rt, st, new LogisticModel, new ADAGradUpdater)()
  	learner.run
  }
}
