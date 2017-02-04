package BIDMat

import Mat._
import MatFunctions._
import org.scalatest._;
import org.scalatest.junit._;
import org.scalatest.prop._;
import org.junit.runner.RunWith

@RunWith(classOf[JUnitRunner])
class DMatTest extends BIDMatSpec {
      Mat.checkMKL(true)

}
