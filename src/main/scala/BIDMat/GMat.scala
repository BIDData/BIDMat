
package BIDMat

import scala.concurrent.Future
import scala.concurrent.ExecutionContext.Implicits.global
import scala.util.hashing.MurmurHash3
import edu.berkeley.bid.MurmurHash3.MurmurHash3_x64_64
import java.io.ObjectOutputStream
import java.io.ObjectInputStream

object GMat {
  
  object BinOp {
  	val op_add=0
  	val op_sub=1
	  val op_mul=2
	  val op_div=3
	  val op_gt=4
	  val op_lt=5
	  val op_eq=6
	  val op_ge=7
	  val op_le=8
	  val op_ne=9
	  val op_max=10
	  val op_min=11
	  val op_atan2=12
	  val op_pow=13
	  val op_ifpos=14
  }  
  
  object TransF {
    val abs=0
    val exp=1 
    val expm1=2
    val sqrt=3
    val ln=4
    val log10=5 
    val log1p=6
    val cos=7
    val sin=8
    val tan=9
    val cosh=10 
    val sinh=11 
    val tanh=12 
    val acos=13 
    val asin=14
    val atan=15 
    val acosh=16 
    val asinh=17 
    val atanh=18
    val erf=19
    val erfinv=20 
    val erfc=21
    val erfcinv=22 
    val gammaln=23
    val gamma=24
    val ceil=25
    val floor=26
    val round=27
    val trunc=28
    val sign=29
    val exppsi=34
    val normcdf=35
    val normcdfinv=36
    val logistic=37    
  }
  
  var multYieldSize = 1e6f;
  
  object TransF2 {
    val atan2=0
    val pow=1 
  }  
}







