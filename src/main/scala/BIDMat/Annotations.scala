
package BIDMat
import java.lang.reflect.Method
import scala.util.Sorting._

class about(x:String) extends scala.annotation.Annotation {}

class Options {
  val ignore = List("what", "wait", "equals", "toString", "hashCode", "getClass", "notify", "notifyAll", "ignore")
  def what:Unit = {
    val a = for (meth <- this.getClass.getMethods 
    		if (!meth.getName.contains("$eq") && !meth.getName.contains("$$methOrdering") && !ignore.contains(meth.getName))) yield meth
    implicit object methOrdering extends Ordering[Method] { 
    	def compare(x: Method, y: Method) = x.getName.toLowerCase.compareTo(y.getName.toLowerCase) 
    }
    stableSort(a)
    println("Type        Value       Field Name")
    println("====        =====       ==========")
    for (meth <- a) {
      val cname = meth.getDeclaringClass().getName
      val params = meth.getParameterTypes();
      if (params.length == 0) {
        val ref = meth.invoke(this)
        val valstr = if (ref != null) ref.toString else "null"
        println("%-10s  %-10s  %s" format (meth.getReturnType.getSimpleName, valstr, meth.getName))
      }
    }    
  }
}