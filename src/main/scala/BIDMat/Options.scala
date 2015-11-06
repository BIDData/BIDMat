package BIDMat
import java.lang.reflect.Method
import scala.util.Sorting._

class about(x: String) extends scala.annotation.Annotation {}

trait Opts {
  def what: Unit = {
    val ignore = List("what", "wait", "equals", "toString", "hashCode", "getClass", "notify", "notifyAll", "ignore", "copyFrom")

    val a = for (
      meth <- this.getClass.getMethods if (!meth.getName.contains("$eq") && !meth.getName.contains("$$methOrdering") && !ignore.contains(meth.getName))
    ) yield meth;
    implicit object methOrdering extends Ordering[Method] {
      def compare(x: Method, y: Method) = x.getName.toLowerCase.compareTo(y.getName.toLowerCase);
    }
    stableSort(a);
    println("Option Name       Type          Value");
    println("===========       ====          =====");
    for (meth <- a) {
      val cname = meth.getDeclaringClass().getName
      val params = meth.getParameterTypes();
      if (params.length == 0) {
        val ref = meth.invoke(this)
        val valstring = if (ref != null) ref.toString else "null"
        val valstr = valstring.substring(0, math.min(20, valstring.length)) + (if (valstring.length > 20) "..." else "")
        println("%-16s  %-12s  %s" format (meth.getName, meth.getReturnType.getSimpleName, valstr))
      }
    }
  }

  def copyFrom(that: Opts): Unit = {
    val ignore = List("what", "wait", "equals", "toString", "hashCode", "getClass", "notify", "notifyAll", "ignore", "copyFrom");
    for (setter <- this.getClass.getMethods if (setter.getName.contains("$eq") && !ignore.contains(setter.getName))) {
      val setname = setter.getName;
      val getname = setname.substring(0, setname.length - 4); // getter name = remove _= suffix from a setter method
      try {
        val getter = that.getClass.getMethod(getname); // Try to find a method with that name in the other class
        val value = getter.invoke(that); // Get the value from that
        setter.invoke(this, value); // Copy into this
      } catch {
        case _: NoSuchMethodException => {} // OK if not found. "that" could be a superclass
      }
    }
  }
}
