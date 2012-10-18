/* Copyright (c) 2012, Regents of the University of California                     */
/* All rights reserved.                                                            */

/* Redistribution and use in source and binary forms, with or without              */
/* modification, are permitted provided that the following conditions are met:     */
/*     * Redistributions of source code must retain the above copyright            */
/*       notice, this list of conditions and the following disclaimer.             */
/*     * Redistributions in binary form must reproduce the above copyright         */
/*       notice, this list of conditions and the following disclaimer in the       */
/*       documentation and/or other materials provided with the distribution.      */
/*     * Neither the name of the <organization> nor the                            */
/*       names of its contributors may be used to endorse or promote products      */
/*       derived from this software without specific prior written permission.     */

/* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND */
/* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED   */
/* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE          */
/* DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY              */
/* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES      */
/* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;    */
/* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND     */
/* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT      */
/* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS   */
/* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                    */

package BIDMat
import scala.util.parsing.combinator._

class MKLdeclarationsParser extends JavaTokenParsers { 
  var otype:Int = 3
  var prefix:String = "LAPACK"
  def fdecl: Parser[String] = ftype~ident~"("~decls~")"~""";*""".r ^^ { 
    case a~b~c~d~e~f => {
      val b0 = b.split("_"); val b1 = b0(b0.length-1)
      val toreturn:Boolean = (a != fvoid)
    	otype match {
    	case 0 => "JNIEXPORT "+a("")+"JNICALL Java_edu_berkeley_bid_"+prefix+"_"+b1+
    	"\n(JNIEnv * env, jobject calling_obj, "+d+") {\n"
    	case 1 => d
    	case 2 => "    "+(if (toreturn) a("")+"retval = " else "") +b+c+d+e+";\n"
    	case 3 => d+(if (toreturn) "    return retval;\n}\n" else "}\n")
    	case 4 => "  public static native "+a(b1)+c+d+e+";\n"
    	}
    }
  } 
  
  def decls: Parser[String] = repsep(adecl, ",") ^^ {
    case a:List[String] => {
      var ll = ""
      	for (el <- a) {
      	  ll = ll + el
      	  if (otype == 0 || otype == 2 || otype == 4) {
      	    if (el != a.last) ll = ll+", "
      	  }
      	}
      ll
      }
    }
  
  def adecl: Parser[String] = ftype~ident ^^ { case a~b => a(b) }
  
  def ftype: Parser[String=>String] = fvoidp | fvoid | flogicalp | flogical | fintp | fint | fuintp | fuint | fchar |
    ffloatp | ffloat | ffcomplexp | ffcomplex | fdoublep | fdouble | fdcomplexp | fdcomplex 
   
  def sstring(y:String, wrap:Int, p1:String, p2:String, p3:String, p4:String) = {
    otype match {
      case 0 => p1 + " " + (if (wrap > 0) "j_"+y else y)
      case 1 => if (wrap == 1) "    "+p2+" "+y+" = (*env)->GetPrimitiveArrayCritical(env, j_"+y+", JNI_FALSE);\n" else 
      	if (wrap == 2) "    char * "+y+" = (char *)(*env)->GetStringUTFChars(env, j_"+y+", JNI_FALSE);\n" else ""
      case 2 => if (y.length == 0) p2+" " else  p3+y
      case 3 => if (wrap == 1) "    (*env)->ReleasePrimitiveArrayCritical(env, j_"+y+", "+y+", 0);\n" else 
        if (wrap == 2) 	"    (*env)->ReleaseStringUTFChars(env, j_"+y+", "+y+");\n" else ""
      case 4 => p4+" "+y
    }
  }
  
  def fvoid: Parser[String=>String] = "void" ^^ (x => (y=>
    sstring(y, 0, "void", "void", "", "void")))
    
  def fenum: Parser[String=>String] = ("const"~"enum"~ident | "enum"~ident)  ^^ {
    case a~b~c => (y => sstring(y, 0, "jint", "jint", "("+c+")", "int"))
    case b~c => (y => sstring(y, 0, "jint", "jint", "("+c+")", "int"))
  }
  
  def flogical: Parser[String=>String] = ("const"~"lapack_logical" | "lapack_logical") ^^ (x => (y=>
    sstring(y, 1, "int", "jint *", "(lapack_logical)", "int")))

  def fchar: Parser[String=>String] = ("const"~"char" | "char") ^^ (x => (y=> 
    sstring(y, 2, "jstring", "char *", "*", "String")))
    
  def fint: Parser[String=>String] = ("const"~"int" | "int" | "const"~"lapack_int" | "lapack_int" | "MKL_INT") ^^ (x => (y=> 
    sstring(y, 0, "jint", "jint", "", "int")))
  
  def fuint: Parser[String=>String] = ("const"~"MKL_UINT" | "MKL_UINT"| "const"~"unsigned"~"int" | "unsigned"~"int")^^ (x => (y =>
    sstring(y, 0, "jint", "jint", "(unsigned int)", "int")))
  
  def ffloat: Parser[String=>String] = ("const"~"float" | "float" | "lapack_float") ^^ (x => (y =>
    sstring(y, 0, "jfloat", "jfloat", "", "float")))
    
  def fdouble: Parser[String=>String] = ("const"~"double" | "double" | "lapack_double") ^^ (x => (y =>
    sstring(y, 0, "jdouble", "jdouble", "", "double")))
    
  def ffcomplex: Parser[String=>String] = ("const"~"lapack_complex_float" | "lapack_complex_float" | "MKL_complex8") ^^ (x => (y =>
    sstring(y, 1, "jfloatArray", "jfloat *", "*", "float []")))
    
  def fdcomplex: Parser[String=>String] = ("const"~"lapack_complex_double" | "lapack_complex_double" | "MKL_complex16") ^^ (x => (y =>
    sstring(y, 1, "jdoubleArray", "jdouble *", "*", "double []")))
  
  def fvoidp: Parser[String=>String] = ("const"~"void"~"*" | "void"~"*") ^^ (x => (y =>
    sstring(y, 1, "jfloatArray", "jfloat *", "(void *)", "float []")))
  
  def flogicalp: Parser[String=>String] = ("const"~"lapack_logical"~"*" | "lapack_logical"~"*") ^^ (x => (y=>
    sstring(y, 1, "jintArray", "jint *", "(lapack_logical *)", "int []")))

  def fintp:  Parser[String=>String] = ("const"~"int"~"*" | "int"~"*" | "MKL_INT"~"*" | "lapack_int"~"*") ^^ (x => (y=> 
    sstring(y, 1, "jintArray", "jint *", "", "int []")))
     
  def fuintp: Parser[String=>String] = ("const"~"MKL_UINT"~"*" | "MKL_UINT"~"*" | "const"~"unsigned"~"int"~"*" | "unsigned"~"int"~"*") ^^ (x => (y =>
    sstring(y, 1, "jintArray ", "jint *", "", "int []")))

  def ffloatp: Parser[String=>String] = ("const"~"float"~"*" | "float"~"*") ^^ (x => (y=> 
    sstring(y, 1, "jfloatArray", "jfloat *", "", "float []")))
  
  def fdoublep: Parser[String=>String] = ("const"~"double"~"*" | "double"~"*") ^^ (x => (y=> 
    sstring(y, 1, "jdoubleArray", "jdouble *", "", "double []")))
  
  def ffcomplexp: Parser[String=>String] = ("const"~"lapack_complex_float"~"*" | "lapack_complex_float"~"*" | "MKL_complex8"~"*") ^^ {
    case x0~x1~x2 => (y => sstring(y, 1, "jfloatArray", "jfloat *", "("+x1+" *)", "float []"))
    case x1~x2 => (y => sstring(y, 1, "jfloatArray", "jfloat *", "("+x1+" *)", "float []"))
    }
    
  def fdcomplexp: Parser[String=>String] = ("const"~"lapack_complex_double"~"*" | "lapack_complex_double"~"*" | "MKL_complex16"~"*") ^^ {
    case x0~x1~x2 => (y=> sstring(y, 1, "jdoubleArray", "jdouble *", "("+x1+" *)", "double []"))
    case x1~x2 => (y=> sstring(y, 1, "jdoubleArray", "jdouble *", "("+x1+" *)", "double []"))
    }

}

object translateTester { 
  def main(args: Array[String]) {
    val p = new MKLdeclarationsParser
    val s = scala.io.Source.fromFile(args(0))
    if (args.length == 1) {
    	s.getLines.foreach((line) => {
    		if (line.length > 1) {
    			p.otype = 0
    			println(p.parseAll(p.fdecl, line).get)
    			p.otype = 1
    			println(p.parseAll(p.fdecl, line).get)
    			p.otype = 2
    			println(p.parseAll(p.fdecl, line).get)
    			p.otype = 3
    			println(p.parseAll(p.fdecl, line).get)
    		}
    	})
    } else {
    	s.getLines.foreach((line) => {
    		if (line.length > 1) {
    			p.otype = 4	  
    			println(p.parseAll(p.fdecl, line).get)
    		}
    	})
    }
  }
}



object Translators { 


}
