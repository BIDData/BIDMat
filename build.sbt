name := "BIDMat"

version := "1.1.0"

organization := "edu.berkeley.bid"

scalaVersion := "2.11.2"

val jcudaVersion = "0.7.0a"

val OSmap = List(("windows", "windows"), ("linux", "linux"), ("Mac", "apple"))

val ARCHmap = List(("amd64", "x86_64"), ("x86_64", "x86_64"));

val OS_ = System.getProperty( "os.name" ).toLowerCase
val OS = OSmap.flatMap( v => if (OS_.contains(v._1)) List(v._2) else List()).head;
val ARCH_ = System.getProperty( "os.arch" ).toLowerCase
val ARCH = ARCHmap.flatMap( v => if (ARCH_.contains(v._1)) List(v._2) else List()).head;

artifactName := { (sv: ScalaVersion, module: ModuleID, artifact: Artifact) =>
  "../../BIDMat.jar"
}

resolvers ++= Seq(
  "Scala Tools Snapshots" at "http://scala-tools.org/repo-snapshots/",
  "Scala Mirror" at "https://oss.sonatype.org/content/repositories/releases/"
)

resolvers += Resolver.bintrayRepo("biddata", "BIDData")

libraryDependencies <<= (scalaVersion, libraryDependencies) { (sv, deps) =>
  deps :+ ("org.scala-lang" % "scala-compiler" % sv)
}

//libraryDependencies += "jline" % "jline" % "2.11"
libraryDependencies += "jline" % "jline" % "2.10"

libraryDependencies += "org.apache.commons" % "commons-math3" % "3.2"

//libraryDependencies += "org.scala-lang.modules" %% "scala-parser-combinators" % "1.0.3"

libraryDependencies += "org.scalatest" %% "scalatest" % "2.2.4" % "test"

libraryDependencies += "org.scalacheck" %% "scalacheck" % "1.11.6" % "test"

libraryDependencies += "junit" % "junit" % "4.11" % "test"

libraryDependencies += "net.jpountz.lz4" % "lz4" % "1.3"

libraryDependencies += "com.cedarsoftware" % "json-io" % "4.5.0"

libraryDependencies += "org.jfree" % "jfreechart" % "1.0.19"

libraryDependencies += "org.jfree" % "jcommon" % "1.0.23"

libraryDependencies += "com.jsuereth" %% "scala-arm" % "1.4"

libraryDependencies += "org.jocl" % "jocl" % "2.0.0"

libraryDependencies += "ptplot" % "ptplot" % "1.0"

libraryDependencies += "jhdf5" % "jhdf5" % "3.2.1"
libraryDependencies += "jhdf5" % "jhdf5" % "3.2.1" classifier (OS + "-" + ARCH)

libraryDependencies += "jcuda" % "jcuda" % jcudaVersion
libraryDependencies += "jcuda" % "jcuda" % jcudaVersion classifier (OS + "-" + ARCH)

libraryDependencies += "jcuda" % "jcublas" % jcudaVersion
libraryDependencies += "jcuda" % "jcublas" % jcudaVersion classifier (OS + "-" + ARCH)

libraryDependencies += "jcuda" % "jcurand" % jcudaVersion
libraryDependencies += "jcuda" % "jcurand" % jcudaVersion classifier (OS + "-" + ARCH)

libraryDependencies += "jcuda" % "jcusparse" % jcudaVersion
libraryDependencies += "jcuda" % "jcusparse" % jcudaVersion classifier (OS + "-" + ARCH)

credentials += Credentials(Path.userHome / ".ivy2" / ".credentials")

parallelExecution in Test := false

javacOptions ++= Seq("-source", "1.7", "-target", "1.7")

scalacOptions ++= Seq("-feature","-deprecation","-target:jvm-1.7")

initialCommands := scala.io.Source.fromFile("lib/bidmat_init.scala").getLines.mkString("\n")

javaOptions += "-Xmx8g"

