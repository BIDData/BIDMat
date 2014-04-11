
name := "BIDMat"

version := "0.1.0"

organization := "edu.berkeley.bid"

scalaVersion := "2.10.3"

resolvers ++= Seq(
  "Scala Tools Snapshots" at "http://scala-tools.org/repo-snapshots/",
  "Scala Mirror" at "https://oss.sonatype.org/content/repositories/releases/"
)

libraryDependencies <<= (scalaVersion, libraryDependencies) { (sv, deps) =>
  deps :+ ("org.scala-lang" % "scala-compiler" % sv)
}

libraryDependencies += "org.scala-lang" % "jline" % "2.10.3"

libraryDependencies += "org.apache.commons" % "commons-math3" % "3.2"

libraryDependencies += "net.jpountz.lz4" % "lz4" % "1.2.0"

libraryDependencies += "org.scala-saddle" % "jhdf5" % "2.9"

libraryDependencies += "org.opt4j" % "opt4j-viewer" % "3.0.1"

libraryDependencies += "org.scalatest" %% "scalatest" % "2.0" % "test"

libraryDependencies += "org.scalacheck" %% "scalacheck" % "1.11.2" % "test"

libraryDependencies += "junit" % "junit" % "4.10" % "test"

credentials += Credentials(Path.userHome / ".ivy2" / ".credentials")

javacOptions ++= Seq("-source", "1.6", "-target", "1.6")

scalacOptions ++= Seq("-feature","-deprecation","-target:jvm-1.6")

initialCommands := scala.io.Source.fromFile("lib/bidmat_init.scala").getLines.mkString("\n")

javaOptions += "-Xmx12g"

