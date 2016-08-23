name := "BIDMat"

version := "1.1.0"

organization := "edu.berkeley.bid"

scalaVersion := "2.11.2"

artifactName := { (sv: ScalaVersion, module: ModuleID, artifact: Artifact) =>
  "../../BIDMat.jar"
}

resolvers ++= Seq(
  "Scala Tools Snapshots" at "http://scala-tools.org/repo-snapshots/",
  "Scala Mirror" at "https://oss.sonatype.org/content/repositories/releases/"
)

libraryDependencies <<= (scalaVersion, libraryDependencies) { (sv, deps) =>
  deps :+ ("org.scala-lang" % "scala-compiler" % sv)
}

//libraryDependencies += "jline" % "jline" % "2.11"
libraryDependencies += "jline" % "jline" % "2.10"

libraryDependencies += "org.apache.commons" % "commons-math3" % "3.2"

//libraryDependencies += "org.scala-lang.modules" %% "scala-parser-combinators" % "1.0.2"

libraryDependencies += "org.scalatest" %% "scalatest" % "2.2.4" % "test"

libraryDependencies += "org.scalacheck" %% "scalacheck" % "1.11.6" % "test"

libraryDependencies += "junit" % "junit" % "4.11" % "test"

libraryDependencies += "net.jpountz.lz4" % "lz4" % "1.3"

libraryDependencies += "com.cedarsoftware" % "json-io" % "2.2.29"

libraryDependencies += "com.jsuereth" %% "scala-arm" % "1.4"

//libraryDependencies += "org.scala-saddle" % "jhdf5" % "2.9"

credentials += Credentials(Path.userHome / ".ivy2" / ".credentials")

parallelExecution in Test := false

javacOptions ++= Seq("-source", "1.7", "-target", "1.7")

scalacOptions ++= Seq("-feature","-deprecation","-target:jvm-1.7")

initialCommands := scala.io.Source.fromFile("lib/bidmat_init.scala").getLines.mkString("\n")

javaOptions += "-Xmx4g"
