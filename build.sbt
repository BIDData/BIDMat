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

credentials += Credentials(Path.userHome / ".ivy2" / ".credentials")

parallelExecution in Test := false

javacOptions ++= Seq("-source", "1.7", "-target", "1.7")

scalacOptions ++= Seq("-feature","-deprecation","-target:jvm-1.7")

initialCommands := scala.io.Source.fromFile("lib/bidmat_init.scala").getLines.mkString("\n")

javaOptions += "-Xmx4g"
