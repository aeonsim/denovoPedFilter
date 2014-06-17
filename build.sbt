import AssemblyKeys._

assemblySettings

jarName in assembly := "pedigreeFilter.jar"

name := "advFilter"

version := "1.0"

scalaVersion := "2.11.1"

scalacOptions ++= Seq(
      "-optimize"
    )
