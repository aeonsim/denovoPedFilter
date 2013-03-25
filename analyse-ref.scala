import java.io._
import scala.collection.mutable.HashMap

val input = new BufferedReader(new FileReader("umd31MT.fa"))

val nucStats = new HashMap[String,Int]
val diNucStats = new HashMap[String,Int]
var lastBase = ""
val bases = Array("A","G","C","T","N")
    
    def permu(str: String, str1: String): List[String] = {
      val result: List[String] = for (i <- str.toList; j <- str1.toList) yield i.toString.concat(j.toString)
      result
    }
    
val diNuc = permu("AGCTN","AGCTN")

for (i <- bases){
nucStats += i -> 0
}

for (i <- diNuc){
diNucStats += i -> 0
}

while (input.ready) {
val line = input.readLine
var pos = 0
if (line(0) != '>'){
while (pos < line.size){
nucStats(line(pos).toString.toUpperCase) += 1
if (lastBase != ""){
diNucStats(lastBase + line(pos).toString.toUpperCase) += 1
lastBase = line(pos).toString.toUpperCase
} else {
lastBase = line(pos).toString.toUpperCase
}//Eif Lastbase
pos += 1
}//E while Line
}//Eif line(0)
}//Ewhile input loop

nucStats.foreach(s => println(s._1 + "\t" + s._2))
diNucStats.foreach(s => println(s._1 + "\t" + s._2))