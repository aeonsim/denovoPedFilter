import java.io._
import scala.collection.mutable.HashMap

val input = new BufferedReader(new FileReader("umd31MT.fa"))

val nucStats = new HashMap[String,Int]
val diNucStats = new HashMap[String,Int]
val revdiNucStats = new HashMap[String,Int]
val triNucStats = new HashMap[String,Int]
val rtriNucStats = new HashMap[String,Int]


var lastBase = ""
val bases = Array("A","G","C","T","N")
val alts = HashMap("A" -> "T","T" -> "A", "G" -> "C", "C" -> "G", "N" -> "N")
    
    def permu(str: String, str1: String): List[String] = {
      val result: List[String] = for (i <- str.toList; j <- str1.toList) yield i.toString.concat(j.toString)
      result
    }
    
	def permu3(str: String, str1: String, str2: String): List[String] = {
      val result: List[String] = for (i <- str.toList; j <- str1.toList; k <- str2.toList) yield i.toString.concat(j.toString.concat(k.toString))
      result
    }
    
val diNuc = permu("AGCTN","AGCTN")
val triNuc = permu3("AGCTN","AGCTN","AGCTN")

for (i <- bases){
nucStats += i -> 0
}

for (i <- diNuc){
diNucStats += i -> 0
revdiNucStats += i -> 0
}

for (i <- triNuc){
triNucStats += i -> 0
rtriNucStats += i -> 0
}

var lastLine = ""


while (input.ready) {
val line = input.readLine
var pos = 0
//println(line)
if ((line.size != 0) && (line(0) != '>')){
while (pos < line.size){
nucStats(line(pos).toString.toUpperCase) += 1 //Count base simple
if (lastBase == ""){ //Not first character
if (pos == (line.size -1)){
 revdiNucStats(line(pos).toString.toUpperCase + line(pos - 1).toString.toUpperCase) += 1
 lastBase = line(pos).toString.toUpperCase
 } else {
 	if (pos == 0){
 	 diNucStats(line(pos).toString.toUpperCase + line(pos + 1).toString.toUpperCase) += 1
 	 } else {
 	  revdiNucStats(line(pos).toString.toUpperCase + line(pos - 1).toString.toUpperCase) += 1
 	  diNucStats(line(pos).toString.toUpperCase + line(pos + 1).toString.toUpperCase) += 1
 	 }
 	
 	}
} else { // LastBase Else
	revdiNucStats(line(pos).toString.toUpperCase + lastBase) += 1
	diNucStats(line(pos).toString.toUpperCase + line(pos + 1).toString.toUpperCase) += 1
	diNucStats(lastBase + line(pos).toString.toUpperCase) += 1
	lastBase = ""
}//Eif Lastbase

if (lastLine == "") {

} else {

} //LastLine EIF

pos += 1
}//E while Line
} else {//Eif line(0)
lastBase = ""
}
}//Ewhile input loop

nucStats.foreach(s => println(s._1 + "\t" + s._2))
println("Forward diNucs")
diNucStats.foreach(s => println(s._1 + "\t" + s._2))
println("Reverse diNucs")
revdiNucStats.foreach(s => println(s._1 + "\t" + s._2))

