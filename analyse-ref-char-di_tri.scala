import java.io._
import scala.collection.mutable.HashMap

object anFasta{

def main(args:Array[String]) : Unit = {

val input = new BufferedReader(new FileReader (args(0)))

val nucStats = new HashMap[String,Int]
val diNucStats = new HashMap[String,Int]
val rdiNucStats = new HashMap[String,Int]
val triNucStats = new HashMap[String,Int]
val rtriNucStats = new HashMap[String,Int]


var lastBase = ""
val bases = Array("A","G","C","T","N")
val alts = HashMap('A' -> 'T','T' -> 'A', 'G' -> 'C', 'C' -> 'G', 'N' -> 'N')
    
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
rdiNucStats += i -> 0
}

for (i <- triNuc){
triNucStats += i -> 0
rtriNucStats += i -> 0
}

var cur, chr1, chr2, chr3, chr4, chr5 = '#'

while (input.ready){
cur = input.read.toChar
if (alts.contains(cur)) nucStats(cur.toString) += 1
if (cur == '\n') cur = input.read.toChar
if (cur == '>' || cur == ' '){
while (cur != '\n') {
cur = input.read.toChar
//\println(cur)
}
chr1 = '#'
chr2 = '#'
chr3 = '#'
chr4 = '#'
chr5 = '#'
} else {
chr5 = chr4
chr4 = chr3
chr3 = chr2
chr2 = chr1
chr1 = cur
if(alts.contains(chr2) && alts.contains(chr1)) {
diNucStats(chr2.toString + chr1.toString) += 1
}
if (alts.contains(chr3) && alts.contains(chr2) && alts.contains(chr1)) {
rdiNucStats(alts(chr2).toString + alts(chr3).toString) += 1
triNucStats(chr3.toString + chr2.toString + chr1.toString) += 1
}
if (alts.contains(chr5)){
rtriNucStats(alts(chr3).toString + alts(chr4).toString + alts(chr5).toString) += 1
}
}
}
input.close
println("Base composition")
nucStats.foreach(s => println(s._1 + "\t" + s._2))
println("Forward (5' -> 3') diNucs")
diNucStats.foreach(s => println(s._1 + "\t" + s._2))
println("Reverse Strand 5' -> 3' diNucs")
rdiNucStats.foreach(s => println(s._1 + "\t" + s._2))
println("Forward 5' -> 3' TriNucs")
triNucStats.foreach(s => println(s._1 + "\t" + s._2))
println("Reverse Strand 5' -> 3' TriNucs")
rtriNucStats.foreach(s => println(s._1 + "\t" + s._2))

}
}