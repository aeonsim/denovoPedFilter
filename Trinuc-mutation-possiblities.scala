
import scala.collection.mutable.HashMap
import scala.collection.mutable.HashSet

def permu3(str: String, str1: String, str2: String): List[String] = {
           val result: List[String] = for (i <- str.toList; j <- str1.toList; k <- str2.toList) yield i.toString.concat(j.toString.concat(k.toString))
            result
          }
          
val alts = HashMap('A' -> 'T','T' -> 'A', 'G' -> 'C', 'C' -> 'G', 'N' -> 'N')

val nucs = permu3("AGCT","AGCT","AGCT").toArray
val bases = Array('A','G','C','T')
val used = new HashSet[String]

for (nuc <- nucs){
for (b <- bases){
var working =  nuc
if (working(1) != b){
val mut = s" ${working(0)}$b${working(2)}"
val revNuc = s"${alts(working(2))}${alts(working(1))}${alts(working(0))}"
val revMut =  s"${alts(working(2))}${alts(b)}${alts(working(0))}"
//if (!(used.contains(nuc)) && !(used.contains(mut)) && !(used.contains(revNuc)) && !(used.contains(revMut))){
if (!(used.contains(nuc)) || !(used.contains(revNuc))){
//if ( !(used.contains(revNuc))){
//used += mut
used += revNuc
//used += revMut
println(s"$nuc -> ${working(0)}$b${working(2)} & ${alts(working(2))}${alts(working(1))}${alts(working(0))} -> ${alts(working(2))}${alts(b)}${alts(working(0))}")
}
}
}
used += nuc
}
