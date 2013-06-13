import scala.collection.mutable.HashMap
import java.io._

val in = new BufferedReader(new FileReader("LIC_ELs_results.txt"))

var current = in.readLine.split("\t")

var markerPos = new HashMap[Int,String]
var markers = new Array[HashMap[String,Int]](current.size -1)

for (i <- (1 to current.size - 1)){
//current(i).split("-").apply(0)
markerPos +=  i -> current(i)
}

while (in.ready){
current = in.readLine.split("\t")
for (i <- (0 to (current.size - 2))){

if (current(i + 1) != ""){
if(markers(i) != null){
if (markers(i).contains(current(i +1))){
markers(i)(current(i+1)) += 1
} else {
markers(i) += current(i+1) -> 1
}//eif 3
} else {

markers(i) = new HashMap[String,Int]
markers(i) += current(i+1) -> 1

}//eif2
}//eif1
}//efor
}//ewhile

for (i <- (0 to markers.size -1)){
print(markerPos(i + 1))
if (markers(i) != null) {
if (markers(i).toArray.size == 1) print("\tXX\t0\tXx\t0")
if (markers(i).toArray.size == 2) print("\tXX\t0")
markers(i).toArray.sortBy(s => s._2).foreach(s => print( "\t" + s._1 + "\t" + s._2))
}
print("\n")
}

