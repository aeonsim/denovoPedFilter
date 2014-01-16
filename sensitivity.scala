import java.io._
import scala.collection.mutable.HashMap
import scala.collection.mutable.HashSet
import java.util.zip._

object sensVCF {

def main(args: Array[String]): Unit = {

if (args.size < 3){
println("snsVCF gsr-genotypes.txt gsr-markers.txt file.vcf optional{minDP maxDP minGQ}")
System.exit(1)
}
val gInput = new BufferedReader(new FileReader(args(0)))
val gDetail =  new BufferedReader(new FileReader(args(1)))
val vcf = new BufferedReader(new FileReader(args(2)))

var minDP: Int = if (args.size >= 4) args(3).toInt else 0
var maxDP: Int = if (args.size >= 5) args(4).toInt else 10000
var minGQ: Int = if (args.size >= 6) args(5).toInt else 0


var vcfAnimals = new HashMap[String,Int]
var vcfLine = vcf.readLine
while ((vcfLine(0) == '#')&(vcfLine(1) == '#')){
vcfLine = vcf.readLine
}

if ((vcfLine(0) == '#')&(vcfLine(1) == 'C')){
var vcfHeader = vcfLine.split("\t")
for (i <- 9 to (vcfHeader.size - 1)){
vcfAnimals += vcfHeader(i) -> i
}
}



var markerSet = new HashMap[String,Tuple3[Int, String, String]]
var markerPos = 1
while (gDetail.ready){
val marker = gDetail.readLine.split(" ")
markerSet += ("Chr" + marker(0) + ":" + marker(2)) -> (markerPos,marker(3),marker(4))
markerPos += 2
}
gDetail.close

val gts = new HashSet[Array[String]]
while (gInput.ready){
gts += gInput.readLine.split(" ")
}
gInput.close
var animalStats = new HashMap[String,Array[Int]]

//Create something to store our stats
for (ans <- gts.toArray){
animalStats += ans(0) -> Array(0,0,0,0,0,0,0,0,0)
}

// Remove any animals from the GSR list that aren't present in the VCF
for (ans <- gts.toArray){
if (!(vcfAnimals.contains(ans(0)))) {
gts.remove(ans)
//println("remove " + ans(0))
} //else println("keep " + ans(0))
}

// Want to store: 0 number missed, 1# agree, 2het -> hom, 3Hom -> Het, 4Het to Hom DP, 5Hom -> Het DP

while (vcf.ready){
val vcfRecord = vcf.readLine.split("\t")
val pos = vcfRecord(0) + ":" + vcfRecord(1)
if (markerSet.contains(pos)){
val curMarker = markerSet(pos)

for (animal <- gts.toArray){
if (vcfAnimals.contains(animal(0))){
val geno = vcfRecord(vcfAnimals(animal(0))).split(":")

// Set index positions of Records we are interested in

val DP = vcfRecord(8).split(":").indexOf("DP")
val GQ = vcfRecord(8).split(":").indexOf("GQ")
val GT = vcfRecord(8).split(":").indexOf("GT")

//Look at Markers
if (animal(curMarker._1) != "0"){
// Make sure we are not looking at . or ./.
if (geno.size > 1) {
// Apply filters
if ((minDP <= geno(DP).toInt) & (maxDP >= geno(DP).toInt) & (minGQ <= geno(GQ).toInt)){
//Evaluate if the GSR values are the same for a marker it's Hom (either Ref/Ref or Alt/Alt) note this when reporting
//println("DP and GQ are:" + geno(DP) + "\t" + geno(GQ))
val gtAllele = if ( animal(curMarker._1) == animal(curMarker._1 + 1)) "hom" else "het"
val vcfAllele = if (geno(GT).apply(0) == geno(GT).apply(2)) "hom" else "het"
// Check Values and Record stats
if (gtAllele == vcfAllele) {
if (gtAllele == "het") animalStats(animal(0)).apply(8) += 1
 animalStats(animal(0)).apply(1) += 1
 animalStats(animal(0)).apply(6) += geno(DP).toInt
} else {
if (gtAllele == "het") {
animalStats(animal(0)).apply(2) += 1
animalStats(animal(0)).apply(4) += geno(DP).toInt
} else {
if (gtAllele ==  "hom") {
animalStats(animal(0)).apply(3) += 1
animalStats(animal(0)).apply(5) += geno(DP).toInt
} else {
 println("ERROR: " + gtAllele + " " +  animal(curMarker._1) +  animal(curMarker._1 + 1) )

}
}

}

} //eIf for DP GQ checks
} else { //Else for VCF ./.
animalStats(animal(0)).apply(7) += 1 
}

} else {//eIf
//Increase Missing count
animalStats(animal(0)).apply(0) += 1 
}
}// eif Animal Exists in VCF
} //efor

} //eif

} //eWhile

vcf.close

println("minDP,maxDP,minGQ\tAnimal\tSNPChip SNPs\tSC no call\t%\tVCF no call\t%\tVCF missing\t%\tSC Agree\t%\tHet Agree\tSC Het->Hom\t%\tSC Hom->Het\t%\tMean DP Agree\tMean DP Het->Hom\tMean DP Hom->Het")

for (animals <- gts){
val total = ((animals.size - 1)/2 )*1.0
//println ("Animal: " + animals(0) + "\tSNPChip SNPS: " + ((animals.size - 1)/2 )) 
val stats = animalStats(animals(0))
println( minDP + "," + maxDP + "," + minGQ + "\t"  + animals(0) + "\t" + ((animals.size - 1)/2 ) + "\t" + stats(0) + "\t" + (stats(0)/total*100) + "\t" +
stats(7) + "\t" + (stats(7)/total*100) + "\t" + (total - (stats(0) + stats(1) + stats(2) + stats(3) + stats(7))) + "\t" + ((total - (stats(0) + stats(1) + stats(2) + stats(3) + stats(7)))/total*100) + "\t" +
stats(1) + "\t" + (stats(1)/total*100) + "\t" + stats(8) + "\t" +  stats(2) + "\t" + (stats(2)/total*100) + "\t" + stats(3) +  "\t" + (stats(3)/total*100) + "\t" +
(if (stats(1) == 0) "-" else (stats(6)/stats(1))) + 
"\t" + (if (stats(2) == 0) "-" else (stats(4)/stats(2))) + 
"\t" + (if (stats(3) == 0) "-" else (stats(5)/stats(3)))
)
}
}
}