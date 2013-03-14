import java.io._
import scala.collection.mutable.HashMap
import scala.collection.mutable.HashSet
import java.util.zip._

object sensVCF {

def main(args: Array[String]): Unit = {

if (args.size < 3){
println("snsVCF gsr-genotypes.txt gsr-markers.txt file.vcf")
System.exit(1)
}
val gInput = new BufferedReader(new FileReader(args(0)))
val gDetail =  new BufferedReader(new FileReader(args(1)))
val vcf = new BufferedReader(new FileReader(args(2)))

//val gInput = new BufferedReader(new FileReader("chhar0-20130311-12230670/hd-gsr-genotypes.txt"))
//val gDetail =  new BufferedReader(new FileReader("chhar0-20130311-12230670/hd-marker-map.txt"))
//val vcf = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream("/Volumes/Travel/PhD/Denovo-Hunting/GATK-LIC-54-indv-VQSR-99.9.vcf.gz"))))
//val vcf = new BufferedReader(new FileReader("/Volumes/Travel/PhD/Denovo-Hunting/GATK-LIC-54-indv-VQSR-99.9.vcf"))


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

for (ans <- gts){
animalStats += ans(0) -> Array(0,0,0,0,0,0,0,0)
}

for (i <- gts){
if (!(vcfAnimals.contains(i(0)))) {
gts.remove(i)
}
}

// Want to store: 0 number missed, 1# agree, 2het -> hom, 3Hom -> Het, 4Het to Hom DP, 5Hom -> Het DP

while (vcf.ready){
val vcfRecord = vcf.readLine.split("\t")
val pos = vcfRecord(0) + ":" + vcfRecord(1)
if (markerSet.contains(pos)){
val curMarker = markerSet(pos)
for (animal <- gts){
if (vcfAnimals.contains(animal(0))){
val geno = vcfRecord(vcfAnimals(animal(0))).split(":")
val DP = vcfRecord(8).split(":").indexOf("DP")
//println(pos + "\t" + vcfRecord(vcfAnimals(animal(0))))
//Assume GT = 0 for now
//
if (animal(curMarker._1) != "0"){
if (geno.size > 1) {
val gtAllele = if ( animal(curMarker._1) == animal(curMarker._1 + 1)) "hom" else "het"
val vcfAllele = if (geno(0).apply(0) == geno(0).apply(2)) "hom" else "het"
if (gtAllele == vcfAllele) {
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

println("Animal\tSNPChip SNPs\tSC no call\t%\tVCF no call\t%\tVCF missing\t%\tSC Agree\t%\tSC Het->Hom\t%\tSC Hom->Het\t%\tMean DP Agree\tMean DP Het->Hom\tMean DP Hom->Het")

for (animals <- gts){
val total = ((animals.size - 1)/2 )*1.0
//println ("Animal: " + animals(0) + "\tSNPChip SNPS: " + ((animals.size - 1)/2 )) 
val stats = animalStats(animals(0))
println( animals(0) + "\t" + ((animals.size - 1)/2 ) + "\t" + stats(0) + "\t" + (stats(0)/total*100) + "\t" +
stats(7) + "\t" + (stats(7)/total*100) + "\t" + (total - (stats(0) + stats(1) + stats(2) + stats(3) + stats(7))) + "\t" + ((total - (stats(0) + stats(1) + stats(2) + stats(3) + stats(7)))/total*100) + "\t" +
stats(1) + "\t" + (stats(1)/total*100) + "\t" +  stats(2) + "\t" + (stats(2)/total*100) + "\t" + stats(3) +  "\t" + (stats(3)/total*100) + "\t" +
(stats(6)/stats(1)) + "\t" + (stats(4)/stats(2)) + "\t" + (stats(5)/stats(3)))

/*
println("Cat\tNumber\tPercentage")
println("SNPChip no Call:\t" + stats(0) + "\t" + (stats(0)/total*100))
println("VCF no Call:\t" + stats(7) + "\t" + (stats(7)/total*100))
println("VCF Missing:\t" + (total - (stats(0) + stats(1) + stats(2) + stats(3) + stats(7))) + "\t" + ((total - (stats(0) + stats(1) + stats(2) + stats(3) + stats(7)))/total*100) )
println("SNPChip Agree:\t" + stats(1) + "\t" + (stats(1)/total*100))
println("SNPChip Error Het -> Hom:\t" + stats(2) + "\t" + (stats(2)/total*100))
println("SNPChip Error Hom -> Het:\t" + stats(3) + "\t" + (stats(3)/total*100))
println("Mean DP Agree:\t" + stats(6) + "\t" + (stats(6)/stats(1)))
println("Mean DP Het->Hom:\t" + stats(4) + "\t" + (stats(4)/stats(2)))
println("Mean DP Hom->Het:\t" + stats(5) + "\t" + (stats(5)/stats(3)))
println("---------\t---------\t---------")
*/
}
}
}