import java.io._
import net.sf.samtools.util._
import scala.collection.mutable.HashMap

object genie{

var GQ, DP = -1

def main(args: Array[String]) : Unit = {

println("VCF sire dam proband children.list\n Gzipped VCF file Sire, Dam & Proband IDs, List of Children 1 per line\n")

val vcf = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream(args(0))))) ///args(0)))))
val children = new BufferedReader(new FileReader(args(4)))

val sireID = args(1) //"21930910" //
val damID = args(2) //"22328328" //
val proID = args(3) //"23934418" //
var haps = new HashMap[String,List[Char]]
//var prohaps = new HashMap[String,List[Char]]
//var chihaps = new HashMap[String,List[Char]]
var anPos = new HashMap[String,Int]

var chlID: List[String] = Nil

while (children.ready) chlID = children.readLine :: chlID


var cline = vcf.readLine.split("\t")

while (cline(0).apply(1) == '#') cline = vcf.readLine.split("\t")

for (idx <- 9 to (cline.size - 1)){
anPos += cline(idx) -> idx
}

var count = 0
// Add blank Haps
haps += (sireID +  "R") -> Nil
haps += (sireID +  "L") -> Nil
haps += (damID +  "R") -> Nil
haps += (damID +  "L") -> Nil
haps += (proID +  "R") -> Nil
haps += (proID +  "L") -> Nil

for (chl <- chlID){
haps += (chl +  "R") -> Nil
haps += (chl +  "L") -> Nil
}

var inher  = false
var fails,numInher = 0
var pfail = ""

while(vcf.ready){
 cline = vcf.readLine.split("\t")
 //skip multiAllelic sites (for now)
 setFields(cline(8))
 
 
if(cline(4).split(",").size > 1){
cline = vcf.readLine.split("\t")
setFields(cline(8))
}

//filter GQ & DP
if ((cline(anPos(sireID)).split(":").apply(GQ).toInt >= 30 && cline(anPos(damID)).split(":").apply(GQ).toInt >= 30 && cline(anPos(proID)).split(":").apply(GQ).toInt >= 30) &&
(cline(anPos(sireID)).split(":").apply(DP).toInt >= 16 && cline(anPos(damID)).split(":").apply(DP).toInt >= 16 && cline(anPos(proID)).split(":").apply(DP).toInt >= 16)){

if (cline(anPos(sireID)).apply(1) == '|' && cline(anPos(damID)).apply(1) == '|' && cline(anPos(proID)).apply(1) == '|'  && (
	cline(anPos(sireID)).apply(0) != '.' && cline(anPos(damID)).apply(0) != '.' && cline(anPos(proID)).apply(0) != '.')){
if (haps(sireID + "L").size == 10){

haps(sireID + "L") =  cline(anPos(sireID)).apply(0) :: haps(sireID + "L").dropRight(1)
haps(sireID + "R") =  cline(anPos(sireID)).apply(2) :: haps(sireID + "R").dropRight(1)

haps(damID + "L") =  cline(anPos(damID)).apply(0) :: haps(damID + "L").dropRight(1)
haps(damID + "R") =  cline(anPos(damID)).apply(2) :: haps(damID + "R").dropRight(1)

haps(proID + "L") =  cline(anPos(proID)).apply(0) :: haps(proID + "L").dropRight(1)
haps(proID + "R") =  cline(anPos(proID)).apply(2) :: haps(proID + "R").dropRight(1)

for (chl <- chlID){
haps(chl + "L") =  cline(anPos(chl)).apply(0) :: haps(chl + "L").dropRight(1)
haps(chl + "R") =  cline(anPos(chl)).apply(2) :: haps(chl + "R").dropRight(1)
}
	
if (!(((haps(proID + "L").sameElements(haps(sireID + "L")) || haps(proID + "L").sameElements(haps(sireID + "R"))) || 
	(haps(proID + "R").sameElements(haps(sireID + "L")) || haps(proID + "R").sameElements(haps(sireID + "R")))) &&
	((haps(proID + "L").sameElements(haps(damID + "L")) || haps(proID + "L").sameElements(haps(damID + "R"))) || 
	(haps(proID + "R").sameElements(haps(damID + "L")) || haps(proID + "R").sameElements(haps(damID + "R")))))){
	for (chl <- chlID){
	if((haps(proID + "L").sameElements(haps(chl + "L")) || haps(proID + "L").sameElements(haps(chl + "R"))) || 
	(haps(proID + "R").sameElements(haps(chl + "L")) || haps(proID + "R").sameElements(haps(chl + "R")))) inher = true
	}
	//if (count <= 0 && inher){
	if (inher){
	if (!(haps(proID + "L").sameElements(haps(sireID + "L")) || haps(proID + "L").sameElements(haps(sireID + "R")) || 
	(haps(proID + "R").sameElements(haps(sireID + "L")) || haps(proID + "R").sameElements(haps(sireID + "R"))))){
	println("Fail Sire")
	}
	if (!((haps(proID + "L").sameElements(haps(damID + "L")) || haps(proID + "L").sameElements(haps(damID + "R"))) || 
	(haps(proID + "R").sameElements(haps(damID + "L")) || haps(proID + "R").sameElements(haps(damID + "R"))))){
	println("Fail Dam")
	}
	//if (count <= 0){
	println(cline(0) + ":" + cline(1))
	println("New Haplotype " + haps(proID + "L") + " " + haps(proID + "R"))
	println("Dam " + haps(damID + "L") + " " + haps(damID + "R"))
	println("Sire " + haps(sireID + "L") + " " + haps(sireID + "L"))
	for (chl <- chlID){
	println(chl + " " + haps(chl + "L") + " " + haps(chl + "R"))
	}
	println("--")
	count = 10
	inher = false
	} else {
	count -= 1
	}
	}

} else { //eif have haps

haps(sireID + "L") =  cline(anPos(sireID)).apply(0) :: haps(sireID + "L")
haps(sireID + "R") =  cline(anPos(sireID)).apply(2) :: haps(sireID + "R")

haps(damID + "L") =  cline(anPos(damID)).apply(0) :: haps(damID + "L")
haps(damID + "R") =  cline(anPos(damID)).apply(2) :: haps(damID + "R")

haps(proID + "L") =  cline(anPos(proID)).apply(0) :: haps(proID + "L")
haps(proID + "R") =  cline(anPos(proID)).apply(2) :: haps(proID + "R")

for (chl <- chlID){
haps(chl + "L") =  cline(anPos(chl)).apply(0) :: haps(chl + "L")
haps(chl + "R") =  cline(anPos(chl)).apply(2) :: haps(chl + "R")
}

} //eelse nohaps

} //eif Phased
}//pass filters
}//ewhile



//x.sameElements(y)

} // main


def setFields(geno: String): Unit = {
	DP = geno.split(":").indexOf("DP")
	GQ = geno.split(":").indexOf("GQ")
}

def shareHap(id1: String, id2: String): Boolean = {
if (haps(id1 + "L").sameElements(haps(id2 + "L")) || haps(id1 + "L").sameElements(haps(id2 + "R")) || haps(id1 + "R").sameElements(haps(id2 + "L")) || haps(id1 + "R").sameElements(haps(id2 + "R"))) true
else false
}

} //object