import java.io._
import net.sf.samtools.util._
import scala.collection.mutable.HashMap

object genie{

var GQ, DP = -1
var haps = new HashMap[String,List[Char]]

def main(args: Array[String]) : Unit = {

var GQmin = 20
var DPmin = 8
var numSNPS = 3

if (args.size < 5) {
println("genie VCF.gz sire dam proband children.list\n Gzipped VCF file Sire, Dam & Proband IDs, List of Children 1 per line\n")
System.exit(1)
}

val vcf = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream(args(0))))) ///args(0)))))
val children = new BufferedReader(new FileReader(args(4)))

val sireID = args(1) //"21930910" //
val damID = args(2) //"22328328" //
val proID = args(3) //"23934418" //
var anPos = new HashMap[String,Int]
var chlID: List[String] = Nil

//Load Children IDs
while (children.ready) chlID = children.readLine :: chlID

var cline = vcf.readLine.split("\t")

//Skip headers
while (cline(0).apply(1) == '#') cline = vcf.readLine.split("\t")

//Load VCF Animal ID's and positions
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

var inher,posEvent  = false
var fails,numInher = 0
var pfail = ""

while(vcf.ready){
cline = vcf.readLine.split("\t")

//skip multiAllelic sites (for now)
setFields(cline(8))
 
//System.err.println(s"Ref & Alt ${cline(3)} ${cline(4)} size: ${cline(4).split(",").size}")
 
if(cline(4).split(",").size > 1){
System.err.println(s"${cline(0)} ${cline(1)} ${cline(2)} ${cline(3)} ${cline(4)}")
cline = vcf.readLine.split("\t")
setFields(cline(8))
}

//filter GQ & DP
if ( (GQ == -1 && DP == -1) || (cline(anPos(sireID)).split(":").apply(GQ).toInt >= GQmin && 
		cline(anPos(damID)).split(":").apply(GQ).toInt >= GQmin && 
			cline(anPos(proID)).split(":").apply(GQ).toInt >= GQmin) &&
				(cline(anPos(sireID)).split(":").apply(DP).toInt >= DPmin && 
					cline(anPos(damID)).split(":").apply(DP).toInt >= DPmin && 
						cline(anPos(proID)).split(":").apply(DP).toInt >= DPmin)){

if (cline(anPos(sireID)).apply(1) == '|' &&
	 	cline(anPos(damID)).apply(1) == '|' && 
	 		cline(anPos(proID)).apply(1) == '|'  &&
	 			(cline(anPos(sireID)).apply(0) != '.' &&
	 				cline(anPos(damID)).apply(0) != '.' && 
	 					cline(anPos(proID)).apply(0) != '.')){

if (haps(sireID + "L").size == numSNPS){

haps(sireID + "L") =  cline(anPos(sireID)).apply(0) :: haps(sireID + "L").dropRight(1)
haps(sireID + "R") =  cline(anPos(sireID)).apply(2) :: haps(sireID + "R").dropRight(1)

haps(damID + "L") =  cline(anPos(damID)).apply(0) :: haps(damID + "L").dropRight(1)
haps(damID + "R") =  cline(anPos(damID)).apply(2) :: haps(damID + "R").dropRight(1)

haps(proID + "L") =  cline(anPos(proID)).apply(0) :: haps(proID + "L").dropRight(1)
haps(proID + "R") =  cline(anPos(proID)).apply(2) :: haps(proID + "R").dropRight(1)

for (chl <- chlID){
if (cline(anPos(chl)).apply(1) == '|'){
haps(chl + "L") =  cline(anPos(chl)).apply(0) :: haps(chl + "L").dropRight(1)
haps(chl + "R") =  cline(anPos(chl)).apply(2) :: haps(chl + "R").dropRight(1)
} else {
haps(chl + "L") =  '.' :: haps(chl + "L").dropRight(1)
haps(chl + "R") =  '.' :: haps(chl + "R").dropRight(1)
}
}
//}

//Checks to see if either parent lacks a shared Haplotype

if (!(shareHap(proID,sireID))){
pfail = sireID
posEvent = true
} else {
if (!(shareHap(proID,damID))){
pfail = damID
posEvent = true
}
}

if (!(shareHap(proID,sireID)) && !(shareHap(proID,damID))){
 posEvent = false
pfail = ""
}

//If a parent failed then check the children	
if (posEvent){
	//Is inherited?
	for (chl <- chlID){
	//if (shareHapSire(proID,chl) && !(shareHapSire(damID,chl)) && !(shareHapSire(sireID,chl))) {
	if (shareHap(proID,chl) && !(shareHap(damID,chl)) && !(shareHap(sireID,chl))) {
	inher = true
	numInher += 1
	}
	}

	if (inher && count <= 0){
	println("Failed: " + pfail  + "    numChildren " + numInher)
	println(cline(0) + ":" + cline(1))
	println("New Hapl " + haps(proID + "L") + " " + haps(proID + "R") + " " + snpHap(proID))
	println("Dam      " + haps(damID + "L") + " " + haps(damID + "R") + " " + snpHap(damID))
	println("Sire     " + haps(sireID + "L") + " " + haps(sireID + "R") + " " + snpHap(sireID))
	for (chl <- chlID){
	println(chl + " " + haps(chl + "L") + " " + haps(chl + "R"))
	}
	println("--")
	count = numSNPS
	inher = false
	posEvent = false
	numInher = 0
	pfail = ""
	} else {
	count -= 1
	numInher = 0
	inher = false
	posEvent = false
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
if (cline(anPos(chl)).apply(1) == '|'){
haps(chl + "L") =  cline(anPos(chl)).apply(0) :: haps(chl + "L")
haps(chl + "R") =  cline(anPos(chl)).apply(2) :: haps(chl + "R")
} else {
haps(chl + "L") =  '.' :: haps(chl + "L")
haps(chl + "R") =  '.' :: haps(chl + "R")
}
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

def shareHapSire(id1: String, id2: String): Boolean = {
if (haps(id1 + "L").sameElements(haps(id2 + "L")) || haps(id1 + "L").sameElements(haps(id2 + "R"))) true
else false
}

def snpHap(id: String): List[Int] ={
(haps(id + "L"),haps(id + "R")).zipped.map(_.toString.toInt + _.toString.toInt)
}

} //object