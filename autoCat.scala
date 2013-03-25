import java.io._
import scala.collection.mutable.HashMap

object autoCat{


def main(args: Array[String]): Unit = {

//Note currently GATK file has BUG with AD values Doubled, so have doubled Cat8 Check fix when new GATK arrives.
var bugfixMultiplier = 1.0

if (args.size < 4){
println("autoCat input.vcf sireID damID proID optional{AvgDPpro minDP minQualityScore}")
println("With the optional filters autoCat will discard Variants 1.8x AvgDPpro & > minQualityScore")
System.exit(1)
}

val sire = args(1) //"21930910"
val dam = args(2) //"22328328"
val proband = args(3) //"23934418"
val inName = args(0) //"test.vcf"
val filterDP = if (args.size >= 5) args(4).toInt * 1.8 else 100000
val minDP = if (args.size >= 5) args(5).toInt  else 0
val filterQual = if (args.size >= 6) args(6).toFloat else 0.0

val pileupFile = new File(inName + ".pileup")
var pileup = new BufferedReader(new FileReader(inName))

if (pileupFile.exists()){
println("Pileup supplied assuming Sire, Dam, Proband Order, using Pileup analysis")
pileup = new BufferedReader(new FileReader(pileupFile))
} else {
println("Using VCF AD")
}

println(s"Filters are:\nDP <= ${filterDP}\nQUAL >\t${filterQual}")

val input = new BufferedReader(new FileReader(inName))

val cat9 = new BufferedWriter(new FileWriter(inName + ".cat9.vcf"))
val cat8 = new BufferedWriter(new FileWriter(inName + ".cat8.vcf"))
val cat7 = new BufferedWriter(new FileWriter(inName + ".cat7.vcf"))
val catX = new BufferedWriter(new FileWriter(inName + ".catX.vcf"))

var head = input.readLine

while (head(0) == '#' & head(1) != 'C'){
cat9.write(head + "\n")
cat8.write(head + "\n")
cat7.write(head + "\n")
catX.write(head + "\n")
head = input.readLine
}

var cline = head.split("\t")

val sirePos = cline.indexOf(sire)
val damPos = cline.indexOf(dam)
val proPos = cline.indexOf(proband)

var cat9score, cat8score, cat7score, catXscore, filterScore = 0

if (sirePos == -1 | damPos == -1 | proPos == -1){
System.err.println("ERROR: VCF file is missing a Trio member (Position of: Sire, Dam, Proband): " + sirePos + " " + damPos + " " + proPos)
System.exit(1)
}

val header =("##INFO=<ID=CAT9,Number=0,Type=Flag,Description=\"High Confidence De novo\">\n" + 
			"##INFO=<ID=CAT8,Number=0,Type=Flag,Description=\"Moderate Confidence De novo\">\n" +
		"##INFO=<ID=CAT7,Number=0,Type=Flag,Description=\"Low Confidence De novo\">\n" +
		"##INFO=<ID=CATX,Number=0,Type=Flag,Description=\"Probably Not De novo\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sire + "\t" + dam + "\t" + proband + "\n")
cat9.write(header)
cat8.write(header)
cat7.write(header)
catX.write(header)

var pileupData = new HashMap[String,Array[String]]

while (pileup.ready){
var pline = pileup.readLine.split("\t")
if ( pline.size == 12) pileupData += s"${pline(0)}:${pline(1)}" -> pline
}

var DP,AD,PL = 0


while (input.ready){
cline = input.readLine.split("\t")

DP = cline(8).split(":").indexOf("DP")
AD = cline(8).split(":").indexOf("AD")
PL = cline(8).split(":").indexOf("PL")

//if (DP == -1 | AD == -1 | PL == -1){
if ( AD == -1 | PL == -1){
System.err.println("ERROR: VCF file is missing Annotation (Position of:  AD, PL): " + AD + " " + PL + " @")
println(cline(0) + ":" + cline(1))
//System.exit(1)
}

if (AD != -1 && cline(damPos).size > 3 && cline(sirePos).size > 3 && cline(proPos).size > 3  &&  cline(5).toFloat >= filterQual && getDP(cline(proPos)) <= filterDP
	&& getDP(cline(proPos)) >= minDP && ( pileupData.contains(s"${cline(0)}:${cline(1)}") || (validAD(cline(proPos)) && validAD(cline(damPos)) &&  validAD(cline(sirePos))))){

if (isDenovo(cline(damPos),cline(sirePos),cline(proPos))){
cat9.write(s"${cline(0)}\t${cline(1)}\t${cline(2)}\t${cline(3)}\t${cline(4)}\t${cline(5)}\tPASS\tCAT9;${cline(7)}\t${cline(8)}\t${cline(sirePos)}\t${cline(damPos)}\t${cline(proPos)}\n")
cat9score += 1
} else { //Eif isDenovo
	if (isCat8(cline(damPos),cline(sirePos),cline(proPos))){
		cat8.write(s"${cline(0)}\t${cline(1)}\t${cline(2)}\t${cline(3)}\t${cline(4)}\t${cline(5)}\tPASS\tCAT8;${cline(7)}\t${cline(8)}\t${cline(sirePos)}\t${cline(damPos)}\t${cline(proPos)}\n")
		cat8score += 1
	} else { //
		if (isCat7(cline(damPos),cline(sirePos),cline(proPos))){
			cat7.write(s"${cline(0)}\t${cline(1)}\t${cline(2)}\t${cline(3)}\t${cline(4)}\t${cline(5)}\tPASS\tCAT7;${cline(7)}\t${cline(8)}\t${cline(sirePos)}\t${cline(damPos)}\t${cline(proPos)}\n")
			cat7score += 1
		} else {
			catX.write(s"${cline(0)}\t${cline(1)}\t${cline(2)}\t${cline(3)}\t${cline(4)}\t${cline(5)}\tFAIL\tCATX;${cline(7)}\t${cline(8)}\t${cline(sirePos)}\t${cline(damPos)}\t${cline(proPos)}\n")
			catXscore += 1
		}
	}
}//isDenovo Else
} else {
filterScore += 1
}
}//Ewhile


def isDenovo(dam: String, sire: String, pro:String) : Boolean = {
//Assumes pileup order Sire, Dam , Proband
var damAD, sireAD, proAD = 0
if (pileupData.contains(s"${cline(0)}:${cline(1)}")){
damAD = pileup2AD(pileupData(s"${cline(0)}:${cline(1)}").apply(7).toUpperCase).split(",").apply(1).toInt
sireAD = pileup2AD(pileupData(s"${cline(0)}:${cline(1)}").apply(4).toUpperCase).split(",").apply(1).toInt
proAD = pileup2AD(pileupData(s"${cline(0)}:${cline(1)}").apply(10).toUpperCase).split(",").apply(1).toInt
} else {
damAD = dam.split(":").apply(AD).split(",").apply(1).toInt
sireAD = sire.split(":").apply(AD).split(",").apply(1).toInt
proAD = pro.split(":").apply(AD).split(",").apply(1).toInt
}
if (( damAD == 0) &
		( sireAD == 0) &
			( proAD != 0)
	) true else false
}

def isCat8(dam: String, sire: String, pro:String) : Boolean = {
var damAD, sireAD, proAD = 0
if (pileupData.contains(s"${cline(0)}:${cline(1)}")){
damAD = pileup2AD(pileupData(s"${cline(0)}:${cline(1)}").apply(7).toUpperCase).split(",").apply(1).toInt
sireAD = pileup2AD(pileupData(s"${cline(0)}:${cline(1)}").apply(4).toUpperCase).split(",").apply(1).toInt
proAD = pileup2AD(pileupData(s"${cline(0)}:${cline(1)}").apply(10).toUpperCase).split(",").apply(1).toInt
} else {
damAD = dam.split(":").apply(AD).split(",").apply(1).toInt
sireAD = sire.split(":").apply(AD).split(",").apply(1).toInt
proAD = pro.split(":").apply(AD).split(",").apply(1).toInt
}
if ((( damAD <= (1 * bugfixMultiplier)) ^ //XOR one may be true but not both
		( sireAD <= (1 * bugfixMultiplier))) &
			( proAD != 0)
	) true else false

}

def isCat7(dam: String, sire: String, pro:String) : Boolean = {
var damAD, sireAD, proAD = Array("","")
if (pileupData.contains(s"${cline(0)}:${cline(1)}")){
damAD = pileup2AD(pileupData(s"${cline(0)}:${cline(1)}").apply(7).toUpperCase).split(",")
sireAD = pileup2AD(pileupData(s"${cline(0)}:${cline(1)}").apply(4).toUpperCase).split(",")
proAD = pileup2AD(pileupData(s"${cline(0)}:${cline(1)}").apply(10).toUpperCase).split(",")
} else {
damAD = dam.split(":").apply(AD).split(",")
sireAD = sire.split(":").apply(AD).split(",")
proAD = pro.split(":").apply(AD).split(",")
}
if ((( (damAD(1).toFloat/damAD(0).toFloat) < 0.1) ^ //XOR one may be true but not both
		( (sireAD(1).toFloat/sireAD(0).toFloat) < 0.1)) &
			( proAD(1) != 0)
	) true else false

}

def getDP (geno: String): Int ={
if (DP == -1){
val ADdata = geno.split(":").apply(AD).split(",")
if (ADdata.size < 2) 0 else ADdata(0).toInt + ADdata(1).toInt
} else {
geno.split(":").apply(DP).toInt
}

}

def validAD (geno: String): Boolean = {
val ADdata = geno.split(":").apply(AD).split(",")
if (ADdata.size == 2) true else false
}

def pileup2AD (bases: String): String = {
val ref = bases.count(_==',') + bases.count(_=='.')
val alt = bases.count(_=='A') + bases.count(_=='T') + bases.count(_=='C') + bases.count(_=='G') + bases.count(_=='-') + bases.count(_=='+') 
s"$ref,$alt"
//$alt
}

val total : Float = cat9score + cat8score + cat7score + catXscore + filterScore
println(s"Cat9\t$cat9score\t${cat9score/total}")
println(s"Cat8\t$cat8score\t${cat8score/total}")
println(s"Cat7\t$cat7score\t${cat7score/total}")
println(s"CatX\t$catXscore\t${catXscore/total}")
println(s"Filtered\t$filterScore\t${filterScore/total}")
println("Total\t" + total)


input.close
cat9.close
cat8.close
cat7.close
catX.close

}//main
}//class