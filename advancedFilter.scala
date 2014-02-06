/* Advanced Pedigree Filter
* Takes in VCF (gz) & Pedigree file, with cutoffs for DP, and a list of Proband ID's
* Use pedigree file to build family tree for each Proband, Proband, Children, Descendants, Parents, Ancestors
* Select variants based on pedigree & min DP check AD for low coverage
* Parents & Ancestors -ve, Proband +ve, Children 1+ positive, Descendants possibly +ve
* Output VCF per Proband ID of possible De novos along with, score file with counts in the population.
* Check to see if each VCF line is the correct length & discard if not.
*/

object advFilter{

import java.io._
import net.sf.samtools.util._
import scala.collection.mutable.HashMap

/* Recursive function to find all parents & ancestors of nominated individual
* returns list of those present in both the PED & VCF files.
*/

def itPed (pop: HashMap[String, Array[String]], vcf: HashMap[String,Int] , rec: String): List[String] = {
if (pop.contains(rec) & vcf.contains(rec)){
return rec :: itPed(pop,vcf,pop(rec).apply(2)) ::: itPed(pop,vcf,pop(rec).apply(3))
} else {
return Nil
}
}

/* Takes Parental Genotypes & Generates all combinations*/

def permu(str: String, str1: String): List[String] = {
val result: List[String] = for (i <- str.toList; j <- str1.toList) yield i.toString.concat(j.toString)
result
}

/* Finds Children for a nominated parent based on Pedigree records
* Only finds ones that are both in the VCF & Ped file
*/

def findChildren (pop: HashMap[String, Array[String]], vcf: HashMap[String,Int] , parent: String) : List[String]={
var result: List[String] = Nil
for (item <- pop){
if (vcf.contains(item._1) && ( (item._2(2) == parent)||(item._2(3) == parent) )){
result = item._1 :: result
}//eif

}//efor
result
}//edef

/* Is this a Het variant? */

def isHet(genotype:String): Boolean = {
if (genotype.size == 3) {
if ( genotype(0) != genotype(2)) true
else false
} else false
}

/* Calc Ratio of Ref/Alt*/

def rc(ref: Int, alt: Int): Float = {
alt/ref.toFloat
}

/* Are there two or more Alt reads? */
def sigAD(alt: Int): Boolean = {
if (alt >= 2 ) true
else false
}

/* Is this a Variant ie not ./. or 0/0 */

def isVar(genotype:String): Boolean = {
if (genotype.size == 3 && (genotype(0) != '.')) {
if ((genotype(0) != '0' && genotype(2) != '0')) {
true
} else {
 false
}
} else {
 false
 }
}

/* Extract the correct Read counts from AD or RO/AO AD format is REF,ALT,ALT2,... */

def selROvAD(indv: Array[String], ADval: Int, ROval: Int, AOval: Int, GTval: Int): Tuple2[Int,Int] = {
//println(s"${indv.toList} ${ADval} ${ROval} ${AOval}")
var alt = -1
var ref = -1
val gtsplit = if (indv(GTval).contains('/')) indv(GTval).split('/') else indv(GTval).split('|')
var refGT = gtsplit(0).toInt
var altGT = gtsplit(1).toInt // initial work to deal with multiple alleles ie GT = 2/3

if (ADval != -1){
ref = indv(ADval).split(",")(refGT).toInt
alt = indv(ADval).split(",")(altGT).toInt
} else {
//println(refGT + " " + altGT)
val alts = indv(AOval).split(",")

if (ROval != -1 && AOval != -1){

if (refGT == 0){
ref = indv(ROval).toInt
//print(s"ref==0 ${ref}\t")
} else {
ref = alts(refGT - 1).toInt
//print(s"ref!=0 ${ref}\t")
}
if (altGT == 0){
alt = alts(0).toInt
//print(s"alt==0 ${alt}\t")
} else {
alt = alts(altGT - 1).toInt
//print(s"alt!=0 ${alt}\t")
}

//alt = indv(AOval).split(",")(0).toInt
}//eif
//if (indv(AOval).split(",").size >= 2 && ref == 0){
//ref = indv(AOval).split(",")(0).toInt
//alt = indv(AOval).split(",")(1).toInt
//} else {
//if (indv(AOval).split(",").size >= 3 && alt == 0){
//ref = indv(AOval).split(",")(1).toInt
//alt = indv(AOval).split(",")(2).toInt
//}
//}
}//eelse
(ref,alt)
}

/* Take a Genotype String and check against DP limits*/

def checkDP (genos: Array[String], DPpos: Int, minDP: Int, maxDP: Int): Boolean = {
if (genos.size >= 2 && DPpos != -1){
val curDP = genos(DPpos).toInt
if (curDP >= minDP && curDP <= maxDP) true
else false
} else {
false
}
}

/* Main body of Program */

def main (args: Array[String]): Unit = {

if (args.size <= 2) {
println("advFilter input.vcf.gz input.ped input_probands.txt minDP minALT reoccurring?\n")
System.exit(1)
}

println(s"Chrom\tPos\tRef\tRefSize\tAlt\tAltSize\tQUAL\tTrio\tGenotype\tAnces\tPars\tChildren\tDesc\tPop\tPopFreq\tProband\tSire\tDam\tWarning")
val in_vcf = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream(args(0)))))
val in_ped = new BufferedReader(new FileReader(args(1)))
val in_pro = new BufferedReader(new FileReader(args(2)))
//val log = new BufferedWriter(new FileWriter(args(2)+".log"))
val minDP = args(3).toInt
val minALT = args(4).toInt
val reoccur = if(List("TRUE","YES","Y","T").contains(args(5).toUpperCase)) true else false

val outname = args(0).split("/")(args(0).split("/").size - 1)

val out_vcf = new BufferedWriter(new FileWriter(outname + ".reoccur-" + reoccur + "-denovos.vcf"))

var pedFile = new HashMap[String, Array[String]]

// Tuple4(Ancestors, Parents, Children, Descendants, TrioDP)

var trios = new HashMap[String, Tuple5[List[String],List[String],List[String],List[String],Int]]
var vcfanimals = new HashMap[String, Int]
var ancestors : List[String] = Nil
var parents : List[String] = Nil
var children : List[String] = Nil
var descendents : List[String] = Nil
var AD = -1
var GT = -1
var DP = -1
var RO = -1
var AO = -1

/*
* Iterate through VCF file to find Column headers and store positions in Map for later recall
*/

var vcfrec = in_vcf.readLine().split("\t")
while (vcfrec(0).apply(1) == '#'){
out_vcf.write(vcfrec.reduceLeft{(a,b) => a + "\t" + b} + "\n")
vcfrec = in_vcf.readLine().split("\t")
}

if (vcfrec(0).apply(0) == '#' && vcfrec(0).apply(1) == 'C'){
out_vcf.write(vcfrec.reduceLeft{(a,b) => a + "\t" + b} + "\n")
for (i <- 9 to (vcfrec.size -1)){
vcfanimals += vcfrec(i) -> i
}
}

val animalIDS = vcfanimals.keys

/*
* Build Map file of Pedigree
*/

while (in_ped.ready){
val temp = in_ped.readLine().split("\t")
pedFile += temp(1) -> temp
}

in_ped.close

/*
* Iterate through our list of probands finding their sequenced Ancestors, Parents, Children
* & other descendants, to use for filtering, and build the Trios Map
*/

while (in_pro.ready){
val curPro = in_pro.readLine.split("\t")

if (pedFile.contains(curPro(0))){
ancestors = itPed(pedFile,vcfanimals,curPro(0)).tail
parents = pedFile(curPro(0))(2) :: pedFile(curPro(0))(3) :: Nil
children = findChildren(pedFile,vcfanimals,curPro(0))

for (kid <- children){
descendents = findChildren(pedFile,vcfanimals,kid) 
}//efor

trios += curPro(0) -> (ancestors, parents, children, descendents, curPro(1).toInt)

}//eif

}//ewhile

in_pro.close

//Closed proband file

/*
*	Iterate through VCF file line by line, at each line load each Trio and count existence of variants in different categories
*	if de novo, flag and output snp detail and variant info, (count in pop, children ancestors etc)
*/

while (in_vcf.ready){
var denovo = false
var line = in_vcf.readLine().split("\t")
val format = line(8).split(":")
AD = format.indexOf("AD")
GT = format.indexOf("GT")
DP = format.indexOf("DP")
AO = format.indexOf("AO")
RO = format.indexOf("RO")

//println(line.size + " " + vcfanimals.size)
if (line.size == (vcfanimals.size + 9)){

for (fam <- trios){
//println(fam + " in fam loop")
val ped = fam._2
var ances, par, kids, desc, popFreq = 0
val maxDP = (ped._5 * 1.7).toInt
/*
* Loop through each family group and record Hets
*/


for (indv <- ped._1){
if (line(vcfanimals(indv)).size > 3){
val curAn = line(vcfanimals(indv)).split(":")
val refAlt = selROvAD(curAn,AD, RO, AO, GT)
if (isVar(curAn(GT)) || sigAD(refAlt._2)){
ances += 1
}
}
}//Efor ansc

/* Parental Test using permutations of Alleles */

val par1 = line(vcfanimals(ped._2.apply(0))).split(":")
val par2 = line(vcfanimals(ped._2.apply(1))).split(":")
val proBand = line(vcfanimals(fam._1)).split(":")

if (par1(GT).size == 3 && par2(GT).size == 3 && proBand(GT).size == 3){
val valGTs = permu(par1(GT)(0).toString + par1(GT)(2),par2(GT)(0).toString + par2(GT)(2))
if (valGTs.contains(proBand(GT)(0).toString + proBand(GT)(2))){
par += 1
}
}


/* After checking that the Parent GT's can produce the Denovo check AD/RO/AO incase GT misscall */

for (indv <- ped._2){
if (line(vcfanimals(indv)).size > 3){
val curAn = line(vcfanimals(indv)).split(":")
val refAlt = selROvAD(curAn,AD, RO, AO, GT)
if (sigAD(refAlt._2)){
par += 1
}
}
}

for (indv <- ped._3){
if (isHet(line(vcfanimals(indv)).split(":").apply(GT))){
kids += 1
}
}

for (indv <- ped._4){
if (isHet(line(vcfanimals(indv)).split(":").apply(GT))){
desc += 1
}
}

for (indv <- animalIDS){
//print(line(vcfanimals(indv)).split(":").apply(GT) + " ")
if (line(vcfanimals(indv)).size > 3){
val curAn = line(vcfanimals(indv)).split(":")
val refAlt = selROvAD(curAn,AD, RO, AO, GT)
if ((isVar(curAn(GT)) || sigAD(refAlt._2)) && (if(reoccur){curAn(DP).toInt >= minDP}else{true}) ){
popFreq += 1
}
}
}

/*
* Check Pedigree segregation pattern
*/

/* -------- Denovo Check ---------- */

val curPro = line(vcfanimals(fam._1)).split(":")

if (((isHet(curPro(GT)) && selROvAD(proBand,AD, RO, AO, GT)._2 >= minALT) || ( isVar(curPro(GT)) && selROvAD(proBand,AD, RO, AO, GT)._2 >= (minALT * 3))) 
	&& checkDP(curPro, DP, minDP, maxDP) && 
	 checkDP(par1,DP,minDP,maxDP) && checkDP(par2,DP,minDP,maxDP) &&
	(ances == 0) && (par == 0) && (kids >= 1) && 
	 (if (reoccur){true}else{(kids + desc + 1) == popFreq})
	){
denovo = true
var adratio = 0.0
/* Check RefAlt Ratio in case of mis call */

for (indv <- ped._1){
if(line(vcfanimals(indv)).size > 3){
val curAn = line(vcfanimals(indv)).split(":")
val refAlt = selROvAD(curAn,AD, RO, AO, GT)
adratio += rc(refAlt._1,refAlt._2)
}
}//Efor ansc

for (indv <- ped._2){
if(line(vcfanimals(indv)).size > 3){
val curAn = line(vcfanimals(indv)).split(":")
val refAlt = selROvAD(curAn,AD, RO, AO, GT)
adratio += rc(refAlt._1,refAlt._2)
}
}

if (reoccur && adratio == 0.0){
println(s"${line(0)}\t${line(1)}\t${line(3)}\t${line(3).size}\t${line(4)}\t${line(4).size}\t${line(5)}\t${fam._1}\t${line(vcfanimals(fam._1)).split(":").apply(GT)}\t${ances}\t${par}\t${kids}\t${desc}\t${popFreq}\t${popFreq.toFloat/(animalIDS.size)}\t" + selROvAD(proBand,AD, RO, AO, GT) + "\t" + selROvAD(par1,AD, RO, AO, GT) + "\t" + selROvAD(par2,AD, RO, AO, GT) + "\t" + proBand(RO) + " " + proBand(AO))
out_vcf.write(line.reduceLeft{(a,b) => a + "\t" + b} + "\n")
}else {
print(s"${line(0)}\t${line(1)}\t${line(3)}\t${line(3).size}\t${line(4)}\t${line(4).size}\t${line(5)}\t${fam._1}\t${line(vcfanimals(fam._1)).split(":").apply(GT)}\t${ances}\t${par}\t${kids}\t${desc}\t${popFreq}\t${popFreq.toFloat/(animalIDS.size)}\t" + selROvAD(proBand,AD, RO, AO, GT) + "\t" + selROvAD(par1,AD, RO, AO, GT) + "\t" + selROvAD(par2,AD, RO, AO, GT))
if (adratio != 0.0) {
line(6) = "LOWQUAL_ADratio"
print ("\t WARNING: Low confidence de novo\n")
}else {
print("\n")
}//eelse
out_vcf.write(line.reduceLeft{(a,b) => a + "\t" + b} + "\n")
}//eif reoccur


}
//log.write(s"${line(0)}\t${line(1)}\t${line(2)}\t${line(3)}\t${line(5)}\t${ances}\t${par}\t${kids}\t${desc}\tPop: ${popFreq} (${popFreq.toFloat/(animalIDS.size)})\t${line(vcfanimals(fam._1)).split(':').apply(GT)}\n")
}//Efor fam <- trios
} else {
println(s"Error ${line(0)} ${line(1)} ${line.size}")
} //else

}// Ewhile

//log.close
in_vcf.close
out_vcf.close
}//eMain

}//eObject