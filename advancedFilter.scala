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

def isHet(genotype:String): Boolean ={
if (genotype.size == 3) {
if ( genotype(0) != genotype(2)) true
else false
} else false
}


def main (args: Array[String]): Unit = {

if (args.size <= 2) {
println("advFilter input.vcf.gz input.ped input_probands.txt minDP")
System.exit(1)
}
println(s"Log file: ${args(2)}.log\n")


val in_vcf = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream(args(0)))))
val in_ped = new BufferedReader(new FileReader(args(1)))
val in_pro = new BufferedReader(new FileReader(args(2)))
val log = new BufferedWriter(new FileWriter(args(2)+".log"))
val minDP = args(3).toInt

var pedFile = new HashMap[String, Array[String]]

// Tuple4(Ancestors, Parents, Children, Descendants)

var trios = new HashMap[String, Tuple4[List[String],List[String],List[String],List[String]]]
var vcfanimals = new HashMap[String, Int]
var ancestors : List[String] = Nil
var parents : List[String] = Nil
var children : List[String] = Nil
var descendents : List[String] = Nil
var AD = -1
var GT = -1
var DP = -1

/*
* Iterate through VCF file to find Column headers and store positions in Map for later recall
*/

//println("Start VCF ID extraction")

var vcfrec = in_vcf.readLine().split("\t")
while (vcfrec(0).apply(1) == '#'){
vcfrec = in_vcf.readLine().split("\t")
}

if (vcfrec(0).apply(0) == '#' && vcfrec(0).apply(1) == 'C'){
for (i <- 0 to (vcfrec.size -1)){
vcfanimals += vcfrec(i) -> i
}
}

val animalIDS = vcfanimals.keys.toList.slice(9,vcfanimals.size)
println(animalIDS + " Got Slice")
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
val curPro = in_pro.readLine()

if (pedFile.contains(curPro)){
ancestors = itPed(pedFile,vcfanimals,curPro).tail
parents = pedFile(curPro)(2) :: pedFile(curPro)(3) :: Nil
children = findChildren(pedFile,vcfanimals,curPro)

for (kid <- children){
descendents = findChildren(pedFile,vcfanimals,kid) 
}//efor

trios += curPro -> (ancestors, parents, children, descendents)

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
val line = in_vcf.readLine().split("\t")

val format = line(8).split(":")
AD = format.indexOf("AD")
GT = format.indexOf("GT")
DP = format.indexOf("DP")

//println(line.size + " " + vcfanimals.size)
if (line.size == vcfanimals.size){

for (fam <- trios){
//println(fam + " in fam loop")
val ped = fam._2
var ances, par, kids, desc, popFreq = 0

/*
* Loop through each family group and record Hets
*/


for (indv <- ped._1){
if (isHet(line(vcfanimals(indv)).split(":").apply(GT))){
ances += 1
}
}//Efor ansc

for (indv <- ped._2){
if (isHet(line(vcfanimals(indv)).split(":").apply(GT))){
par += 1
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
if (isHet(line(vcfanimals(indv)).split(":").apply(GT))){
popFreq += 1
}
}

/*
* Check Pedigree segregation pattern
*/
log.write(s"${line(0)}\t${line(1)}\t${ances}\t${par}\t${kids}\t${desc}\n")
log.write(line(vcfanimals(fam._1)).split(":").apply(GT) + "\t" + isHet(line(vcfanimals(fam._1)).split(":").apply(GT)) + "\n")

if (isHet(line(vcfanimals(fam._1)).split(":").apply(GT)) && ances == 0 && par == 0 && kids >= 1 && ((kids + desc + 1) == popFreq)){
denovo = true
println(s"Denovo\t${line(0)}\t${line(1)}\tTrio:\t${fam._1}\t${line(vcfanimals(fam._1)).split(":").apply(GT)}\t${ances}\t${par}\t${kids}\t${desc}\n")
}

}//Efor fam <- trios

} else {
println(s"Error ${line(0)} ${line(1)} ${line.size}")
} //else

}// Ewhile

log.close
in_vcf.close
}//eMain

}//eObject