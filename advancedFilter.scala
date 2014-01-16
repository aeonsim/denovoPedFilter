/* Advanced Pedigree Filter
* Takes in VCF (gz) & Pedigree file, with cutoffs for DP, and a list of Proband ID's
* Use pedigree file to build family tree for each Proband, Proband, Children, Descendants, Parents, Ancestors
* Select variants based on pedigree & min DP check AD for low coverage
* Parents & Ancestors -ve, Proband +ve, Children 1+ positive, Descendants possibly +ve
* Output VCF per Proband ID of possible De novos along with, score file with counts in the population.
* Check to see if each VCF line is the correct length & discard if not.
*/

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
if (vcf.contains(item) && ( (item._2(2) == parent)||(item._2(3) == parent) )){
result = parent :: result
}//eif

}//efor
result
}//edef


def main (args: Array[String]): Unit = {

if (args.size <= 2) println("advFilter input.vcf.gz input.ped input_probands.txt minDP")

val in_vcf = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream(args(0)))))
val in_ped = new BufferedReader(new FileReader(args(1)))
val in_pro = new BufferedReader(new FileReader(args(2)))
val minDP = args(3).toInt

var pedFile = new HashMap[String, Array[String]]

// Tuple4(Ancestors, Parents, Children, Descendants)

var trios = new HashMap[String, Tuple4[List[String],List[String],List[String],List[String]]]
var vcfanimals = new HashMap[String, Int]
var ancestors : List[String] = Nil
var parents : List[String] = Nil
var children : List[String] = Nil
var descendents : List[String] = Nil

/*
* Iterate through VCF file to find Column headers and store positions in Map for later recall
*/


var vcfrec = in_vcf.readLine().split("\t")
while (vcfrec(0).apply(1) == '#'){
vcfrec = in_vcf.readLine().split("\t")
}

if (vcfrec(0).apply(0) == '#' && vcfrec(0).apply(1) == 'C'){
for (i <- 0 to (vcfrec.size -1)){
vcfanimals += vcfrec(i) -> i
}
}

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
ancestors = itPed(pedFile,vcfanimals,curPro)
parents = pedFile(curPro)(2) :: pedFile(curPro)(3) :: Nil
children = findChildren(pedFile,vcfanimals,curPro)

for (kid <- children){
descendents = findChildren(kid) 
}//efor

trios += curPro -> (ancestors, parents, children, descendents)

}//eif


}//ewhile

in_pro.close

//Closed




}//eMain

