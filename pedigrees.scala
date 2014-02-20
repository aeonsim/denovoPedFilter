/* Advanced Pedigree Filter
* Takes in VCF (gz) & Pedigree file, with cutoffs for DP, and a list of Proband ID's
* Use pedigree file to build family tree for each Proband, Proband, Children, Descendants, Parents, Ancestors
* Select variants based on pedigree & min DP check AD for low coverage
* Parents & Ancestors -ve, Proband +ve, Children 1+ positive, Descendants possibly +ve
* Output VCF per Proband ID of possible De novos along with, score file with counts in the population.
* Check to see if each VCF line is the correct length & discard if not.
*/

object pedFind{

import java.io._
import scala.collection.mutable.HashMap

/* Global Var for VCF Type*/
var vcfType = ""


/* Recursive function to find all parents & ancestors of nominated individual
* returns list of those present in both the PED & VCF files.
*/



	def itPed (pop: HashMap[String, Array[String]], rec: String): List[String] = {
		if (pop.contains(rec)){
			return rec :: itPed(pop,pop(rec).apply(2)) ::: itPed(pop,pop(rec).apply(3))
		} else {
			return Nil
		}
	}

/* Finds Children for a nominated parent based on Pedigree records
* Only finds ones that are both in the VCF & Ped file
*/

	def findChildren (pop: HashMap[String, Array[String]] , parent: String) : List[String]={
		var result: List[String] = Nil
		for (item <- pop){
			if ((item._2(2) == parent)||(item._2(3) == parent) ){
				result = item._1 :: result
			}//eif

		}//efor
		result
	}//edef

/* Main body of Program */

def main (args: Array[String]): Unit = {
println("PED File, List of Trios")

	val in_ped = new BufferedReader(new FileReader(args(0)))
	val in_pro = new BufferedReader(new FileReader(args(1)))

	var pedFile = new HashMap[String, Array[String]]

	// Tuple5(Ancestors, Parents, Children, Descendants, TrioDP)

	var trios = new HashMap[String, Tuple5[List[String],List[String],List[String],List[String],Int]]
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
			ancestors = itPed(pedFile,curPro(0)).tail
			parents = pedFile(curPro(0))(2) :: pedFile(curPro(0))(3) :: Nil
			children = findChildren(pedFile,curPro(0))
			for (kid <- children){
				descendents = findChildren(pedFile,kid) 
			}//efor
			trios += curPro(0) -> (ancestors, parents, children, descendents, curPro(1).toInt)

		}//eif
	}//ewhile

	in_pro.close

//Closed proband file

/* Report identified Trios & there Pedigree Structure*/
	println(s"Proband\tSire\tDam\nGrandparents\nChildren\nDescendents")

	for (fam <- trios){
		println(s"TRIO:\t${fam._1}\t${fam._2._2(0)}\t${fam._2._2(1)}")
		fam._2._1.foreach(s => print(s + "\t"))
		print("Grandparents\t")
		print(fam._2._1.size + "\n")
		print("Children\t")		
		fam._2._3.foreach(s => print(s + "\t"))
		print(fam._2._3.size + "\n")
		print("Descendents\t")
		fam._2._4.foreach(s => print(s + "\t"))
		print(fam._2._4.size + "\n")
	}

/*
*	Iterate through VCF file line by line, at each line load each Trio and count existence of variants in different categories
*	if de novo, flag and output snp detail and variant info, (count in pop, children ancestors etc)
*/
// Ewhile

}//eMain

}//eObject