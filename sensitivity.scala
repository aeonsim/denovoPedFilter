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

	def itPed (pop: HashMap[String, Array[String]], vcf: HashMap[String,Int] , rec: String, lim: Int): List[String] = {
		if (lim < 2){
		if (pop.contains(rec) & vcf.contains(rec)){
			return rec :: itPed(pop,vcf,pop(rec).apply(2),lim + 1) ::: itPed(pop,vcf,pop(rec).apply(3), lim + 1)
		} else {
			return Nil
		}
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
		List("1/0","0/1","1|0","0|1").contains(genotype)
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
		var alt = -1
		var ref = -1
		val gtsplit = if (indv(GTval).contains('/')) indv(GTval).split('/') else indv(GTval).split('|')
		var refGT = gtsplit(0).toInt
		var altGT = gtsplit(1).toInt // initial work to deal with multiple alleles ie GT = 2/3
		if (ADval != -1){
			ref = indv(ADval).split(",")(refGT).toInt
			if (altGT == 0) {
				alt = indv(ADval).split(",")(altGT + 1).toInt
				}
			else {
				alt = indv(ADval).split(",")(altGT).toInt
			}
		} else {
			val alts = indv(AOval).split(",")
			if (ROval != -1 && AOval != -1){
				if (refGT == 0){
					ref = indv(ROval).toInt
				} else {
					ref = alts(refGT - 1).toInt
				}
				if (altGT == 0){
					alt = alts(0).toInt
				} else {
					alt = alts(altGT - 1).toInt
				}
			}//eif
		}//eelse
			(ref,alt)
	}

/* Take a Genotype String and check against DP limits*/

	def checkDP (genos: Array[String], DPpos: Int, minDP: Int, maxDP: Int): Boolean = {
		if (genos.size >= 2){
		 if (DPpos != -1){
			val curDP = genos(DPpos).toInt
			if (curDP >= minDP && curDP <= maxDP) true
			else false
		} else {
			false
		}
		} else{
			false
		}
	}

/* Main body of Program */

def main (args: Array[String]): Unit = {

	if (args.size <= 2) {
		println("advFilter input.vcf.gz input.ped input_probands.txt minDP minALT reoccurring?\n")
		System.exit(1)
	}

	println(s"Chrom\tPos\tRef\tRefSize\tAlt\tAltSize\tQUAL\tTrio\tGenotype\tAnces\tPars\tChildren\tDesc\tPop\tPopFreq\tSupport Ratio\tProband\tSire\tDam\tWarning")
	val in_vcf = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream(args(0)))))
	val in_ped = new BufferedReader(new FileReader(args(1)))
	val in_pro = new BufferedReader(new FileReader(args(2)))
	val minDP = args(3).toInt
	val minALT = args(4).toInt
	val reoccur = if(List("TRUE","YES","Y","T").contains(args(5).toUpperCase)) true else false

	val outname = args(0).split("/")(args(0).split("/").size - 1)
	val out_vcf = new BufferedWriter(new OutputStreamWriter(new BlockCompressedOutputStream(outname + ".reoccur-" + reoccur + "-denovos.vcf.gz")))

	var pedFile = new HashMap[String, Array[String]]

	// Tuple5(Ancestors, Parents, Children, Descendants, TrioDP)

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

	val animalIDS = vcfanimals.keys.toArray

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

		if (pedFile.contains(curPro(0)) && animalIDS.contains(curPro(0))){
			ancestors = itPed(pedFile,vcfanimals,curPro(0), 0).tail
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

	var badHets, goodHets = 0

	while (in_vcf.ready){
		var denovo = false
		var line = in_vcf.readLine().split("\t")
		val format = line(8).split(":")
		AD = format.indexOf("AD")
		GT = format.indexOf("GT")
		DP = format.indexOf("DP")
		AO = format.indexOf("AO")
		RO = format.indexOf("RO")
		
        if (format.contains("NV")){
        	AO = format.indexOf("NV")
            RO = format.indexOf("NR")
        }
        
        if (DP == -1 && format.contains("NR")){
        	DP = format.indexOf("NR")
        }
        
                
		if (line.size == (vcfanimals.size + 9)){

			for (fam <- trios){
				val ped = fam._2
				var ances, par = 0
					
				val curPro = line(vcfanimals(fam._1)).split(":")
				
	/*	Check State in Grandparents*/

					for (indv <- ped._1){
						if (line(vcfanimals(indv)).size > 3){
							val curAn = line(vcfanimals(indv)).split(":")
							if (isHet(curAn(GT))){
								ances += 1
							}
						}
					}				

			/* After checking that the Parent GT's can produce the Denovo check AD/RO/AO incase GT misscall */

					for (indv <- ped._2){
						if (line(vcfanimals(indv)).size > 3){
							val curAn = line(vcfanimals(indv)).split(":")
							if (isHet(curAn(GT))){
								par += 1
							}
						}
					}

			if ( isHet(curPro(GT)) && ances == 1 && par == 1){
				goodHets += 1
			}
			if (isHet(curPro(GT)) && ances == 1 && par == 0){
				badHets += 1
			}

			}//Efor fam <- trios
		} else {
			println(s"Error ${line(0)} ${line(1)} ${line.size}")
	} //else

	}// Ewhile
	
	println("Total Good Hets\t= " + goodHets)
	println("Total Bad Hets\t= " + badHets)
	println("Total Tested Hets\t= " + (goodHets + badHets))
	
}//eMain

}//eObject