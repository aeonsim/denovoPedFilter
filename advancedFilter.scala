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

/* Global Var for VCF Type*/
var vcfType = ""
val errors = System.err

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
				if (altGT == 0){
					alt = alts(0).toInt
				} else {
					alt = alts(altGT - 1).toInt
				}
				if (refGT == 0){
					ref = if (vcfType == "platypus") indv(ROval).split(",")(0).toInt else indv(ROval).toInt
				} else {
					ref = if (vcfType == "platypus") (indv(ROval).split(",")(refGT - 1).toInt - alt) else alts(refGT - 1).toInt
				}
			}//eif
		}//eelse
			(ref,alt)
	}

/* Take a Genotype String and check against DP limits*/

	def checkDP (genos: Array[String], DPpos: Int, minDP: Int, maxDP: Int): Boolean = {
		if (genos.size >= 2 && DPpos != -1){
			val curDP = if (vcfType == "platypus") genos(DPpos).split(",")(0).toInt else genos(DPpos).toInt
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

	var trios = new HashMap[String, Tuple6[List[String],List[String],List[String],List[String],Int,List[String]]]
	var vcfanimals = new HashMap[String, Int]
	var ancestors : List[String] = Nil
	var parents : List[String] = Nil
	var children : List[String] = Nil
	var descendents : List[String] = Nil
	var population : List[String] = Nil
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
			parents = pedFile(curPro(0))(2) :: pedFile(curPro(0))(3) :: Nil
			ancestors = itPed(pedFile,vcfanimals,curPro(0)).tail.filterNot(x => parents.contains(x))
			children = findChildren(pedFile,vcfanimals,curPro(0))
			for (kid <- children){
				descendents = findChildren(pedFile,vcfanimals,kid) 
			}//efor
			population = animalIDS.toList.filterNot(x => (children.contains(x) || ancestors.contains(x) || descendents.contains(x) || parents.contains(x)))
			if (animalIDS.contains(parents(0)) && animalIDS.contains(parents(1)) && animalIDS.contains(curPro(0)) && children.size != 0){
			trios += curPro(0) -> (ancestors, parents, children, descendents, curPro(1).toInt, population)
			}
		}//eif
	}//ewhile

	in_pro.close

//Closed proband file

/* Report identified Trios & there Pedigree Structure*/

	for (fam <- trios){
		println(s"TRIO:\t${fam._1}\t${fam._2._2(0)}\t${fam._2._2(1)}")
		print("Grandparents\t")
		fam._2._1.foreach(s => print(s + "\t"))
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
	println(s"Chrom\tPos\tRef\tRefSize\tAlt\tAltSize\tQUAL\tTrio\tGenotype\tAnces\tPars\tChildren\tDesc\tPop\tPopFreq\tSupport Ratio\tScore\tProband\tSire\tDam\tWarning")

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
			DP = format.indexOf("NR")
			vcfType = "platypus"
		}

		if (line.size == (vcfanimals.size + 9)){

			for (fam <- trios){
				var altsPar = 0
				val ped = fam._2
				var ances, par, kids, desc, popFreq = 0
				val maxDP = (ped._5 * 1.7).toInt
				var adratio = 0.0
				
			/* Parental Test using permutations of Alleles */

				val par1 = line(vcfanimals(ped._2.apply(0))).split(":")
				val par2 = line(vcfanimals(ped._2.apply(1))).split(":")
				val proBand = line(vcfanimals(fam._1)).split(":")

				if (par1(GT)(0) != '.' && par2(GT)(0) != '.' && proBand(GT)(0) != '.'){
					val valGTs = permu(par1(GT)(0).toString + par1(GT)(2),par2(GT)(0).toString + par2(GT)(2))
					if (valGTs.contains(proBand(GT)(0).toString + proBand(GT)(2))){
						par += 1
					}

			/* After checking that the Parent GT's can produce the Denovo check AD/RO/AO incase GT misscall */

					for (indv <- ped._2){
						if (line(vcfanimals(indv))(0) != '.'){
							val curAn = line(vcfanimals(indv)).split(":")
							val refAlt = selROvAD(curAn,AD, RO, AO, GT)
							altsPar += refAlt._2
							adratio += rc(refAlt._1,refAlt._2)
							if (sigAD(refAlt._2)){
								par += 1
							}
						}
					}


			/* Loop through each family group and record Hets */
			//Ansc
					for (indv <- ped._1){
						if (line(vcfanimals(indv))(0) != '.'){
							val curAn = line(vcfanimals(indv)).split(":")
							val refAlt = selROvAD(curAn,AD, RO, AO, GT)
							adratio += rc(refAlt._1,refAlt._2)
							if (isVar(curAn(GT)) || sigAD(refAlt._2)){
								ances += 1
							}
						}
					}//Efor ansc

			//Children
					for (indv <- ped._3){
						if (line(vcfanimals(indv))(0) != '.'){
					//print(line(vcfanimals(indv)).split(":").apply(GT) + " " + isVar(line(vcfanimals(indv)).split(":").apply(GT)) + "\t")
							val curAn = line(vcfanimals(indv)).split(":")
							val refAlt = selROvAD(curAn,AD, RO, AO, GT)
							if (isVar(curAn(GT)) || sigAD(refAlt._2)){
								kids += 1
							}
						}
					}

			//Desec
					for (indv <- ped._4){
						if (line(vcfanimals(indv))(0) != '.'){
					//print(line(vcfanimals(indv)).split(":").apply(GT) + " " + isVar(line(vcfanimals(indv)).split(":").apply(GT)) + "\t")
							val curAn = line(vcfanimals(indv)).split(":")
							val refAlt = selROvAD(curAn,AD, RO, AO, GT)
							if (isVar(curAn(GT)) || sigAD(refAlt._2)){
								desc += 1
							}
						}
					}

			/* Population Calc */

					for (indv <- ped._6){
						if (line(vcfanimals(indv))(0) != '.'){
							val curAn = line(vcfanimals(indv)).split(":")
							val refAlt = selROvAD(curAn,AD, RO, AO, GT)
							if ((isVar(curAn(GT)) || sigAD(refAlt._2)) && checkDP(curAn,DP,minDP,maxDP)){
								popFreq += 1
							}
						}
					}
					
					errors.println(s"${line(0)}\t${line(1)}\tAnces\t${ances}\tPar\t${par}\tkids\t${kids}\tdesc\t${desc}\t\t${popFreq - (1 + kids + desc)}")

			/* Check Pedigree segregation pattern */

			/* -------- Denovo Check ---------- */

					val curPro = line(vcfanimals(fam._1)).split(":")
	
					if (((!valGTs.contains(proBand(GT)(0).toString + proBand(GT)(2)) && selROvAD(proBand,AD, RO, AO, GT)._2 >= minALT) || 
							(selROvAD(proBand,AD, RO, AO, GT)._2 >= (minALT * 3))
						) 	&& checkDP(curPro, DP, minDP, maxDP) &&  checkDP(par1,DP,minDP,maxDP) && checkDP(par2,DP,minDP,maxDP) &&
							(ances == 0) && (par == 0) && (kids >= 1) &&  (if (reoccur){true}else{popFreq > 0})
					){
						denovo = true
						
						var rank = 0
						altsPar match {
							case 0 => rank = 9
							case 1 => rank = 8
							case _ => rank = 7
						}

						val proGT = proBand(GT)
						val proRatio = selROvAD(proBand,AD, RO, AO, GT)
						
						if ((proRatio._1 + proRatio._2) <= minDP) {
							errors.println(s"minDP == ${minDP}\t${proRatio._1} + ${proRatio._2}\t ${proBand(DP)}")
							proBand.foreach(er => errors.print(er))
							errors.print("\n")
							}
						
						if (reoccur && adratio == 0.0){
							println(s"${line(0)}\t${line(1)}\t${line(3)}\t${line(3).size}\t${line(4)}\t${line(4).size}\t${line(5)}\t${fam._1}\t${proGT}\t${ances}\t${par}\t${kids}\t${desc}\t${popFreq}\t${popFreq.toFloat/(animalIDS.size)}\t${proRatio._2/proRatio._1.toFloat}\t${rank}\t" + proRatio + "\t" + selROvAD(par1,AD, RO, AO, GT) + "\t" + selROvAD(par2,AD, RO, AO, GT) + "\t")
							out_vcf.write(line.reduceLeft{(a,b) => a + "\t" + b} + "\n")
						}else {
							print(s"${line(0)}\t${line(1)}\t${line(3)}\t${line(3).size}\t${line(4)}\t${line(4).size}\t${line(5)}\t${fam._1}\t${proGT}\t${ances}\t${par}\t${kids}\t${desc}\t${popFreq}\t${popFreq.toFloat/(animalIDS.size)}\t${proRatio._2/proRatio._1.toFloat}\t${rank}\t" + proRatio + "\t" + selROvAD(par1,AD, RO, AO, GT) + "\t" + selROvAD(par2,AD, RO, AO, GT))
							if (adratio != 0.0) {
								line(6) = "LOWQUAL_ADratio"
								print ("\t WARNING: Low confidence de novo\n")
							}else {
								print("\n")
							}//eelse
							out_vcf.write(line.reduceLeft{(a,b) => a + "\t" + b} + "\n")
						}//eif reoccur

					} //eif is Denovo

				}//eisVAR
			}//Efor fam <- trios
		} else {
			println(s"Error ${line(0)} ${line(1)} ${line.size}")
	} //else

	}// Ewhile

	in_vcf.close
	out_vcf.close
}//eMain

}//eObject