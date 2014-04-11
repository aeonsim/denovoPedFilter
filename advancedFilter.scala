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


	def itPed (pop: HashMap[String, Array[String]] , rec: String): List[String] = {
		if (pop.contains(rec)){
			return rec :: itPed(pop,pop(rec).apply(2)) ::: itPed(pop,pop(rec).apply(3))
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
	
	def rtGTs(genotype: String) : List[String] = {
		val gts = if (genotype.size == 3 && genotype.contains('/')) genotype.split('/') else genotype.split('|')
		gts.toList
	}

/* Main body of Program */

def main (args: Array[String]): Unit = {

	var settings = new HashMap[String,String]
	
	for (items <- args){
		val keyVal = items.split("=")
		settings += keyVal(0).toUpperCase -> keyVal(1) 
	}
	

	if ((! settings.contains("VCF")) && (! settings.contains("PED")) & (! settings.contains("TRIOS"))) {
		println("advFilter VCF=input.vcf.gz PED=input.ped TRIOS=input_probands.txt { minDP=0 minALT=0 RECUR=F/T minKIDS=1 PLGL=0,0,0 QUAL=0 }")
		println("Trios = txtfile per line: AnimalID\tavgDepth")
		println("{} Optional arguments, Values shown are default")
		System.exit(1)
	}

	val in_vcf = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream(settings("VCF")))))
	val in_ped = new BufferedReader(new FileReader(settings("PED")))
	val in_pro = new BufferedReader(new FileReader(settings("TRIOS")))
	
	
	/* Default settings*/
	val minDP = if (settings.contains("MINDP")) settings("MINDP").toInt else 0
	val minALT = if (settings.contains("MINALT")) settings("MINALT").toInt else 0
	val reoccur = if (settings.contains("RECUR") && List("TRUE","YES","Y","T").contains(settings("RECUR"))) true else false
	val QUAL = if (settings.contains("QUAL")) settings("QUAL").toInt else 0
	val PLGLS = if (settings.contains("PLGL")) settings("PLGL").split(",") else Array(0,0,0)
	val minKids = if (settings.contains("MINKIDS")) settings("MINKIDS").toInt else 1

	val outname = args(0).split("/")(args(0).split("/").size - 1)
	val out_vcf = new BufferedWriter(new OutputStreamWriter(new BlockCompressedOutputStream(outname + ".mutations-" + reoccur + "-denovos.vcf.gz")))
	val out_somatic = new BufferedWriter(new OutputStreamWriter(new BlockCompressedOutputStream(outname + ".mutations-" + reoccur + "-somatic.vcf.gz")))
	
	var pedFile = new HashMap[String, Array[String]]

	// Tuple6(Ancestors, Parents, Children, Descendants, TrioDP, population, extended-Fam)

	var trios = new HashMap[String, Tuple7[List[String],List[String],List[String],List[String],Int,List[String],List[String]]]
	var vcfanimals = new HashMap[String, Int]
	var ancestors : List[String] = Nil
	var parents : List[String] = Nil
	var children : List[String] = Nil
	var descendents : List[String] = Nil
	var extFam : List[String] = Nil
	var population : List[String] = Nil
	var AD = -1
	var GT = -1
	var DP = -1
	var RO = -1
	var AO = -1
	var PL = -1
	var PLexist = false
	

/*
* Iterate through VCF file to find Column headers and store positions in Map for later recall
*/

	var vcfrec = in_vcf.readLine().split("\t")
	while (vcfrec(0).apply(1) == '#'){
		out_vcf.write(vcfrec.reduceLeft{(a,b) => a + "\t" + b} + "\n")
		out_somatic.write(vcfrec.reduceLeft{(a,b) => a + "\t" + b} + "\n")
		vcfrec = in_vcf.readLine().split("\t")
	}

	out_vcf.write(s"##${args.reduceLeft{(a,b) => a + " " + b}}\n")
	out_somatic.write(s"##${args.reduceLeft{(a,b) => a + " " + b}}\n")

	if (vcfrec(0).apply(0) == '#' && vcfrec(0).apply(1) == 'C'){
		out_vcf.write(vcfrec.reduceLeft{(a,b) => a + "\t" + b} + "\n")
		out_somatic.write(vcfrec.reduceLeft{(a,b) => a + "\t" + b} + "\n")
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
		descendents = Nil
		extFam = Nil
		if (pedFile.contains(curPro(0)) && animalIDS.contains(curPro(0))){
			parents = pedFile(curPro(0))(2) :: pedFile(curPro(0))(3) :: Nil
			ancestors = itPed(pedFile,curPro(0)).tail.filterNot(x => parents.contains(x) || (! vcfanimals.contains(x)))
			children = findChildren(pedFile,vcfanimals,curPro(0))
			for (indv <- animalIDS.toList){ //.filterNot(x => ( (x == curPro(0)) || children.contains(x) || ancestors.contains(x) || parents.contains(x)))){
				val tempPed = itPed(pedFile,indv)
				if (tempPed.contains(curPro(0))){
					descendents = indv :: descendents
				}
				if (pedFile.contains(indv) && ( pedFile(indv)(2) == parents(0) || pedFile(indv)(3) == parents(0) )  && pedFile(indv)(1) != curPro(0)){
					extFam = indv :: extFam
				}
			}//efor
			var tmpdesc = descendents.filterNot(x => ((x == curPro(0)) || children.contains(x) || ancestors.contains(x) || parents.contains(x)|| (! vcfanimals.contains(x))))
			
			population = animalIDS.toList.filterNot(x => ( (x == curPro(0)) || children.contains(x) || ancestors.contains(x) || descendents.contains(x) || parents.contains(x)|| (! vcfanimals.contains(x))))
			
			if (animalIDS.contains(parents(0)) && animalIDS.contains(parents(1)) && animalIDS.contains(curPro(0)) && children.size != 0){
				trios += curPro(0) -> (ancestors, parents, children, tmpdesc, curPro(1).toInt, population, extFam)
			}
		}//eif
	}//ewhile

	in_pro.close

println("Built Pedigrees")
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
		print("EXFAM\t")
		fam._2._7.foreach(s => print(s + "\t"))
		print(fam._2._7.size + "\n\n")
	}

/*
*	Iterate through VCF file line by line, at each line load each Trio and count existence of variants in different categories
*	if de novo, flag and output snp detail and variant info, (count in pop, children ancestors etc)
*/
	println(s"Chrom\tPos\tRef\tRefSize\tAlt\tAltSize\tQUAL\tTrio\tGenotype\tAnces\tPars\tChildren\tDesc\tExFam\tPop\tPopFreq\tSupport Ratio\tScore\tClass\tProband\tSire\tDam\tPopRefCount\tPopAltCount\tWarning")

	while (in_vcf.ready){
		var denovo = false
		var line = in_vcf.readLine().split("\t")
		val format = line(8).split(":")
		if (format.contains("PL") || format.contains("GL")){
			PL = if (format.contains("PL")) format.indexOf("PL") else format.indexOf("GL")
		}
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
		println(line(5).toFloat + " " + QUAL + " " + (line(5).toFloat >= QUAL))
		if (line.size == (vcfanimals.size + 9) && (line(5).toFloat >= QUAL)){

			for (fam <- trios){
				var altsPar = 0
				val ped = fam._2
				var ances, par, kids, desc, popFreq, exFamFreq = 0
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
					var grands: List[String] = Nil
					for (indv <- ped._1){
						if (line(vcfanimals(indv))(0) != '.'){
							val curAn = line(vcfanimals(indv)).split(":")
							val refAlt = selROvAD(curAn,AD, RO, AO, GT)
							adratio += rc(refAlt._1,refAlt._2)						
							grands = grands ++ rtGTs(curAn(GT))
							if (sigAD(refAlt._2)){
								ances += 1
							}
						}
					}
					
					//println(grands + " " + proBand(GT))
					
					if (proBand(GT)(0) != '.' && grands.contains(proBand(GT)(0)) && grands.contains(proBand(GT)(1))){
						ances += 1
					}

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
			var popRef, popALT = 0

					for (indv <- ped._6){
						if (line(vcfanimals(indv))(0) != '.'){
							val curAn = line(vcfanimals(indv)).split(":")
							val refAlt = selROvAD(curAn,AD, RO, AO, GT)
							if ((isVar(curAn(GT)) || sigAD(refAlt._2)) && checkDP(curAn,DP,minDP,maxDP)){
								popRef += refAlt._1
								popALT += refAlt._2
								val curIndv = pedFile(indv)
								if (vcfanimals.contains(curIndv(2)) && vcfanimals.contains(curIndv(3))){
								val sire = line(vcfanimals(curIndv(2))).split(":") 
								val dam = line(vcfanimals(curIndv(3))).split(":")
								if (sire(0).apply(0) != '.' && dam(0).apply(0) != '.'){
									if (isVar(sire(GT)) || isVar(dam(GT))) {
										popFreq += 1
									}
								} else {
									popFreq += 1
								}
							} else {
									popFreq += 1
								}
							}	
						}
					}
					
			/* Extended family Calc */

					for (indv <- ped._7){
						if (line(vcfanimals(indv))(0) != '.'){
							val curAn = line(vcfanimals(indv)).split(":")
							val refAlt = selROvAD(curAn,AD, RO, AO, GT)
							if ((isVar(curAn(GT)) || sigAD(refAlt._2)) && checkDP(curAn,DP,minDP,maxDP)){
								exFamFreq += 1
							}
						}
					}
					//errors.println(s"${line(0)}\t${line(1)}\tAnces\t${ances}\tPar\t${par}\tkids\t${kids}\tdesc\t${desc}\t\t${popFreq - (1 + kids + desc)}")

			/* Check Pedigree segregation pattern */

			/* -------- Denovo Check ---------- */

					val curPro = line(vcfanimals(fam._1)).split(":")
					var rank = 0
					altsPar match {
						case 0 => rank = 9
						case 1 => rank = 8
						case _ => rank = 7
					}
					val proRatio = selROvAD(proBand,AD, RO, AO, GT)
					val proGT = proBand(GT)
					
					if (((!valGTs.contains(proBand(GT)(0).toString + proBand(GT)(2)) && proRatio._2 >= minALT) || 
							(proRatio._2 >= (minALT * 3))) 
								&& checkDP(curPro, DP, minDP, maxDP) &&  checkDP(par1,DP,minDP,maxDP) && checkDP(par2,DP,minDP,maxDP) &&
							(ances == 0) && (par == 0) && (kids >= minKids) &&  (if (reoccur){true}else{popFreq == 0})
					){
						denovo = true
						
						if ((proRatio._1 + proRatio._2) <= minDP) {
							errors.println(s"minDP == ${minDP}\t${proRatio._1} + ${proRatio._2}\t ${proBand(DP)}")
							proBand.foreach(er => errors.print(er))
							errors.print("\n")
							}
						
						if (reoccur && adratio == 0.0){
							println(s"${line(0)}\t${line(1)}\t${line(3)}\t${line(3).size}\t${line(4)}\t${line(4).size}\t${line(5)}\t${fam._1}\t${proGT}\t${ances}\t${par}\t${kids}\t${desc}\t${exFamFreq}\t${popFreq}\t${popFreq.toFloat/(animalIDS.size)}\t${proRatio._2/proRatio._1.toFloat}\t${rank}\tdenovo\t" + proRatio + "\t" + selROvAD(par1,AD, RO, AO, GT) + "\t" + selROvAD(par2,AD, RO, AO, GT) + s"\t${popRef}\t${popALT}\t")
							out_vcf.write(line.reduceLeft{(a,b) => a + "\t" + b} + "\n")
						}else {
							print(s"${line(0)}\t${line(1)}\t${line(3)}\t${line(3).size}\t${line(4)}\t${line(4).size}\t${line(5)}\t${fam._1}\t${proGT}\t${ances}\t${par}\t${kids}\t${desc}\t${exFamFreq}\t${popFreq}\t${popFreq.toFloat/(animalIDS.size)}\t${proRatio._2/proRatio._1.toFloat}\t${rank}\tdenovo\t" + proRatio + "\t" + selROvAD(par1,AD, RO, AO, GT) + "\t" + selROvAD(par2,AD, RO, AO, GT) + s"\t${popRef}\t${popALT}")
							if (adratio != 0.0) {
								line(6) = "LOWQUAL_ADratio"
								print ("\t WARNING: Low confidence de novo\n")
							}else {
								print("\n")
							}//eelse
							out_vcf.write(line.reduceLeft{(a,b) => a + "\t" + b} + "\n")
						}//eif reoccur

					} //eif is Denovo
					
					
					if(!valGTs.contains(proBand(GT)(0).toString + proBand(GT)(2)) && (ances == 0) && (par == 0) && (kids == 0) 
						&& (popFreq == 0) && checkDP(curPro, DP, minDP, maxDP) && checkDP(par1,DP,minDP,maxDP) && checkDP(par2,DP,minDP,maxDP)
						&& proRatio._2 >= minALT && adratio == 0.0){
							println(s"${line(0)}\t${line(1)}\t${line(3)}\t${line(3).size}\t${line(4)}\t${line(4).size}\t${line(5)}\t${fam._1}\t${proGT}\t${ances}\t${par}\t${kids}\t${desc}\t${exFamFreq}\t${popFreq}\t${popFreq.toFloat/(animalIDS.size)}\t${proRatio._2/proRatio._1.toFloat}\t${rank}\tSomatic\t" + proRatio + "\t" + selROvAD(par1,AD, RO, AO, GT) + "\t" + selROvAD(par2,AD, RO, AO, GT) + s"\t${popRef}\t${popALT}\t")
							out_somatic.write(line.reduceLeft{(a,b) => a + "\t" + b} + "\n")
						}
					
					

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