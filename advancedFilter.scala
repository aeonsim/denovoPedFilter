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
var PLexist = false
var PL = -1

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
			if ((genotype(0) == '0' && genotype(2) == '0')) {
				false
			} else {
 				true
			}
		} else {
 			false
 		}
	}
	
	/* Return Genotypes as List*/
	
	def rtGTs(genotype: String) : List[String] = {
		val gts = if (genotype.size == 3 && genotype.contains('/')) genotype.split('/') else genotype.split('|')
		gts.toList
	}
	
	/* Select method for returning Ref & Alt allele counts*/

	def selROvAD(indv: Array[String], ADval: Int, ROval: Int, AOval: Int, GTval: Int): Tuple2[Int,Int] = {
		var refAlt = (1,-1)
		vcfType match {
			case "gatk" => refAlt = gatkAD(indv,ADval,GTval)
			case "platypus" => refAlt = platypusRefAlt(indv,ROval, AOval, GTval)
			case "freebayes" => refAlt = fbRefAlt(indv,ROval, AOval, GTval)
			case _ => println("\n\n\nUnknown VCF type in SelROvAD"); System.exit(1)
		}
		refAlt
	}

	/* GATK return Ref & Alt counts*/
	def gatkAD(indv: Array[String], ADval: Int, GTval: Int): Tuple2[Int,Int] = {
		if (ADval != -1 && indv.size > ADval && indv(ADval) != "."){
			val RefAlt = indv(ADval).split(",")
			val GT = rtGTs(indv(GTval)).sorted
			(RefAlt(GT(0).toInt).toInt,RefAlt(GT(1).toInt).toInt)
		} else {
		if (PLexist){
			//Check other Genotypes are atleast Phred 10 unlikely
			if (indv(PL).split(",").sorted.tail.forall(_.toInt > 10)) (1,0) else (1,-1)
		}else {
			(1,-1)
			}
		}
	}
	
	/* Platy return Ref & Alt counts*/
	def platypusRefAlt(indv: Array[String], ROval: Int, AOval: Int, GTval: Int): Tuple2[Int,Int] = {
		var alt = -1
		var ref = 1
		
		if (indv.size > AOval){
			val GT = rtGTs(indv(GTval)).sorted
			var refGT = GT(0).toInt
			var altGT = GT(1).toInt
			val alts = indv(AOval).split(",")
			if (ROval != -1 && AOval != -1){
				if (altGT == 0){
					alt = alts(0).toInt
				} else {
					alt = alts(altGT - 1).toInt
				}
				if (refGT == 0){
					ref =  indv(ROval).split(",")(0).toInt - alt
				} else {
					ref = (indv(ROval).split(",")(refGT - 1).toInt - alt) 
				}
			}
		}
		(ref,alt)
	}

	/* Freebayes return Ref & Alt counts*/
	def fbRefAlt(indv: Array[String], ROval: Int, AOval: Int, GTval: Int): Tuple2[Int,Int] = {
		var alt = -1
		var ref = 1
		
		if (indv.size > AOval){
			val GT = rtGTs(indv(GTval)).sorted
			var refGT = GT(0).toInt
			var altGT = GT(1).toInt
			val alts = indv(AOval).split(",")
			if (ROval != -1 && AOval != -1){
				if (altGT == 0){
					alt = alts(0).toInt
				} else {
					alt = alts(altGT - 1).toInt
				}
				if (refGT == 0){
					ref =  indv(ROval).toInt
				} else {
					ref = alts(refGT - 1).toInt
				}
			}

		} 
		(ref,alt)
	}



/* Take a Genotype String and check against DP limits*/

	def checkDP (genos: Array[String], DPpos: Int, minDP: Int, maxDP: Int): Boolean = {
		if (genos.size >= 2 && DPpos != -1 && genos(DPpos) != "."){
			val curDP = if (vcfType == "platypus") genos(DPpos).split(",")(0).toInt else genos(DPpos).toInt
			if (curDP >= minDP && curDP <= maxDP) true
			else false
		} else {
			false
		}
	}
	
	/* Check PL is sufficient, if no PL's then aways ok*/
	
	def checkPL(min: Int, indv: Array[String]): Boolean = {
		if (PLexist){
			val pls = indv(PL).split(",").map(num => scala.math.abs(num.toFloat)).sorted.tail
			if (pls(0) >= min) true 
			else false
		} else {
			true
		}	
	}
	
	
/* Phase Code, return format is (sireAllele,damAllele)*/

def phase(indv: Array[String], sire: Array[String], dam: Array[String]) : Tuple2[String, String] = {
if (sire.size >= PL && dam.size >= PL && indv.size >= PL){
	val minPLval = if (vcfType == "gatk") 40 else 5
	val indvGT = indv(0)
	val sireGT = sire(0)
	val damGT = dam(0)
	val indvPL = indv(PL).split(",").map(num => scala.math.abs(num.toFloat)).sorted.tail
	val sirePL = sire(PL).split(",").map(num => scala.math.abs(num.toFloat)).sorted.tail
	val damPL = dam(PL).split(",").map(num => scala.math.abs(num.toFloat)).sorted.tail
if (indvPL(0) >= minPLval && sirePL(0) >= minPLval && damPL(0) >= minPLval){
	if(sireGT == "0/0" && damGT == "1/1" && (indvGT == "0/1" || indvGT == "1/0")){
		("0","1")
	} else {
	 if (sireGT == "1/1" && damGT == "0/0" && (indvGT == "0/1" || indvGT == "1/0")){
	 	("1","0")
	 } else {
		if(sireGT == "0/1" && (indvGT == "0/1" || indvGT == "1/0") && damGT == "0/0") {
			("1","0")
		} else {
				if(sireGT == "0/1" && (indvGT == "0/1" || indvGT == "1/0") && damGT == "1/1") {
					("0","1")
				} else {
						if (damGT == "0/1" && (indvGT == "0/1" || indvGT == "1/0") && sireGT == "0/0" ) {
							("0","1")
						} else {
							if (damGT == "0/1" && (indvGT == "0/1" || indvGT == "1/0") && sireGT == "1/1") {
								("1","0")
							} else {
							("x","x")
							}
					}
			}
		}		
	}
	
	}
	} else ("x","x")
	} else ("x","x")

}

/*
def advPhase(curPhase: Tuple2[String,String], child: Array[String], family: Array[String], pedigree: HashMap[String, Array[String]]): String ={
	if (pedigree.contains())
	}
*/

/* Take phased site from parents and use Homozygous SNPs in Child to drop phase down*/

def childPhase(curPhase: Tuple2[String,String], child: Array[String]): String ={
	val childGT = child(0)
	val minPLval = if (vcfType == "gatk") 40 else 5
	val childPL = child(PL).split(",").map(num => scala.math.abs(num.toFloat)).sorted.tail
	if ((curPhase._1 != "x") && (! childPL.exists(_.toInt <= minPLval)) && (childGT == "1/1" || childGT == "0/0")){
		if (curPhase._1 == childGT(0).toString || curPhase._1 == childGT(2).toString) "S" else "D"
		} else {
		"U"
		}
}
	

/* Main body of Program */

def main (args: Array[String]): Unit = {

	var settings = new HashMap[String,String]
	
	for (items <- args){
		val keyVal = items.split("=")
		settings += keyVal(0).toUpperCase -> keyVal(1) 
	}
	

	if ((! settings.contains("VCF")) && (! settings.contains("PED")) & (! settings.contains("TRIOS")) & (! settings.contains("TYPE")) & (! settings.contains("OUT"))) {
		println("advFilter VCF=input.vcf.gz PED=input.ped TRIOS=input_probands.txt OUT=denovo-stats.txt type=gatk,plat,fb { minDP=0 minALT=0 RECUR=F/T minKIDS=1 minPL=0 QUAL=0 minRAFQ=0.2 }")
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
	val reoccur = if (settings.contains("RECUR") && List("TRUE","YES","Y","T").contains(settings("RECUR").toUpperCase)) true else false
	val QUAL = if (settings.contains("QUAL")) settings("QUAL").toInt else 0
	val minPL = if (settings.contains("MINPL")) settings("MINPL").toInt else 0
	val minKids = if (settings.contains("MINKIDS")) settings("MINKIDS").toInt else 1
	val minRAFreq = if (settings.contains("MINRAFQ")) settings("MINRAFQ").toFloat else 0.2
	val debug = if (settings.contains("DEBUG")) true else false
	settings("TYPE").toUpperCase match {
		case "GATK" => vcfType = "gatk"
		case "PLAT" => vcfType = "platypus"
		case "FB" => vcfType = "freebayes"
		case _ => println("\n\nUnknown VCF type, VCF type must be on of gatk, plat (platypus) or fb (Freebayes)"); System.exit(1)
	}

	val outname = settings("VCF").split("/")(settings("VCF").split("/").size - 1)
	val out_vcf = new BufferedWriter(new OutputStreamWriter(new BlockCompressedOutputStream(outname + ".mutations-" + reoccur + "-denovos.vcf.gz")))
	val out_somatic = new BufferedWriter(new OutputStreamWriter(new BlockCompressedOutputStream(outname + ".mutations-" + reoccur + "-somatic.vcf.gz")))
	val statsOut = new BufferedWriter(new FileWriter(settings("OUT")))
	
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
	var lastPhase = new HashMap[String,Tuple2[String,String]]
	var childState = new HashMap[String,String]
	var allChildren = new HashMap[String,String]
	var childHaps = new HashMap[String, List[String]]
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
			parents = List(pedFile(curPro(0))(2),pedFile(curPro(0))(3))
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
				if (! pedFile.contains(indv)){
					pedFile += indv -> Array("noFAM",indv,"0","0","0")
				}
			}//efor
			var tmpdesc = descendents.filterNot(x => ((x == curPro(0)) || children.contains(x) || ancestors.contains(x) || parents.contains(x)|| (! vcfanimals.contains(x))))
			
			population = animalIDS.toList.filterNot(x => ( (x == curPro(0)) || children.contains(x) || ancestors.contains(x) || descendents.contains(x) || parents.contains(x)|| (! vcfanimals.contains(x))))
			
			if (animalIDS.contains(parents(0)) && animalIDS.contains(parents(1)) && animalIDS.contains(curPro(0)) && children.size >= minKids){
				trios += curPro(0) -> (ancestors, parents, children, tmpdesc, curPro(1).toInt, population, extFam)
				for (child <- children) {
				allChildren += child -> ""
				childHaps += child ->  ("chr1\t1\t1\tU" :: Nil)
				}
			}
		}//eif
	}//ewhile

	in_pro.close

println("Built Pedigrees\n")
//Closed proband file replace println with statsOut.write

/* Report identified Trios & there Pedigree Structure*/

	for (fam <- trios){
		println(s"TRIO:\t${fam._1}\t${fam._2._2(0)}\t${fam._2._2(1)}")
		println("Grandparents\t")
		fam._2._1.foreach(s => println(s + "\t"))
		println(fam._2._1.size)
		println("Children\t")		
		fam._2._3.foreach(s => println(s + "\t"))
		println(fam._2._3.size)
		println("Descendents\t")
		fam._2._4.foreach(s => println(s + "\t"))
		println(fam._2._4.size + "\n")
		println("EXFAM\t")
		fam._2._7.foreach(s => println(s + "\t"))
		println(fam._2._7.size)
		lastPhase += fam._1 -> ("u","u")
		childState += fam._1 -> ""
	}

/*
*	Iterate through VCF file line by line, at each line load each Trio and count existence of variants in different categories
*	if de novo, flag and output snp detail and variant info, (count in pop, children ancestors etc)
*/
	print(s"Chrom\tPos\tRef\tRefSize\tAlt\tAltSize\tQUAL\tTrio\tGenotype\tPLs\tPhase\t Vars S|D Haps S|D\tAnces\tPars\tChildren\tDesc\tExFam\tPop\tPopFreq\tSupport Ratio\tScore\tClass\tProband\tSire\tDam\tPopRefCount\tPopAltCount\tWarning\tPhaseInfo\n")
	var lastChr = ""
	while (in_vcf.ready){
		PL = -1
		PLexist = false
		var denovo = false
		val formatDetails = new HashMap[String,Int]
		var line = in_vcf.readLine().split("\t")
		val format = line(8).split(":")
		if (format.contains("PL") || format.contains("GL")){
			PL = if (format.contains("PL")) format.indexOf("PL") else format.indexOf("GL")
			PLexist = true
		}
		if (lastChr != line(0)){
			lastChr = line(0)
			allChildren.keys.foreach(s => allChildren(s) = "")
			errors.print(line(0))
		}
		
		AD = format.indexOf("AD")
		GT = format.indexOf("GT")
		DP = if (format.contains("NV")) format.indexOf("NR") else format.indexOf("DP")
		AO = if (format.contains("NV")) format.indexOf("NV") else format.indexOf("AO")
		RO = if (format.contains("NV")) format.indexOf("NR") else format.indexOf("RO")
				
		/*To be considered the VCF record must be ok, the Qual score >= Min & no more than 3 alternative alleles*/
//try {		
		if (line.size == (vcfanimals.size + 9) && (line(5).toFloat >= QUAL) && (line(4).split(",").size < 3)){
			var trioPos = 0
			val triosArray = trios.keys.toArray
			
			for (fam <- trios.toArray.par){
			//while(trioPos < triosArray.size){
				//var fam = trios.toArray.apply(trioPos)
				//trioPos += 1
				var indv = ""
				var kidsPhase : List[String] = Nil
	//try {
				var altsPar = 0
				val ped = fam._2
				var ances, par, kids, desc, popFreq, exFamFreq = 0
				val maxDP = (ped._5 * 1.7).toInt
				var adratio = 0.0
				var sirePhase, damPhase = 0
				var allChildrenState = ""
				var parPos, grandPos, childPos, descPos, popRef, popALT, popPos, exfPos = 0
				
				val sire = ped._2.apply(0)
				val dam = ped._2.apply(1)
				
			/* Parental Test using permutations of Alleles */

				val par1 = line(vcfanimals(sire)).split(":")
				val par2 = line(vcfanimals(dam)).split(":")
				val proBand = line(vcfanimals(fam._1)).split(":")
			if (isVar(proBand(GT))){

				if (par1(GT)(0) != '.' && par2(GT)(0) != '.' && proBand(GT)(0) != '.'){
					var phasVal = if (checkDP(proBand, DP, minDP, maxDP) && checkDP(par1,DP,minDP,maxDP) && checkDP(par2,DP,minDP,maxDP)) phase(proBand,par1, par2) else ("x","x")
					val valGTs = permu(par1(GT)(0).toString + par1(GT)(2),par2(GT)(0).toString + par2(GT)(2))
					if (valGTs.contains(proBand(GT)(0).toString + proBand(GT)(2))){
						par += 1
					}

			/* After checking that the Parent GT's can produce the Denovo check AD/RO/AO incase GT misscall */
				while(parPos < ped._2.size){
					indv = ped._2.apply(parPos)
					parPos += 1
					//for (indv <- ped._2){
						if (line(vcfanimals(indv))(0) != '.'){
							val curAn = line(vcfanimals(indv)).split(":")
							val refAlt = selROvAD(curAn,AD, RO, AO, GT)
							altsPar += refAlt._2
							if (refAlt._2 != -1 ) adratio = rc(refAlt._1,refAlt._2)
							if (sigAD(refAlt._2)){
								par += 1
							}
						}
					}

			/* Loop through each family group and record Hets */
				var grands: List[String] = Nil
				while(grandPos < ped._1.size){
					indv = ped._1.apply(grandPos)
					grandPos += 1
					//for (indv <- ped._1){
						if (line(vcfanimals(indv))(0) != '.'){
							val curAn = line(vcfanimals(indv)).split(":")
							val refAlt = selROvAD(curAn,AD, RO, AO, GT)
							if (refAlt._2 != -1 ) adratio = rc(refAlt._1,refAlt._2) 
							grands = grands ++ rtGTs(curAn(GT))
							if (sigAD(refAlt._2)){
								ances += 1
							}
						}
					}
					
					if (proBand(GT)(0) != '.' && grands.contains(proBand(GT)(0)) && grands.contains(proBand(GT)(1))){
						ances += 1
					}

			//Children
				while(childPos < ped._3.size){
					indv = ped._3.apply(childPos)
					childPos += 1
					//for (indv <- ped._3){
						val curAn = line(vcfanimals(indv)).split(":")
						if (line(vcfanimals(indv))(0) != '.'){
							var inherited = childPhase(phasVal,curAn)
							val sireid = pedFile(indv).apply(2)
							val damid = pedFile(indv).apply(3)
							if (vcfanimals.contains(sireid) && vcfanimals.contains(damid)) {
								val cSire = line(vcfanimals(sireid)).split(":")
								val cDam = line(vcfanimals(damid)).split(":")
								//val fullPhase = if (checkDP(curAn, DP, minDP, maxDP) && checkDP(cSire,DP,minDP,maxDP) && checkDP(cDam,DP,minDP,maxDP)) phase(curAn,cSire, cDam) else ("x","x")
								val fullPhase = phase(curAn, cSire, cDam)
								if (fullPhase._1 != "x") { 
								if (fullPhase._1 == phasVal._1) inherited = "S" else inherited = "D" 							
								} else {
									inherited = childPhase(phasVal,curAn)
								}
							} else {
								inherited = childPhase(phasVal,curAn)
							}
							val refAlt = selROvAD(curAn,AD, RO, AO, GT)					
							if (isVar(curAn(GT)) || sigAD(refAlt._2)){
								if (inherited != "U" && (checkDP(curAn,DP,minDP,maxDP))) {
								val parID = if (inherited == "S") sire else dam
									allChildren(indv) = parID
									childHaps(indv) = s"${line(0)}\t${line(1)}\t${line(1)}\t${parID}" :: childHaps(indv)
								}
								
								kids += 1
								kidsPhase = allChildren(indv) :: kidsPhase
							}
						}
						if (allChildren(indv) == sire) sirePhase += 1 else damPhase += 1
						allChildrenState = allChildrenState + s"${indv}:${if (curAn.size >= PL ) curAn(PL) else 0}:${curAn(0)}:${allChildren(indv)}\t"
					}

			//Desec
				while(descPos < ped._4.size){
					indv = ped._4.apply(descPos)
					descPos += 1
			//		for (indv <- ped._4){
						if (line(vcfanimals(indv))(0) != '.'){
							val curAn = line(vcfanimals(indv)).split(":")
							val refAlt = selROvAD(curAn,AD, RO, AO, GT)
							if (isVar(curAn(GT)) || sigAD(refAlt._2)){
								desc += 1
							}
						}
					}

			/* Population Calc */
					
					while (popPos < ped._6.size){
						val indv = ped._6.apply(popPos)
						popPos += 1
					//for (indv <- ped._6){
						if (line(vcfanimals(indv))(0) != '.'){
							val curAn = line(vcfanimals(indv)).split(":")
							val refAlt = selROvAD(curAn,AD, RO, AO, GT)
							if (isVar(curAn(GT)) || sigAD(refAlt._2)){
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

					while (exfPos < ped._7.size){
						val indv = ped._7.apply(exfPos)
						exfPos += 1
					//for (indv <- ped._7){
						if (line(vcfanimals(indv))(0) != '.'){
							val curAn = line(vcfanimals(indv)).split(":")
							val refAlt = selROvAD(curAn,AD, RO, AO, GT)
							if ((isVar(curAn(GT)) || sigAD(refAlt._2)) && checkDP(curAn,DP,minDP,maxDP)){
								exFamFreq += 1
							}
						}
					}

			/* -------- Denovo Check ---------- */

					val curPro = line(vcfanimals(fam._1)).split(":")
					var rank = ""
					altsPar match {
						case 0 => rank = "probable"
						case 1 => rank = "possible"
						case _ => rank = "unlikely"
					}
					
					val proRatio = selROvAD(proBand,AD, RO, AO, GT)
					val proGT = proBand(GT)
					
					if (debug && popFreq == 0 & par == 0) {
					proBand.foreach(s => print(s"${s} "))
					 print("\n")
					 print(s"GT ${!valGTs.contains(proBand(GT)(0).toString + proBand(GT)(2))} ALTGood(${proRatio._2 >= minALT} ${proRatio._2 == -1}) PL${checkPL(minPL, curPro)} SPL${checkPL(minPL, par1)} DPL${checkPL(minPL, par2)}\n")
					 print(s"DP${checkDP(curPro, DP, minDP, maxDP)}  SDP${checkDP(par1, DP, minDP, maxDP)} DDP${checkDP(par2, DP, minDP, maxDP)} AN${ances == 0} PAR${par == 0} KID${kids >= minKids} RATIO${(proRatio._2/proRatio._1.toFloat) >= minRAFreq} pop${popFreq} \n")
					 }
					 
					/*
					* De novo Identification logic!
					* Proband is Variant, has > min num of Alt alleles. DP is between minDP & maxDP
					* Parental DP is good, variant not present in Ancestors, Parents & is present in at least
					* minKids kids & min ratio between Ref & Alt is fine
					* If recurring allowed don't check population otherwise not present in population
					*/
					
					if (
							(
								(
									!valGTs.contains(proBand(GT)(0).toString + proBand(GT)(2)
								) && 	(
										proRatio._2 >= minALT || proRatio._2 == -1) || 
											
												proRatio._2 >= (minALT * 3)
											)
										) && checkPL(minPL, curPro) && checkPL(minPL, par1) && checkPL(minPL, par2)
											&& checkDP(curPro, DP, minDP, maxDP) &&  checkDP(par1,DP,minDP,maxDP) && checkDP(par2,DP,minDP,maxDP) &&
											(ances == 0) && (par == 0) && (kids >= minKids) && 
											(
											 	(proRatio._2/proRatio._1.toFloat) >= minRAFreq
											 ) && 
											 (
											 	if (reoccur){true}else{popFreq == 0}
											 )
							){
						//denovo = true
						
						var phaseQual = ""

						if (kidsPhase.exists(_ == sire) && kidsPhase.exists(_ == dam)){
							phaseQual = "Bad\t" + kidsPhase.count(_ == sire) + "|" + kidsPhase.count(_ == dam) + " " + sirePhase + "|" + damPhase
						} else {
							if ((kidsPhase.exists(_ == sire) && (kidsPhase.count(_ == sire) != sirePhase)) || (kidsPhase.exists(_ == dam) && (kidsPhase.count(_ == dam) != damPhase))){
								phaseQual = "Partial\t" + kidsPhase.count(_ == sire) + "|" + kidsPhase.count(_ == dam) + " " + sirePhase + "|" + damPhase
							} else {
								phaseQual = "Good\t" + kidsPhase.count(_ == sire) + "|" + kidsPhase.count(_ == dam) + " " + sirePhase + "|" + damPhase
							}
						}
						
						/*
						if ((proRatio._1 + proRatio._2) <= minDP) {
							errors.println(s"minDP == ${minDP}\t${proRatio._1} + ${proRatio._2}\t ${proBand(DP)}")
							proBand.foreach(er => errors.print(er + " "))
							errors.print("\n")
							}
							*/
						
						if (reoccur && adratio == 0.0){
							print(s"${line(0)}\t${line(1)}\t${line(3)}\t${line(3).size}\t${line(4)}\t${line(4).size}\t${line(5)}\t${fam._1}\t'${proGT}\t${if (PLexist) proBand(PL) else -1}\t${phaseQual}\t${ances}\t${par}\t${kids}\t${desc}\t${exFamFreq}\t${popFreq}\t${popFreq.toFloat/(animalIDS.size)}\t${proRatio._2/proRatio._1.toFloat}\t${rank}\tdenovo\t" + 
								proRatio + "\t" + selROvAD(par1,AD, RO, AO, GT) + " " + (if (PLexist) par1(PL) else "0,0,0") + "\t" + selROvAD(par2,AD, RO, AO, GT) + " " + (if (PLexist) par2(PL) else "0,0,0") + s"\t${popRef}\t${popALT}\t\t${allChildrenState}\n")
							out_vcf.write(line.reduceLeft{(a,b) => a + "\t" + b} + "\n")
						}else {
							print(s"${line(0)}\t${line(1)}\t${line(3)}\t${line(3).size}\t${line(4)}\t${line(4).size}\t${line(5)}\t${fam._1}\t'${proGT}\t${if (PLexist) proBand(PL) else -1}\t${phaseQual}\t${ances}\t${par}\t${kids}\t${desc}\t${exFamFreq}\t${popFreq}\t${popFreq.toFloat/(animalIDS.size)}\t${proRatio._2/proRatio._1.toFloat}\t${rank}\tdenovo\t" + 
							proRatio + "\t" + selROvAD(par1,AD, RO, AO, GT) + " " + (if (PLexist) par1(PL) else "0,0,0") + "\t" + selROvAD(par2,AD, RO, AO, GT) + " " + (if (PLexist) par2(PL) else "0,0,0") +  s"\t${popRef}\t${popALT}")
							if (adratio != 0.0) {
								line(6) = "LOWQUAL_ADratio"
								print(s"\t WARNING: Low confidence de novo\t${allChildrenState}\n")
							}else {
								print("\n")
							}//eelse
							out_vcf.write(line.reduceLeft{(a,b) => a + "\t" + b} + "\n")
						}//eif reoccur

					} //eif is Denovo
					
					
					if(!valGTs.contains(proBand(GT)(0).toString + proBand(GT)(2)) && (ances == 0) && (par == 0) && (kids == 0) 
						&& (popFreq == 0) && checkDP(curPro, DP, minDP, maxDP) && checkDP(par1,DP,minDP,maxDP) && checkDP(par2,DP,minDP,maxDP)
						&& proRatio._2 >= minALT && adratio == 0.0){
							print(s"${line(0)}\t${line(1)}\t${line(3)}\t${line(3).size}\t${line(4)}\t${line(4).size}\t${line(5)}\t${fam._1}\t'${proGT}\t${if (PLexist) proBand(PL) else -1}\t\t${ances}\t${par}\t${kids}\t${desc}\t${exFamFreq}\t${popFreq}\t${popFreq.toFloat/(animalIDS.size)}\t${proRatio._2/proRatio._1.toFloat}\t${rank}\tSomatic\t" +
							 proRatio + "\t" + selROvAD(par1,AD, RO, AO, GT) + " " + (if (PLexist) par1(PL) else "0,0,0") + "\t" + selROvAD(par2,AD, RO, AO, GT) + " " + (if (PLexist) par2(PL) else "0,0,0") + s"\t${popRef}\t${popALT}\t\n")
							out_somatic.write(line.reduceLeft{(a,b) => a + "\t" + b} + "\n")
						}
					
					

				}//eisVAR
			} // IS VAR
			}//Efor fam <- trios
		}//IF Valid VCF line
	}// Ewhile
	in_vcf.close
	out_vcf.close
	
	for (child <- childHaps){
		val childOrigin = child._2.reverse.toArray
		errors.print(childOrigin.size + "\n")
		val out = new BufferedWriter(new FileWriter(child._1 + "-origin.bed"))
		var count = 1
		var lastPos, lastChr, lastOrig = "" 
		while (count < childOrigin.size){
			val cline = childOrigin(count).split("\t")
			if (lastChr == ""){
				lastChr = cline(0)
				lastPos = cline(1)
				lastOrig = cline(3)
				out.write(s"${cline(0)}\t${cline(1)}\t")
			}
				if ((cline(3) != lastOrig) || (lastChr != cline(0))){
					out.write(s"${lastPos}\t${lastOrig}\n")
					lastPos = cline(1)
					lastChr = cline(0)
					lastOrig = cline(3)
					out.write(s"${cline(0)}\t${cline(1)}\t")
				} else {
					lastPos = cline(1)
				}
			
			count += 1

		}
		out.write(s"${lastPos}\t${lastOrig}\n")
		//child._2.reverse.foreach(s => out.write(s + "\n"))
		out.close
	}
	
}//eMain

}//eObject