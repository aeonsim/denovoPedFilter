/* Advanced Pedigree Filter
* Takes in VCF (gz) & Pedigree file, with cutoffs for DP, and a list of Proband ID's
* Use pedigree file to build family tree for each Proband, Proband, Children, Descendants, Parents, Ancestors
* Select variants based on pedigree & min DP check AD for low coverage
* Parents & Ancestors -ve, Proband +ve, Children 1+ positive, Descendants possibly +ve
* Output VCF per Proband ID of possible De novos along with, score file with counts in the population.
* Check to see if each VCF line is the correct length & discard if not.
*/

/* phaseTracker takes in phase blocks and tracks position within the block for each query*/
class phaseTracker (phaseData: Array[Tuple4[String,Int,Int,String]]){
	var curPhaseBlock = 0
	
	def getPhase(chrom: String, position: Int) : String ={
		println("Block " + phaseData(curPhaseBlock)._1 + " "  + phaseData(curPhaseBlock)._2 + " " + phaseData(curPhaseBlock)._3 + " " + phaseData(curPhaseBlock)._4 + " Looking For:" + chrom + " " + position)
		
		 if (phaseData(curPhaseBlock)._1 != chrom) {
		 	val chr = if (chrom(0).toString.toUpperCase == "C") chrom.substring(3,chrom.size) else chrom
		 	val chrnum = if (chr.toUpperCase == "X") 999 else if (List("M","MT").contains(chr.toUpperCase)) 1000 else chr.toInt
		 	while (phaseData(curPhaseBlock)._1 != chrom){
				val bchr = if( phaseData(curPhaseBlock)._1(0).toString.toUpperCase == "C" ) phaseData(curPhaseBlock)._1.substring(3, phaseData(curPhaseBlock)._1.size) else  phaseData(curPhaseBlock)._1
				val bchrnum = if (bchr.toUpperCase == "X") 999 else if (List("M","MT").contains(bchr.toUpperCase)) 1000 else bchr.toInt
				if (chrnum >= bchrnum ) curPhaseBlock += 1 else curPhaseBlock -= 1
				System.err.println("looping through Chrs currently " + phaseData(curPhaseBlock)._1 + " Going to " + chrom )
			}
		}
		val nextPhaseBlock = curPhaseBlock + 1
		val curStart = phaseData(curPhaseBlock)._2
		val curEnd = phaseData(curPhaseBlock)._3
		
		if (nextPhaseBlock < phaseData.size){
			val nextStart = phaseData(nextPhaseBlock)._2
			val nextEnd = phaseData(nextPhaseBlock)._3
			val nextChrom = phaseData(nextPhaseBlock)._1
			
			if (position <= curEnd && position >= curStart){
				return phaseData(curPhaseBlock)._4
			} else {
				if (position < curStart) {
					return "BOTH"
				} else {
					if (position > curEnd && position < nextStart){
						return "BOTH"
					} else {
						if (position <= nextEnd && position >= nextStart){
							curPhaseBlock += 1
							return phaseData(curPhaseBlock)._4
						} else {
							if (position > nextEnd && nextChrom == chrom){
								curPhaseBlock += 1
								println("FWD Block")
								getPhase(chrom, position)
							} else {
								//curPhaseBlock += 1
								return "BOTH"
							}
						}
					}
				}				
			}
		} else {
			if (position < curStart){
				return "BOTH"
			} else {
				if (position <= curEnd && position >= curStart){
					return phaseData(curPhaseBlock)._4
				} else {
					return "BOTH"
				}
			}
		}
	}
}


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


	def itPed (pop: HashMap[String, Array[String]] , rec: String, maxIT: Int): List[String] = {
		if (pop.contains(rec) && maxIT >= 0){
			return rec :: itPed(pop,pop(rec).apply(2),maxIT -1) ::: itPed(pop,pop(rec).apply(3),maxIT -1)
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
		try{
		if (DPpos != -1 && genos.size > DPpos  && genos(DPpos) != "."){
			val curDP = if (vcfType == "platypus") genos(DPpos).split(",")(0).toInt else genos(DPpos).toInt
			if (curDP >= minDP && curDP <= maxDP) true
			else false
		} else {
			false
		}
		} catch {
			 case e: Exception => System.err.println(e + " DPpos " + DPpos + " Geno " + genos.reduceLeft{(a,b) => a + ":" + b})
			 false
		}
	}
	
	/* Check PL is sufficient, if no PL's then aways ok*/
	
	def checkPL(min: Int, indv: Array[String]): Boolean = {
		if (PLexist){
			val pls = indv(PL).split(",").sorted.tail
			if (pls(0).toInt >= min) true 
			else false
		} else {
			true
		}	
	}
	
	
/* Phase Code, return format is (sireAllele,damAllele)*/

def phase(indv: Array[String], sire: Array[String], dam: Array[String]) : Tuple2[String, String] = {
if (sire.size >= PL && dam.size >= PL && indv.size >= PL) {
	val minPLval = if (vcfType == "gatk") 40 else 5
	val indvGT = indv(0)
	val sireGT = sire(0)
	val damGT  = dam(0)
	val indvPL = indv(PL).split(",").sorted.tail 
	val sirePL = sire(PL).split(",").sorted.tail
	val damPL  = dam(PL).split(",").sorted.tail
if ((indvPL(0).toInt >= minPLval && sirePL(0).toInt >= minPLval && damPL(0).toInt >= minPLval)){
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

/* Take phased site from parents and use Homozygous SNPs in Child to drop phase down*/

def childPhase(curPhase: Tuple2[String,String], child: Array[String]): String ={
	if (child.size >= PL){
	val childGT = child(0)
	val minPLval = if (vcfType == "gatk") 40 else 5
	val childPL = child(PL).split(",").sorted.tail 
	if ((curPhase._1 != "x") && (! childPL.exists(_.toInt <= minPLval)) && (childGT == "1/1" || childGT == "0/0")){
		if (curPhase._1 == childGT(0).toString || curPhase._1 == childGT(2).toString) "S" else "D"
		} else {
		"U"
		}
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
	

	if ((! settings.contains("VCF")) && (! settings.contains("PED")) && (! settings.contains("TRIOS")) && (! settings.contains("REF")) && (! settings.contains("TYPE")) && (! settings.contains("OUT"))) {
		println("advFilter VCF=input.vcf.gz PED=input.ped TRIOS=input_probands.txt OUT=denovo-stats.txt type=gatk,plat,fb ref=genome.fasta { phaseVCF=SNPChip.vcf.gz minDP=0 minALT=0 RECUR=F/T minKIDS=1 PLGL=0 QUAL=0 minRAFQ=0.2 }")
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
	val minPL = if (settings.contains("minPL")) settings("minPL").toInt else 0
	val minKids = if (settings.contains("MINKIDS")) settings("MINKIDS").toInt else 1
	val minRAFreq = if (settings.contains("MINRAFQ")) settings("MINRAFQ").toFloat else 0.2
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
	//var childHaps = new HashMap[String, List[String]]
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
	//out_somatic.write(s"##${args.reduceLeft{(a,b) => a + " " + b}}\n")

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
			ancestors = itPed(pedFile,curPro(0),2).tail.filterNot(x => parents.contains(x) || (! vcfanimals.contains(x)))
			children = findChildren(pedFile,vcfanimals,curPro(0))
			for (indv <- animalIDS.toList){ //.filterNot(x => ( (x == curPro(0)) || children.contains(x) || ancestors.contains(x) || parents.contains(x)))){
				val tempPed = itPed(pedFile,indv,99)
				if (tempPed.contains(parents(0)) || tempPed.contains(parents(1))){
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
			
			if (animalIDS.contains(parents(0)) && animalIDS.contains(parents(1)) && animalIDS.contains(curPro(0)) && children.size != 0){
				trios += curPro(0) -> (ancestors, parents, children, tmpdesc, curPro(1).toInt, population, extFam)
				for (child <- children) {
				allChildren += child -> ""
				//childHaps += child ->  ("chr1\t1\t1\tU" :: Nil)
				}
			}
		}//eif
	}//ewhile

	in_pro.close

println("Built Pedigrees\n")
//Closed proband file replace println with statsOut.write

/* Report identified Trios & there Pedigree Structure*/

	for (fam <- trios){
		println(s"TRIO:\t${fam._1}\t${fam._2._2(0)}\t${fam._2._2(1)}\n")
		println("Grandparents\t")
		fam._2._1.foreach(s => println(s + "\t"))
		println(fam._2._1.size + "\n")
		println("Children\t")		
		fam._2._3.foreach(s => println(s + "\t"))
		println(fam._2._3.size + "\n")
		println("Descendents\t")
		fam._2._4.foreach(s => println(s + "\t"))
		println(fam._2._4.size + "\n")
		println("EXFAM\t")
		fam._2._7.foreach(s => println(s + "\t"))
		println(fam._2._7.size + "\n\n")
		lastPhase += fam._1 -> ("u","u")
		childState += fam._1 -> ""
	}

/*
* PrePhase Data using SNPChip data if available
*/

var phaseTracking = new HashMap[String, phaseTracker]
var readDepths = new HashMap[String,Array[Int]]
var phaseBlock = new HashMap[String,List[Tuple4[String,Int,Int,String]]]
//var tmpPhase = new HashMap[String,List[Tuple3[String,Int,String]]]

	val phaseInfo = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream(settings("VCF")))))
	System.err.println("Have VCF, beginning Pre-phasing")
	var phaseLine = phaseInfo.readLine.split("\t")
	while (phaseLine(0).apply(1) == '#') phaseLine = phaseInfo.readLine.split("\t")
	for (col <- 9 to (phaseLine.size -1)){
		readDepths += phaseLine(col) -> new Array[Int](52)
		phaseBlock += phaseLine(col) -> Nil
	}
	
	/* Will output RAW data to disk rather than trying to store in memory */
	
	var rawOutput = new HashMap[String,BufferedWriter]
	for (fam <- trios){
	/* Record Sites that are decent in Trio */
		readDepths += s"Trio_${fam._1}" -> new Array[Int](4)
		for (kid <- fam._2._3){
			//tmpPhase += kid -> (("",0,"") :: Nil)
			rawOutput += kid -> new BufferedWriter(new FileWriter(kid + "-rawPhase.txt"))
		} 
	}
	
	while(phaseInfo.ready){
	
		phaseLine = phaseInfo.readLine.split("\t")
	
		val format = phaseLine(8).split(":")
		val GT = format.indexOf("GT")
		DP = if (format.contains("NV")) format.indexOf("NR") else format.indexOf("DP")
		if (format.contains("PL") || format.contains("GL")){
			PL = if (format.contains("PL")) format.indexOf("PL") else format.indexOf("GL")
			PLexist = true
		}
		
		for (fam <- trios.par){
			/* Family = (ancestors, parents, children, tmpdesc, curPro(1).toInt, population, extFam) */
			val family = fam._2
			val maxDP = (family._5 * 1.7).toInt
				val proband = phaseLine(vcfanimals(fam._1)).split(":")
				val sire = phaseLine(vcfanimals(family._2(0))).split(":")
				val dam = phaseLine(vcfanimals(family._2(1))).split(":")
				
				//println(proband.reduceLeft{(a,b) => a + ":" + b})
				//println(phaseLine(8) + " DP=" + DP)
				//val pDP = if (proband.size > DP && proband(DP) != ".") proband(DP).toInt else 0
				val pDP = if (proband.size > DP && DP != -1 && proband(DP) != ".") proband(DP).toInt else 0
				val sDP = if (sire.size > DP && DP != -1 && sire(DP) != "." ) sire(DP).toInt else 0
				val dDP = if (dam.size > DP && DP != -1 && dam(DP) != ".") dam(DP).toInt else 0
				
				if (pDP >= 51) readDepths(fam._1)(51) += 1 else readDepths(fam._1)(pDP) += 1
				if (sDP >= 51) readDepths(family._2(0))(51) += 1 else readDepths(family._2(0))(sDP) += 1
				if (dDP >= 51) readDepths(family._2(1))(51) += 1 else readDepths(family._2(1))(dDP) += 1
								
				/* Phase Trio at Site 0 = Good Autozome, 1 = bad Autozome, 2 = Good X, 3 = Bad X */
				val phaseVal = if (checkDP(proband, DP, minDP, maxDP) && checkDP(sire,DP,minDP,maxDP) && checkDP(dam,DP,minDP,maxDP)) {
						if (proband(0).toUpperCase != "CHRX") readDepths(s"Trio_${fam._1}")(0) += 1 else readDepths(s"Trio_${fam._1}")(2) += 1 
						phase(proband, sire, dam) 
					} else { 
						if (proband(0).toUpperCase != "CHRX") readDepths(s"Trio_${fam._1}")(1) += 1 else readDepths(s"Trio_${fam._1}")(3) += 1 
						("x","x")
					}
						
				for (kid <- family._3){
					var curKid = phaseLine(vcfanimals(kid)).split(":")
					val kidDP = if (curKid.size > DP && DP != -1 && curKid(DP) != ".") curKid(DP).toInt else 0
					if (kidDP >= 51) readDepths(kid)(51) += 1 else readDepths(kid)(kidDP) += 1
					
					var inherited = childPhase(phaseVal,curKid)
					val sireid = pedFile(kid).apply(2)
					val damid = pedFile(kid).apply(3)
					if (vcfanimals.contains(sireid) && vcfanimals.contains(damid)) {
					
						val cSire = phaseLine(vcfanimals(sireid)).split(":")
						val cDam = phaseLine(vcfanimals(damid)).split(":")
						
						val fullPhase = phase(curKid, cSire, cDam)
						
						if (fullPhase._1 != "x") { 
							if (fullPhase._1 == phaseVal._1) inherited = "S" else inherited = "D" 							
						}
					} 
					
					if (inherited != "U") {
						val parID = if (inherited == "S") family._2(0) else family._2(1)
						rawOutput(kid).write(s"${phaseLine(0)}\t${phaseLine(1)}\t${parID}\n")
						//tmpPhase(kid) =  (phaseLine(0),phaseLine(1).toInt,parID) :: tmpPhase(kid)
					}
				}
		}
	} //While Phasing
	
	for (fileout <- rawOutput){
		fileout._2.close
	}
	
	//for (fam <- trios){
		//for (kid <- fam._2._3){
			//tmpPhase(kid) = tmpPhase(kid).reverse.tail
		//rawOutput.foreach(child => child._2.close)
		//} 
	//}
	
	phaseInfo.close
	
	for (indv <- readDepths){
		val out = new BufferedWriter(new FileWriter(indv._1 + "-depths.txt"))
		if (indv._2.size > 5){
			val sum = indv._2.sum
			for (num <- 0 to 51){
				out.write(s"${num}\t${indv._2(num)}\t${indv._2(num)/sum.toFloat * 100}\n")
			}
		} else {
			out.write(s"Autozome good: ${indv._2(0)}\t${indv._2(0)/(indv._2(0) + indv._2(1) + 1) * 100}\n")
			out.write(s"Autozome bad : ${indv._2(1)}\t${indv._2(1)/(indv._2(0) + indv._2(1) + 1) * 100}\n")
			out.write(s"ChrX good    : ${indv._2(2)}\t${indv._2(2)/(indv._2(2) + indv._2(3) + 1) * 100}\n")
			out.write(s"ChrX bad     : ${indv._2(3)}\t${indv._2(3)/(indv._2(2) + indv._2(3) + 1) * 100}\n")
		}
		out.close
	}
	
		//for (child <- tmpPhase.par){
		for(child <- rawOutput){
			val in = new BufferedReader(new FileReader(child._1 + "-rawPhase.txt"))
			val out = new BufferedWriter(new FileWriter(child._1 + "-origin.bed"))	
			
			var phased : List[Tuple3[String,Int,String]] = Nil 
			
			while (in.ready){
				/* Read file into Array of tuples */
				val line = in.readLine.split("\t")
				phased = (line(0),line(1).toInt,line(2)) :: phased
			}
			in.close
						
			val childOrigin = phased.reverse.toArray
			
			var count = 0			
			var chrom = childOrigin(count)._1
			var parent = childOrigin(count)._3
			var start = childOrigin(count)._2
			var end = childOrigin(count)._2

			
			phaseBlock += child._1 -> Nil
			
			while (count < (childOrigin.size - 1)){
				count +=1
				if (parent != childOrigin(count)._3 | chrom != childOrigin(count)._1){
					phaseBlock(child._1) = Tuple4(chrom,start,end,parent) :: phaseBlock(child._1)
					out.write(s"${chrom} ${start} ${end} ${parent}\n")
					chrom = childOrigin(count)._1
					parent = childOrigin(count)._3
					start = childOrigin(count)._2
					end = childOrigin(count)._2
				} else {
					end = childOrigin(count)._2
				}				
			}
			phaseBlock(child._1) = Tuple4(chrom,start,end,parent) :: phaseBlock(child._1)
			out.write(s"${chrom} ${start} ${end} ${parent}\n")
			out.close
			
			phaseTracking += child._1 -> new phaseTracker(phaseBlock(child._1).reverse.toArray)
		}

readDepths = new HashMap[String,Array[Int]]
phaseBlock = new HashMap[String,List[Tuple4[String,Int,Int,String]]]
//tmpPhase = new HashMap[String,List[Tuple3[String,Int,String]]]

/* Build Reference Array, Note all positions are exact S added as a Guard for 0 */

System.err.println("Reading Reference")
import scala.collection.mutable.ListBuffer
val ref = new BufferedReader(new FileReader(settings("REF")))
//var refTable = new HashMap[String,ListBuffer[Char]]
var refTable = new HashMap[String,Array[Char]]
var line, refLastChr = ""
var refSeq = ListBuffer("Init")

while (ref.ready){
	line = ref.readLine
	if (line(0) == '>'){
		System.err.println(line)
		refTable += refLastChr -> refSeq.reduceLeft{(a,b) => a + b}.toArray
		refSeq = new ListBuffer[String]
		refLastChr = if (line.indexOf(' ') == -1) line.substring(1,line.size) else line.substring(1).split("\\s+")(0)
	} else {
		refSeq.append(line)
	}
}

refSeq = new ListBuffer[String]
ref.close

/*
*	Iterate through VCF file line by line, at each line load each Trio and count existence of variants in different categories
*	if de novo, flag and output snp detail and variant info, (count in pop, children ancestors etc)
*/
	//statsOut.write(s"Chrom\tPos\tRef\tRefSize\tAlt\tAltSize\tQUAL\tTrio\tGenotype\tPLs\tPhase\t Vars S|D Haps S|D\tAnces\tPars\tChildren\tDesc\tExFam\tPop\tPopFreq\tSupport Ratio\tScore\tClass\tProband\tSire\tDam\tPopRefCount\tPopAltCount\tWarning\tPhaseInfo\n")
	println(s"Chrom\tPos\tTRI-NUC\tRef\tRefSize\tAlt\tAltSize\tQUAL\tTrio\tGenotype\tPLs\tPhase\t Vars S|D Haps S|D\tAnces\tPars\tChildren\tDesc\tExFam\tPop\tPopFreq\tSupport Ratio\tScore\tClass\tProband\tSire\tDam\tPopRefCount\tPopAltCount\tWarning\tPhaseInfo")
	var lastChr = ""
	
	System.err.println("Phasing Complete, beginning De Novo identification & Characterisation")
	
	/* Store Output data in List so can loop through after and fix phase */
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
				var sirePhase, damPhase = 0.0
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
					//var phasVal = if (checkDP(proBand, DP, minDP, maxDP) && checkDP(par1,DP,minDP,maxDP) && checkDP(par2,DP,minDP,maxDP)) phase(proBand,par1, par2) else ("x","x")
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
			var varSirePhase, varDamPhase = 0.0
				while(childPos < ped._3.size){
					indv = ped._3.apply(childPos)
					childPos += 1
					//for (indv <- ped._3){
						val curAn = line(vcfanimals(indv)).split(":")
						//if (line(vcfanimals(indv))(0) != '.'){
						if(curAn(0).apply(0) != '.'){	
							/* Remove Phase code belongs in own prePhase block now
							//var inherited = childPhase(phasVal,curAn)
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
							*/
							val refAlt = selROvAD(curAn,AD, RO, AO, GT)			
							if (isVar(curAn(GT)) || sigAD(refAlt._2)){
							allChildren(indv) = phaseTracking(indv).getPhase(line(0),line(1).toInt)	// Maybe should be in isVar test?	
							
								/* More Phase Code
								if (inherited != "U" && (checkDP(curAn,DP,minDP,maxDP))) {
								val parID = if (inherited == "S") sire else dam
									allChildren(indv) = parID
									childHaps(indv) = s"${line(0)}\t${line(1)}\t${line(1)}\t${parID}" :: childHaps(indv)
								}
								*/
								kids += 1
								allChildren(indv) match {
									case `sire` => varSirePhase += 1
									case `dam` => varDamPhase += 1
									case "BOTH" => varSirePhase += 0.01 ; varDamPhase += 0.01
									case _ => System.err.println("##############\n Error in Phase Identification\n" + allChildren(indv) + "\n" + curAn.reduceLeft{(a,b) => a + ":" + b} +"\n############")
								}
								//kidsPhase = allChildren(indv) :: kidsPhase
							}
						}
						allChildren(indv) match {
							case `sire` => sirePhase += 1
							case `dam` => damPhase += 1
							case "BOTH" => sirePhase += 0.01 ; damPhase += 0.01
							case _ => System.err.println("##############\n Error in Phase Identification\n" + allChildren(indv) + "\n" + curAn.reduceLeft{(a,b) => a + ":" + b} +"\n############")
						}
						//if (allChildren(indv) == sire) sirePhase += 1 else damPhase += 1
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
					
					//println("Parents " + par + "Ancestors " + ances + " Kids " + kids + " Desc " + desc + " Pop " + popFreq)
					val debug = false
					
					if (debug && popFreq == 0 && par == 0) {
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
									) && (
											proRatio._2 >= minALT || proRatio._2 == -1) ||
											proRatio._2 >= (minALT * 3)
										)
										) && checkPL(minPL, curPro) && checkPL(minPL, par1) && checkPL(minPL, par2)
											&& checkDP(curPro, DP, minDP, maxDP) && checkDP(par1,DP,minDP,maxDP) && checkDP(par2,DP,minDP,maxDP) &&
												(ances == 0) && (par == 0) && (kids >= minKids) &&
												(
													(proRatio._2/proRatio._1.toFloat) >= minRAFreq
													) &&
														(
															if (reoccur){true}else{popFreq == 0}
															)
						){
	
						var phaseQual = ""
						
						if (varSirePhase >= 1 && varDamPhase >= 1){
							phaseQual = "Bad\t" + varSirePhase + "|" + varDamPhase + " " + sirePhase + "|" + damPhase
						} else {
							if ((varSirePhase > 0 && varSirePhase != sirePhase) || (varDamPhase > 0 & varDamPhase != damPhase)){
								phaseQual = "Partial\t" + varSirePhase + "|" + varDamPhase + " " + sirePhase + "|" + damPhase
							} else {
								phaseQual = "Good\t" + varSirePhase + "|" + varDamPhase + " " + sirePhase + "|" + damPhase
							}
						}

						println(line(0) + " chr size " + refTable(line(0)).size)
						var triNuc = refTable(line(0)).apply(line(1).toInt -1).toString + refTable(line(0)).apply(line(1).toInt).toString + refTable(line(0)).apply(line(1).toInt +1).toString
						var triNucAlt = triNuc(0) + line(4) + triNuc(2)
						
						if (reoccur && adratio == 0.0){
							print(s"${line(0)}\t${line(1)}\t${triNuc}>${triNucAlt}\t${line(3)}\t${line(3).size}\t${line(4)}\t${line(4).size}\t${line(5)}\t${fam._1}\t'${proGT}\t${if (PLexist) proBand(PL) else -1}\t${phaseQual}\t${ances}\t${par}\t${kids}\t${desc}\t${exFamFreq}\t${popFreq}\t${popFreq.toFloat/(animalIDS.size)}\t${proRatio._2/(proRatio._1.toFloat + proRatio._2)}\t${rank}\tdenovo\t" + 
								proRatio + "\t" + selROvAD(par1,AD, RO, AO, GT) + " " + (if (PLexist) par1(PL) else "0,0,0") + "\t" + selROvAD(par2,AD, RO, AO, GT) + " " + (if (PLexist) par2(PL) else "0,0,0") + s"\t${popRef}\t${popALT}\t\t${allChildrenState}\n")
							out_vcf.write(line.reduceLeft{(a,b) => a + "\t" + b} + "\n")
						}else {
							print(s"${line(0)}\t${line(1)}\t${triNuc}>${triNucAlt}\t${line(3)}\t${line(3).size}\t${line(4)}\t${line(4).size}\t${line(5)}\t${fam._1}\t'${proGT}\t${if (PLexist) proBand(PL) else -1}\t${phaseQual}\t${ances}\t${par}\t${kids}\t${desc}\t${exFamFreq}\t${popFreq}\t${popFreq.toFloat/(animalIDS.size)}\t${proRatio._2/(proRatio._1.toFloat + proRatio._2)}\t${rank}\tdenovo\t" + 
							proRatio + "\t" + selROvAD(par1,AD, RO, AO, GT) + " " + (if (PLexist) par1(PL) else "0,0,0") + "\t" + selROvAD(par2,AD, RO, AO, GT) + " " + (if (PLexist) par2(PL) else "0,0,0") +  s"\t${popRef}\t${popALT}")
							if (adratio != 0.0) {
								line(6) = "LOWQUAL_ADratio"
								print("\t WARNING: Low confidence de novo\t${allChildrenState}\n")
							}else {
								//statsOut.write("\n")
								print("\n")
							}//eelse
							out_vcf.write(line.reduceLeft{(a,b) => a + "\t" + b} + "\n")
						}//eif reoccur

					} //eif is Denovo
					
					
				/*	if(!valGTs.contains(proBand(GT)(0).toString + proBand(GT)(2)) && (ances == 0) && (par == 0) && (kids == 0) 
						&& (popFreq == 0) && checkDP(curPro, DP, minDP, maxDP) && checkDP(par1,DP,minDP,maxDP) && checkDP(par2,DP,minDP,maxDP)
						&& proRatio._2 >= minALT && adratio == 0.0){
							print(s"${line(0)}\t${line(1)}\t${line(3)}\t${line(3).size}\t${line(4)}\t${line(4).size}\t${line(5)}\t${fam._1}\t'${proGT}\t${if (PLexist) proBand(PL) else -1}\t\t${ances}\t${par}\t${kids}\t${desc}\t${exFamFreq}\t${popFreq}\t${popFreq.toFloat/(animalIDS.size)}\t${proRatio._2/proRatio._1.toFloat}\t${rank}\tSomatic\t" +
							 proRatio + "\t" + selROvAD(par1,AD, RO, AO, GT) + " " + (if (PLexist) par1(PL) else "0,0,0") + "\t" + selROvAD(par2,AD, RO, AO, GT) + " " + (if (PLexist) par2(PL) else "0,0,0") + s"\t${popRef}\t${popALT}\t\n")
							out_somatic.write(line.reduceLeft{(a,b) => a + "\t" + b} + "\n")
						}
						*/

				}//eisVAR
			} // IS VAR
			}//Efor fam <- trios
		}//IF Valid VCF line
	}// Ewhile
	in_vcf.close
	out_vcf.close
/*	
	for (child <- childHaps){
		val childOrigin = child._2.reverse.toArray
		errors.print(childOrigin.size + "\n")
		val out = new BufferedWriter(new FileWriter(child._1 + "-origin.bed"))
		var count = 1
		var lastPos, lastChr, lastOrig, startPos = "" 
		while (count < childOrigin.size){
			val cline = childOrigin(count).split("\t")
			if (lastChr == ""){
				lastChr = cline(0)
				lastPos = cline(1)
				lastOrig = cline(3)
				startPos = cline(1)
				out.write(s"${cline(0)}\t${cline(1)}\t")
			}
				if ((cline(3) != lastOrig) || (lastChr != cline(0))){
					out.write(s"${lastPos}\t${lastOrig}\n")
					if (phaseBlock.contains(child._1)){
						phaseBlock(child._1) =  (startPos.toInt,lastPos.toInt,lastChr,lastOrig) :: phaseBlock(child._1)
					} else {
						phaseBlock += child._1 -> (startPos.toInt,lastPos.toInt,lastChr,lastOrig) :: Nil
					}
					
					startPos = cline(1)
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
		out.close
	}
	*/
	//}
	
}//eMain

}//eObject
