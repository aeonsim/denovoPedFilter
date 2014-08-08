package pedFilter

/* Advanced Pedigree Filter
* Takes in VCF (gz) & Pedigree file, with cutoffs for DP, and a list of Proband ID's
* Use pedigree file to build family tree for each Proband, Proband, Children, Descendants, Parents, Ancestors
* Select variants based on pedigree & min DP check AD for low coverage
* Parents & Ancestors -ve, Proband +ve, Children 1+ positive, Descendants possibly +ve
* Output VCF per Proband ID of possible De novos along with, score file with counts in the population.
* Check to see if each VCF line is the correct length & discard if not.
*/

/* phaseTracker takes in phase blocks and tracks position within the block for each query*/
import scala.collection.mutable.HashMap

class phaseTracker (phaseData: HashMap[String,Array[Tuple5[String,Int,Int,String,Int]]]){
	
	var curPhaseBlocks = new HashMap[String,Int]
	for (c <- phaseData.keys){
		curPhaseBlocks += c -> 0
	}
	
	
	def getPhase(chrom: String, position: Int) : String ={
	
	//Simple block ID
		val block = phaseData(chrom).filter(pos => (pos._2 <= position && pos._3 >= position))
		if (block.size == 1){
			block(0)._4
		} else {
			"BOTH"
		}

		
	/*	if (phaseData(chrom).size == 0){
			//System.err.println("Empty Chrom " + chrom)
			return "NONE"
		} else {
			if (curPhaseBlocks(chrom) >= phaseData(chrom).size){
				return "END"
			} else {
				val curStart = phaseData(chrom)(curPhaseBlocks(chrom))._2
				val curEnd = phaseData(chrom)(curPhaseBlocks(chrom))._3
				if (position > curEnd){
					curPhaseBlocks(chrom) += 1
					return getPhase(chrom, position)
				} else {
					if (position >= curStart){
						return phaseData(chrom)(curPhaseBlocks(chrom))._4
					} else {
						return "BOTH"
					}
				}
			}
		}	*/
			
	} //ed def
	
	def getNearestBlock(chrom: String, position: Int) : String = {
		if (curPhaseBlocks(chrom) >= phaseData(chrom).size){
			return phaseData(chrom)(curPhaseBlocks(chrom) - 1)._4
		} else {
			if (curPhaseBlocks(chrom) >= 1 && phaseData(chrom).size != 0){
				val cur = scala.math.abs(position - phaseData(chrom)(curPhaseBlocks(chrom))._2)
				val prev = scala.math.abs(position - phaseData(chrom)(curPhaseBlocks(chrom) - 1)._2)
				if (prev <= cur) return phaseData(chrom)(curPhaseBlocks(chrom) - 1)._4 else return phaseData(chrom)(curPhaseBlocks(chrom))._4				
			} else {
				"NONE"
			}
		}
	
	}
}




object pedigreeFilter {

import java.io._
import net.sf.samtools.util._
import scala.collection.mutable.HashMap
import org.apache.commons.math3.distribution.BinomialDistribution

/* Global Var for VCF Type*/
var vcfType = ""
val errors = System.err
var PLexist = false
var PL = -1
var ped4gen: List[Tuple2[Array[Int],Int]] = Nil
var ped3_5dam: List[Tuple2[Array[Int],Int]] = Nil 
var ped3_5sire: List[Tuple2[Array[Int],Int]] = Nil 
for (p <- Array("RR","RA","AA"); s <- Array("RR","RA","AA"); d <- Array("RR","RA","AA"); sgs <- Array("RR","RA","AA"); sgd <- Array("RR","RA","AA"); dgs <- Array("RR","RA","AA"); dgd <- Array("RR","RA","AA")){
	val ped = Array(gts2Ints(p),gts2Ints(s),gts2Ints(sgs),gts2Ints(sgd),gts2Ints(d),gts2Ints(dgs),gts2Ints(dgd))
	ped4gen = (ped, (mendelian(p,s,d) + mendelian(s,sgs,sgd) + mendelian(d,dgs,dgd))) :: ped4gen
}

for (p <- Array("RR","RA","AA"); s <- Array("RR","RA","AA"); d <- Array("RR","RA","AA"); sgs <- Array("RR","RA","AA"); sgd <- Array("RR","RA","AA")){
	val ped = Array(gts2Ints(p),gts2Ints(s),gts2Ints(sgs),gts2Ints(sgd),gts2Ints(d),0,0)
	ped3_5sire = (ped, (mendelian(p,s,d) + mendelian(s,sgs,sgd))) :: ped3_5sire
}

for (p <- Array("RR","RA","AA"); s <- Array("RR","RA","AA"); d <- Array("RR","RA","AA"); dgs <- Array("RR","RA","AA"); dgd <- Array("RR","RA","AA")){
	val ped = Array(gts2Ints(p),gts2Ints(s),0,0,gts2Ints(d),gts2Ints(dgs),gts2Ints(dgd))
	ped3_5dam = (ped, (mendelian(p,s,d) + mendelian(d,dgs,dgd))) :: ped3_5dam
}


/* Return a probability from a Phred Score */

def phredProb(phred: Double) : Double = {
	if (vcfType != "gatk") {
		scala.math.pow(10.0,(-(scala.math.abs(phred))))
	} else {
		scala.math.pow(10.0,(-(scala.math.abs(phred)/10.0)))
	}
}

/* Extracts PL/GL score and returns as an Array of Doubles*/
def extctPL(geno: String): Array[Double] = {
  if(geno != "." && geno(0) != '.'){
	  val tmp = geno.split(":").apply(PL)
	  if (tmp != "."){
	   return tmp.split(",").map(_.toDouble)
	  } else {
	    return Array(0.0,0.0,0.0)
	  }
  }else {
    return Array(0.0,0.0,0.0)
  }
}

/* returns 0 if Mendelian 1 if not Mendelian */

def mendelian(child: String, sire: String, dam: String) : Int ={
var result: Set[String] = Set()
for (i <- sire.toList; j <- dam.toList) result = result + (i.toString + j.toString) + (j.toString + i.toString)
if (result.contains(child)) 0 else 1
}

/* Converts RR/RA/AA Genotypes to Alt Allele Counts for PL Index*/

def gts2Ints(gt: String): Int = {
if (gt == "RR"){
	return 0
} else {
	if (gt == "RA") return 1 else return 2
}
}

/* Estimates posterior probablities for a mendelian site, denovo site, and bad site (no inheritance)*/
def denovoPostProb4gen( proPL: Array[Double], sirePL: Array[Double], damPL: Array[Double], sgsPL: Array[Double], sgdPL: Array[Double], dgsPL: Array[Double], dgdPL: Array[Double],combos: List[Tuple2[Array[Int],Int]]) : Tuple3[Double,Double,Double] = {

val mendelianGTs = combos.filter(flt => flt._2 == 0)
val denovos = combos.filter(flt => flt._2 == 1)
val badGTs = combos.filter(flt => flt._2 > 1)

val L_R = (phredProb(proPL(0)) * phredProb(sirePL(0)) * phredProb(damPL(0)))
val P_denovo = 4E-8
val P_bad = 1E-4
val P_mendelian = 1E-3
val P_R = (1.0 - P_denovo - P_mendelian)
var L_mendelian, L_denovo, L_bad = 0.0

for (mend <- mendelianGTs){
	L_mendelian += (phredProb(proPL(mend._1(0))) * phredProb(sirePL(mend._1(1))) * phredProb(damPL(mend._1(4))) * phredProb(sgsPL(mend._1(2))) * phredProb(sgdPL(mend._1(3))) * phredProb(dgsPL(mend._1(5))) * phredProb(dgdPL(mend._1(6))) )
}
for (not <- denovos){
	L_denovo += (phredProb(proPL(not._1(0))) * phredProb(sirePL(not._1(1))) * phredProb(damPL(not._1(4))) * phredProb(sgsPL(not._1(2))) * phredProb(sgdPL(not._1(3))) * phredProb(dgsPL(not._1(5))) * phredProb(dgdPL(not._1(6))) )
}
for (bad <- badGTs){
	L_bad += (phredProb(proPL(bad._1(0))) * phredProb(sirePL(bad._1(1))) * phredProb(damPL(bad._1(4))) * phredProb(sgsPL(bad._1(2))) * phredProb(sgdPL(bad._1(3))) * phredProb(dgsPL(bad._1(5))) * phredProb(dgdPL(bad._1(6)))  )
}

val mendelian_posterior = ((L_mendelian * P_mendelian)/( (L_mendelian * P_mendelian) + (L_denovo * P_denovo) + (L_bad * P_bad) + (L_R * P_R)))
val denovo_posterior = ((L_denovo * P_denovo)/( (L_mendelian * P_mendelian) + (L_denovo * P_denovo) + (L_bad * P_bad) + (L_R * P_R)))
val bad_posterior = ((L_bad * P_bad)/( (L_mendelian * P_mendelian) + (L_denovo * P_denovo) + (L_bad * P_bad) + (L_R * P_R)))

(mendelian_posterior, denovo_posterior, bad_posterior )
}

/* Returns a Tuple containing the Posterior Prob of Mendelian & Denovo states */

def denovoPostProb(proPL: Array[Double], sirePL: Array[Double], damPL: Array[Double]) : Tuple2[Double,Double] = {
		val mendelian3G = Array(Array(2,2,2),Array(2,2,1),Array(2,1,2),Array(2,1,1),Array(1,2,1),Array(1,2,0),Array(1,1,2),Array(1,1,1),Array(1,1,0),Array(1,0,2),Array(1,0,1),Array(0,1,1),Array(0,1,0),Array(0,0,1),Array(0,0,0))
		val nonmend3G = Array(Array(0,0,2),Array(0,2,0),Array(0,2,2),Array(0,2,1),Array(0,1,2),Array(2,0,0),Array(2,0,2),Array(2,0,1),Array(2,2,0),Array(2,1,0),Array(1,0,0),Array(1,2,2))
		val L_R = (phredProb(proPL(0)) * phredProb(sirePL(0)) * phredProb(damPL(0)))
		val P_denovo = 4E-8
		val P_mendelian = 1E-3
		val P_R = (1.0 - P_denovo - P_mendelian)
		var L_mendelian, L_denovo = 0.0

		for (mend <- mendelian3G){
			L_mendelian += (phredProb(proPL(mend(0))) * phredProb(sirePL(mend(1))) * phredProb(damPL(mend(2))))
		}
		for (not <- nonmend3G){
			L_denovo += (phredProb(proPL(not(0))) * phredProb(sirePL(not(1))) * phredProb(damPL(not(2))))
		}

		val mendelian_posterior = ((L_mendelian * P_mendelian)/( (L_mendelian * P_mendelian) + (L_denovo * P_denovo) + (L_R * P_R)))
		val denovo_posterior = ((L_denovo * P_denovo)/( (L_mendelian * P_mendelian) + (L_denovo * P_denovo) + (L_R * P_R)))

		(mendelian_posterior, denovo_posterior)
}



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
		var result: Set[String] = Set()
		for (i <- str.toList; j <- str1.toList) result = result + (i.toString + j.toString) + (j.toString + i.toString)
		result.toList
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

/* Is this a Variant ie not 0/0 */

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
	
	
	/* Is this a Variant ie not ./. */

	def isValid(genotype:String): Boolean = {
		if (genotype(0) == '.') {
			false
		} else {
 			true
 		}
	}
	
	/* Return Genotypes as List*/
	
	def rtGTs(genotype: String) : List[String] = {
		val gts = if ( genotype.contains('/')) genotype.split('/') else genotype.split('|')
		gts.toList
	}
	
	/* Select method for returning Ref & Alt allele counts*/

	def selROvAD(indv: Array[String], ADval: Int, ROval: Int, AOval: Int, GTval: Int): Tuple2[Int,Int] = {
		var refAlt = (1,-1)
		if (indv(GTval)(0) != '.') vcfType match {
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
			//val GTlist = rtGTs(indv(GTval)).map(_.toInt).sorted
			//(RefAlt(GTlist(0)).toInt,RefAlt(GTlist(1)).toInt)
			(RefAlt(0).toInt,RefAlt(1).toInt)
		} else {
		if (PLexist){
			//Check other Genotypes are atleast Phred 10 unlikely
			if (indv(PL).split(",").sorted.tail.forall(num => scala.math.abs(num.toDouble) > 10)) (1,0) else (1,-1)
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
			if (scala.math.abs(pls(0).toDouble) >= min) true 
			else false
		} else {
			true
		}	
	}
	
	
	
/* Phase Code, return format is (sireAllele,damAllele)*/

def phase(indv: Array[String], sire: Array[String], dam: Array[String]) : Tuple2[String, String] = {
if (sire.size >= PL && dam.size >= PL && indv.size >= PL) {
	val minPLval = if (vcfType == "gatk") 80 else 5
	val indvGT = indv(0)
	val sireGT = sire(0)
	val damGT  = dam(0)
	val indvPL = indv(PL).split(",").map(num => scala.math.abs(num.toDouble)).sorted.tail 
	val sirePL = sire(PL).split(",").map(num => scala.math.abs(num.toDouble)).sorted.tail
	val damPL  = dam(PL).split(",").map(num => scala.math.abs(num.toDouble)).sorted.tail
if ((indvPL(0) >= minPLval && sirePL(0) >= minPLval && damPL(0) >= minPLval)){
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
	val minPLval = if (vcfType == "gatk") 75 else 5
	val childPL = child(PL).split(",").map(num => scala.math.abs(num.toDouble))
	if ((curPhase._1 != "x") && (childPL(1) >= minPLval) && (childGT == "1/1" || childGT == "0/0")){
		if (curPhase._1 == childGT(0).toString || curPhase._1 == childGT(2).toString) "S" else "D"
		} else {
		"U"
		}
	} else {
		"U"
	}
}

def findRecomb() : Unit = {
  
}

def findGConv() : Unit = {
  
}


def writeDepths(depths: HashMap[String, Array[Int]]): Unit = {
  	for (indv <- depths){
		val out = new BufferedWriter(new FileWriter(indv._1 + "-depths.txt"))
		if (indv._2.size > 5){
			val sum = indv._2.sum
			for (num <- 0 to 51){
				out.write(s"${num}\t${indv._2(num)}\t${indv._2(num)/sum.toFloat * 100}\n")
			}
		} else {
			out.write(s"Autozome good: ${indv._2(0)}\t${indv._2(0)/(indv._2(0) + indv._2(1) + 1.0) * 100}\n")
			out.write(s"Autozome bad : ${indv._2(1)}\t${indv._2(1)/(indv._2(0) + indv._2(1) + 1.0) * 100}\n")
			out.write(s"ChrX good    : ${indv._2(2)}\t${indv._2(2)/(indv._2(2) + indv._2(3) + 1.0) * 100}\n")
			out.write(s"ChrX bad     : ${indv._2(3)}\t${indv._2(3)/(indv._2(2) + indv._2(3) + 1.0) * 100}\n")
		}
		out.close
	}
}
	

/* Main body of Program */

def main (args: Array[String]): Unit = {

	var settings = new HashMap[String,String]
	
	for (items <- args){
		val keyVal = items.split("=")
		settings += keyVal(0).toUpperCase -> keyVal(1) 
	}
	

	if ((! settings.contains("VCF")) && (! settings.contains("PED")) && (! settings.contains("TRIOS"))  && (! settings.contains("TYPE"))) {
		println("advFilter VCF=input.vcf.gz PED=input.ped TRIOS=input_probands.txt phase=phase.vcf.gz type=gatk,plat,fb ref=genome.fasta {maxDPx=2.0 phaseVCF=SNPChip.vcf.gz minDP=0 minALT=0 RECUR=F/T minKIDS=1 minPL=0 QUAL=0 minRAFQ=0.2 }")
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
	val phaseVCF = if (settings.contains("PHASE")) settings("PHASE") else settings("VCF")
	val useRef = if (settings.contains("REF")) true else false
	val gpdb = if (settings.contains("GPDP")) true else false
	val maxDPmulti = if (settings.contains("MAXDPX")) settings("MAXDPX").toDouble else 2.0
	//val region = if (settings.contains("REGION")) true else false
	settings("TYPE").toUpperCase match {
		case "GATK" => vcfType = "gatk"
		case "PLAT" => vcfType = "platypus"
		case "FB" => vcfType = "freebayes"
		case _ => println("\n\nUnknown VCF type, VCF type must be on of gatk, plat (platypus) or fb (Freebayes)"); System.exit(1)
	}

	val outname = settings("VCF").split("/")(settings("VCF").split("/").size - 1)
	val out_vcf = new BufferedWriter(new OutputStreamWriter(new BlockCompressedOutputStream(outname + ".mutations-" + reoccur + "-denovos.vcf.gz")))
	val out_ped = new BufferedWriter(new FileWriter(outname + "-pedigree.tab"))
	//val out_somatic = new BufferedWriter(new OutputStreamWriter(new BlockCompressedOutputStream(outname + ".mutations-" + reoccur + "-somatic.vcf.gz")))
	//val statsOut = new BufferedWriter(new FileWriter(settings("OUT")))
	
	var pedFile = new HashMap[String, Array[String]]

	/* Tuple6(Ancestors, Parents, Children, Descendants, TrioDP, population, extended-Fam)*/

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
	var allChildren = new HashMap[String,String]

/*
* Iterate through VCF file to find Column headers and store positions in Map for later recall
*/

	var vcfrec = in_vcf.readLine().split("\t")
	while (vcfrec(0).apply(1) == '#'){
		out_vcf.write(vcfrec.reduceLeft{(a,b) => a + "\t" + b} + "\n")
		vcfrec = in_vcf.readLine().split("\t")
	}

	out_vcf.write(s"##${args.reduceLeft{(a,b) => a + " " + b}}\n")

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
		descendents = Nil
		extFam = Nil
		if (pedFile.contains(curPro(0)) && animalIDS.contains(curPro(0))){
			parents = List(pedFile(curPro(0))(2),pedFile(curPro(0))(3))
			ancestors = itPed(pedFile,curPro(0),2).tail.filterNot(x => parents.contains(x) || (! vcfanimals.contains(x)))
			children = findChildren(pedFile,vcfanimals,curPro(0))
			for (indv <- animalIDS.toList){
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
				allChildren += (child + "_" + curPro(0)) -> ""
				}
			}
		}//eif
	}//ewhile

	in_pro.close
	trios.foreach(tri => System.err.println(s"${tri._1} <- ${tri._2._2(0)} ${tri._2._2(1)} => Kids: ${tri._2._3}"))
	System.err.println(s"Built ${trios.size} Pedigrees\n")
	
	var otherTrios = new HashMap[String, Tuple3[String, String,List[String]]]
	for (anml <- animalIDS){
	  val cSire = pedFile(anml).apply(2)
	  val cDam = pedFile(anml).apply(3)
	  if (pedFile.contains(anml) && animalIDS.contains(cSire) && animalIDS.contains(cDam)){
	    var kids : List[String] = findChildren(pedFile,vcfanimals,anml)
	    if (kids.size != 0 && !trios.contains(anml)) otherTrios += anml -> (cSire, cDam, kids) //else System.err.println(s"${anml} <- ${cSire} ${cDam} : No Kids??")
	  }
	}
	//otherTrios = otherTrios.filterNot(data => trios.contains(data._1)) //.filter(data => data._2._3.size > 0)
	
	System.err.println(s"Built ${otherTrios.size} additional Trios\n")
	otherTrios.foreach(tri => System.err.println(s"${tri._1} <- ${tri._2._1} ${tri._2._2} => Kids: ${tri._2._3}"))

/* Report identified Trios & there Pedigree Structure*/

	for (fam <- trios){
		out_ped.write(s"TRIO:\t${fam._1}\t${fam._2._2(0)}\t${fam._2._2(1)}\n")
		out_ped.write("Grandparents\t" + fam._2._1.size + "\n")
		fam._2._1.foreach(s => out_ped.write(s + "\t\n"))
		out_ped.write("Children\t" + fam._2._3.size + "\n")		
		fam._2._3.foreach(s => out_ped.write(s + "\t\n"))
		out_ped.write("Descendents\t" + fam._2._4.size + "\n")
		fam._2._4.foreach(s =>out_ped.write(s + "\t\n"))
		out_ped.write("EXFAM\t" + fam._2._7.size + "\n")
		fam._2._7.foreach(s => out_ped.write(s + "\t\n"))
	}
	
	out_ped.close
	/* Build Reference Array, Note all positions are exact _ added as a Guard for 0 Done Async */
	
/*
* PrePhase Data using SNPChip data if available
*/

var phaseTracking = new HashMap[String, phaseTracker]
var readDepths = new HashMap[String,Array[Int]]
var phaseBlock = new HashMap[String,List[Tuple5[String,Int,Int,String,Int]]]

	val phaseInfo = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream(phaseVCF))))
	System.err.println("Have VCF, beginning Pre-phasing")
	var phaseLine = phaseInfo.readLine.split("\t")
	while (phaseLine(0).apply(1) == '#') phaseLine = phaseInfo.readLine.split("\t")
	for (col <- 9 to (phaseLine.size -1)){
		readDepths += phaseLine(col) -> new Array[Int](52)
	}
	
/* Will output RAW data to disk rather than trying to store in memory */
	
	var rawOutput = new HashMap[String,BufferedWriter]
	for (fam <- trios){
		/* Record Sites that are decent in Trio */
		readDepths += s"Trio_${fam._1}" -> new Array[Int](4)
		for (kid <- fam._2._3){
			rawOutput += (kid + "_" + fam._1) -> new BufferedWriter(new FileWriter(kid + "_" + fam._1 + "-rawPhase.txt"))
		} 
	}
	
	for (fam <- otherTrios){
		//readDepths += s"Trio_${fam._1}" -> new Array[Int](4)
		for (kd <- fam._2._3){
		  if (! rawOutput.contains(kd + "_" + fam._1) ){
			 rawOutput += (kd + "_" + fam._1) -> new BufferedWriter(new FileWriter(kd + "_" + fam._1 + "-rawPhase.txt"))
		  }
		}
	}
	
	var vcfChrs : List[String] = Nil
	//rawOutput.foreach(kd => System.err.println(kd._1))
	
	while(phaseInfo.ready){
	
		phaseLine = phaseInfo.readLine.split("\t")
		if (!(vcfChrs.contains(phaseLine(0)))) vcfChrs = phaseLine(0) :: vcfChrs
		
		val format = phaseLine(8).split(":")
		val GT = format.indexOf("GT")
		DP = if (format.contains("NV")) format.indexOf("NR") else format.indexOf("DP")
		if (format.contains("PL") || format.contains("GL")){
			PL = if (format.contains("PL")) format.indexOf("PL") else format.indexOf("GL")
			PLexist = true
		}
	  	if (phaseLine(3).size == 1 && phaseLine(4).size == 1 && phaseLine(5).toDouble > 1000){
		for (fam <- trios.par){
			/* Family = (ancestors, parents, children, tmpdesc, curPro(1).toInt, population, extFam) */
			val family = fam._2
			val maxDP = (family._5 * maxDPmulti).toInt
				val proband = phaseLine(vcfanimals(fam._1)).split(":")
				val sire = phaseLine(vcfanimals(family._2(0))).split(":")
				val dam = phaseLine(vcfanimals(family._2(1))).split(":")
				
				val pDP = if (proband.size > DP && DP != -1 && proband(DP) != ".") proband(DP).toInt else 0
				val sDP = if (sire.size > DP && DP != -1 && sire(DP) != "." ) sire(DP).toInt else 0
				val dDP = if (dam.size > DP && DP != -1 && dam(DP) != ".") dam(DP).toInt else 0
				
				if (pDP >= 51) readDepths(fam._1)(51) += 1 else readDepths(fam._1)(pDP) += 1
				if (sDP >= 51) readDepths(family._2(0))(51) += 1 else readDepths(family._2(0))(sDP) += 1
				if (dDP >= 51) readDepths(family._2(1))(51) += 1 else readDepths(family._2(1))(dDP) += 1
								
				/* Phase Trio at Site 0 = Good Autozome, 1 = bad Autozome, 2 = Good X, 3 = Bad X */
				val phaseVal = if (checkDP(proband, DP, minDP, maxDP) && checkDP(sire,DP,minDP,maxDP) && checkDP(dam,DP,minDP,maxDP)) {
						if (phaseLine(0).toUpperCase != "CHRX") readDepths(s"Trio_${fam._1}")(0) += 1 else readDepths(s"Trio_${fam._1}")(2) += 1 
						phase(proband, sire, dam) 
					} else { 
						if (phaseLine(0).toUpperCase != "CHRX") readDepths(s"Trio_${fam._1}")(1) += 1 else readDepths(s"Trio_${fam._1}")(3) += 1 
						("x","x")
					}
						
				for (kid <- family._3){
					var inherited = ""
					var curKid = phaseLine(vcfanimals(kid)).split(":")
					val kidDP = if (curKid.size > DP && DP != -1 && curKid(DP) != ".") curKid(DP).toInt else 0
					if (kidDP >= 51) readDepths(kid)(51) += 1 else readDepths(kid)(kidDP) += 1
					
					val sireid = pedFile(kid).apply(2)
					val damid = pedFile(kid).apply(3)
					if (vcfanimals.contains(sireid) && vcfanimals.contains(damid)) {
						val cSire = phaseLine(vcfanimals(sireid)).split(":")
						val cDam = phaseLine(vcfanimals(damid)).split(":")
						val fullPhase = if (checkDP(curKid, DP, minDP, maxDP) && checkDP(cSire,DP,minDP,maxDP) && checkDP(cDam,DP,minDP,maxDP)) phase(curKid, cSire, cDam) else ("x","x")
						if (fullPhase._1 != "x") { 
							if (fullPhase._1 == phaseVal._1) inherited = "S" else inherited = "D" 							
						} else {
							inherited = if (checkDP(curKid, DP, 6, maxDP)) childPhase(phaseVal,curKid) else "U"
						}
					} else {
						inherited = if (checkDP(curKid, DP, 6, maxDP)) childPhase(phaseVal,curKid) else "U"
					}
					
					if (inherited != "U") {
						val parID = if (inherited == "S") family._2(0) else family._2(1)
						rawOutput(kid + "_" + fam._1).write(s"${phaseLine(0)}\t${phaseLine(1)}\t${parID}\t${parID}\tS: ${sire(0)}\tD: ${dam(0)}\tP: ${proband(0)} \tK: ${curKid(0)} \t ${phaseLine(vcfanimals(family._2(0)))} ${phaseLine(vcfanimals(family._2(1)))} ${phaseLine(vcfanimals(fam._1))} ${phaseLine(vcfanimals(kid))}\t${family._2(0)} ${family._2(1)} ${fam._1}\n")
					}
				}
		}
		
		for (fam <- otherTrios.par){
		  
		  	val sireDamKids = fam._2
			val maxDP = 50
			val probandOther = phaseLine(vcfanimals(fam._1)).split(":")
			val sireOther = phaseLine(vcfanimals(sireDamKids._1)).split(":")
			val damOther = phaseLine(vcfanimals(sireDamKids._2)).split(":")
			
			val phaseVal = if (checkDP(probandOther, DP, minDP, maxDP) && checkDP(sireOther,DP,minDP,maxDP) && checkDP(damOther,DP,minDP,maxDP)) { 
						phase(probandOther, sireOther, damOther) 
					} else {
						("x","x")
					}
		  	
			for (kid <- sireDamKids._3){
				var curKid = phaseLine(vcfanimals(kid)).split(":")
				var inherited = childPhase(phaseVal,curKid)
				if (inherited != "U") {
					val parID = if (inherited == "S") sireDamKids._1 else sireDamKids._2
					rawOutput(kid + "_" + fam._1).write(s"${phaseLine(0)}\t${phaseLine(1)}\t${parID}\tS: ${sireOther(0)}\tD: ${damOther(0)}\tP: ${probandOther(0)} \tK: ${curKid(0)} \t ${phaseLine(vcfanimals(sireDamKids._1))} ${phaseLine(vcfanimals(sireDamKids._2))} ${phaseLine(vcfanimals(fam._1))} ${phaseLine(vcfanimals(kid))}\t${sireDamKids._1} ${sireDamKids._2} ${fam._1}\n")
				}
			}
			
		}
		
		
		}
	  	
	  	
	} //While Phasing
	
	for (fileout <- rawOutput){
		fileout._2.close
	}
	
	phaseInfo.close
	
	writeDepths(readDepths)
	
		for(child <- rawOutput){
			val in = new BufferedReader(new FileReader(child._1 + "-rawPhase.txt"))
			val out = new BufferedWriter(new FileWriter(child._1 + "-origin.bed"))
			var blocks = new HashMap[String,Array[Tuple5[String, Int, Int, String,Int]]]
			var phased : List[Tuple3[String,Int,String]] = Nil 
			
		while (in.ready){
				/* Read file into Array of tuples */
				val line = in.readLine.split("\t")
				if (line.size >= 3) phased = (line(0),line(1).toInt,line(2)) :: phased
			}
			in.close
			
			if(phased.size > 0) {
						
			val childOrigin = phased.reverse.toArray
			
			var count, score = 0			
			var chrom = childOrigin(count)._1
			var parent = childOrigin(count)._3
			var start = childOrigin(count)._2
			var end = childOrigin(count)._2
			
			if (!phaseBlock.contains(child._1) ) phaseBlock += child._1 -> Nil
			
			while (count < (childOrigin.size -1)){
				count +=1
				if (parent != childOrigin(count)._3 || chrom != childOrigin(count)._1){
					phaseBlock(child._1) = Tuple5(chrom,start,end,parent,score) :: phaseBlock(child._1)
					out.write(s"${chrom} ${start} ${end} ${parent} ${score}\n")
					chrom = childOrigin(count)._1
					parent = childOrigin(count)._3
					start = childOrigin(count)._2
					end = childOrigin(count)._2
					score = 1
				} else {
					end = childOrigin(count)._2
					score += 1
				}				
			}
			phaseBlock(child._1) = Tuple5(chrom,start,end,parent,score) :: phaseBlock(child._1)
			out.write(s"${chrom} ${start} ${end} ${parent} ${score}\n")
			out.close
			
			for (ch <- vcfChrs){
			  /*
				var goodblocks: List[Tuple5[String, Int, Int, String, Int]] = Nil
				val tmpblocks = phaseBlock(child._1).filter(blk => blk._5 > 2 && blk._1 == ch)
				var cntr = 0
				while(cntr < tmpblocks.size){
				  if (cntr + 1 < tmpblocks.size){
					  if (tmpblocks(cntr)._4 == tmpblocks(cntr + 1)._4){
					    goodblocks = (tmpblocks(cntr + 1)._1, tmpblocks(cntr + 1)._2, tmpblocks(cntr)._3, tmpblocks(cntr)._4, tmpblocks(cntr + 1)._5 + tmpblocks(cntr)._5) :: goodblocks
					  } else {
						  goodblocks = tmpblocks(cntr) :: goodblocks
					  }
				  }
				}
				*/
				blocks += ch -> phaseBlock(child._1).filter(blk => blk._5 > 2).filter(blk => blk._1 == ch).reverse.toArray
			}		
			phaseTracking += child._1 -> new phaseTracker(blocks)
		} else {
			for (ch <- vcfChrs){
				blocks += ch -> Array()
			}
			phaseTracking += child._1 -> new phaseTracker(blocks)
		}
	}

readDepths = new HashMap[String,Array[Int]]
phaseBlock = new HashMap[String,List[Tuple5[String,Int,Int,String,Int]]]

/* Load Ref into memory for Tri-Nucs*/

 //TEMP Disabled for speed reasons with testing
System.err.println("Reading Reference")

var refTable = new HashMap[String,Array[Char]]
if (useRef == true){

val ref = new BufferedReader(new FileReader(settings("REF")))
var line, refLastChr = ""
var refSeq = List('_')

while (ref.ready){
	line = ref.readLine
	if (line(0) == '>'){
		refTable += refLastChr -> refSeq.reverse.toArray
		refSeq = List('_')
		refLastChr = if (line.indexOf(' ') == -1) line.substring(1,line.size) else line.substring(1).split("\\s+")(0)
	} else {
		refSeq = line.toList.reverse ::: refSeq
	}
}
refSeq = Nil
ref.close
}

/*
*	Iterate through VCF file line by line, at each line load each Trio and count existence of variants in different categories
*	if de novo, flag and output snp detail and variant info, (count in pop, children ancestors etc)
*/
	//statsOut.write(s"Chrom\tPos\tRef\tRefSize\tAlt\tAltSize\tQUAL\tTrio\tGenotype\tPLs\tPhase\t Vars S|D Haps S|D\tAnces\tPars\tChildren\tDesc\tExFam\tPop\tPopFreq\tSupport Ratio\tScore\tClass\tProband\tSire\tDam\tPopRefCount\tPopAltCount\tWarning\tPhaseInfo\n")
	println(s"Chrom\tPos\trsID\tTRI-NUC\tRef\tAlt\tQUAL\tTrio\tGenotype\tPLs\tWorst Parent PL\tDenovo Posterior\t Probs (mend, denovo, bad) \tSource\tPhase\tVars S|D Haps S|D\tAnces\tPars\tChildren\tDesc\tExFam\tPop\tPopFreq\tSupport Ratio\tOffspring Support Ratio\tClass\tProband\tSire\tgSire\tgDam\tDam\tgSire\tgDam\tPopRefCount\tPopAltCount\tPhaseInfo\tRefSize\tAltSize\t")
	var lastChr = ""
	
	System.err.println("Phasing Complete, beginning De Novo identification & Characterisation")

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
			errors.print(line(0) + " ")
		}
		
		AD = format.indexOf("AD")
		GT = format.indexOf("GT")
		DP = if (format.contains("NV")) format.indexOf("NR") else format.indexOf("DP")
		AO = if (format.contains("NV")) format.indexOf("NV") else format.indexOf("AO")
		RO = if (format.contains("NV")) format.indexOf("NR") else format.indexOf("RO")
				
		/*To be considered the VCF record must be ok, the Qual score >= Min & no more than 3 alternative alleles*/
		if (line.size == (vcfanimals.size + 9) && (line(5).toFloat >= QUAL) && (line(4).split(",").size < 3)){
			var trioPos = 0
			val triosArray = trios.keys.toArray
			
			for (fam <- trios.toArray){
				var indv = ""
				var kidsPhase : List[String] = Nil
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

				val sDetails = line(vcfanimals(sire)).split(":")
				val dDetails = line(vcfanimals(dam)).split(":")
				val proBand = line(vcfanimals(fam._1)).split(":")
				val sGT = sDetails(GT)(0).toString + sDetails(GT)(2)
				val dGT = dDetails(GT)(0).toString + dDetails(GT)(2)
				val pGT = proBand(GT)(0).toString + proBand(GT)(2)
			
			val valGTs = permu(sGT,dGT)
				
	/*	Begin Denovo Logic 	*/
				
		//Is the site a Variant?
		if (isValid(proBand(GT)) && isValid(sDetails(GT)) && isValid(dDetails(GT)) && checkDP(proBand,DP,minDP,maxDP) && checkDP(sDetails,DP,minDP,maxDP) && checkDP(dDetails,DP,minDP,maxDP)){
			var allowDam, allowSire = false
			var sgSire, sgDam, dgSire, dgDam, tDam, tSire, gtgSire, gtgDam, dgSireAD, dgDamAD, sgDamAD, sgSireAD = ""
			var denovoParent = ""
			var sgsDetails, sgdDetails, dgsDetails, dgdDetails: Array[String] = Array()
			var gpdepth = true
			var gpAB = 0.0
			
			for (parent <- ped._2){
				var numGPs, gpVars = 0
				val parentDetails = line(vcfanimals(parent)).split(":")
				
				//GrandSire
				if (vcfanimals.contains(pedFile(parent)(2))){
					val gSire = line(vcfanimals(pedFile(parent)(2))).split(":")
					var tmpAD = selROvAD(gSire, AD,RO,AO,GT)
					if (tmpAD._2 != -1) gpAB += (tmpAD._2 / (tmpAD._1 + tmpAD._2.toFloat))
					tSire = tmpAD.toString
					if (checkDP(gSire,DP,minDP,maxDP)) {
						if (isVar(gSire(GT))) {
						  gpVars += 1
						} 
						gtgSire = s"${gSire(GT)(0)}${gSire(GT)(2)}"
						numGPs += 1
					} else {
						  gpdepth = false
						}
				}
					
				//GrandDam
				if (vcfanimals.contains(pedFile(parent)(3))){
					val gDam = line(vcfanimals(pedFile(parent)(3))).split(":")
					var tmpAD = selROvAD(gDam, AD,RO,AO,GT)
					if (tmpAD._2 != -1) gpAB += (tmpAD._2 / (tmpAD._1 + tmpAD._2.toFloat))
					tDam = tmpAD.toString
					if (checkDP(gDam,DP,minDP,maxDP)){
						if (isVar(gDam(GT))) {
						  gpVars += 1
						} 
						gtgDam = s"${gDam(GT)(0)}${gDam(GT)(2)}"
						numGPs += 1
					} else {
						  gpdepth = false
						}
				}
				
				if (isVar(parentDetails(GT))) {
					par += 1
					denovoParent = denovoParent + parent 
				}
				
				if (numGPs == 2) {
					if (parent == sire){
						allowSire = true
						sgSireAD = tSire
						sgDamAD = tDam
						sgSire = pedFile(parent)(2)
						sgDam = pedFile(parent)(3)
					} 
					if (parent == dam)	{
						allowDam = true
						dgSireAD = tSire
						dgDamAD = tDam
						dgSire = pedFile(parent)(2)
						dgDam = pedFile(parent)(3)
					  }
				} else {
				  if (parent == dam) {
					dgSireAD = tSire
					dgDamAD = tDam
					dgSire = pedFile(parent)(2)
					dgDam = pedFile(parent)(3)
				  } else {
				    sgSireAD = tSire
					sgDamAD = tDam
					sgSire = pedFile(parent)(2)
					sgDam = pedFile(parent)(3)
				  } 
				  
				}
				ances += gpVars
				//if (permu(gtgSire,gtgDam).contains(pGT)) ances += gpVars	
			}
			
		//println("Parents " + par + "Ancestors " + ances + " Kids " + kids + " Desc " + desc + " Pop " + popFreq)
			
		if (!(valGTs.contains(pGT)) || (ances == 0 && par == 1 && ((allowSire && isVar(sDetails(GT)) && !(isVar(dDetails(GT))) ) || (allowDam && isVar(dDetails(GT)) && !(isVar(sDetails(GT))) )))) {

			/* Loop through each family group and record Hets */
			
			var kidAlts, kidRefs = 0
			kids = 0
			popFreq = 0
			//Children
			var varChildren : List[String] = Nil
			for (indv <- ped._3){
						val curAn = line(vcfanimals(indv)).split(":")
						//if(curAn(0).apply(0) != '.' && curAn.size >= PL){
						if(isValid(curAn(GT))){
							val refAlt = selROvAD(curAn,AD, RO, AO, GT)		
							//System.err.println(indv + " " + curAn(AD) + " " + curAn(GT) + " isVar? " + isVar(curAn(GT)) + " refALT: " + refAlt)
							if (isVar(curAn(GT)) || refAlt._2 >= 1){
								kidRefs += refAlt._1
								kidAlts += refAlt._2
								kids += 1
								varChildren = indv :: varChildren
							}
						}
					}

			//Desec
				while(descPos < ped._4.size){
					indv = ped._4.apply(descPos)
					descPos += 1
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
						if (line(vcfanimals(indv))(0) != '.'){
							val curAn = line(vcfanimals(indv)).split(":")
							val refAlt = selROvAD(curAn,AD, RO, AO, GT)
							//System.err.println(indv + " " + curAn(AD) + " " + curAn(GT) + " isVar? " + isVar(curAn(GT)) + " refALT: " + refAlt)
							if ((isVar(curAn(GT)) || sigAD(refAlt._2)) && checkDP(curAn,DP,minDP,maxDP)){
								popRef += refAlt._1
								popALT += refAlt._2
								val curIndv = pedFile(indv)
								if (vcfanimals.contains(curIndv(2)) && vcfanimals.contains(curIndv(3))){
								val psire = line(vcfanimals(curIndv(2))).split(":") 
								val pdam = line(vcfanimals(curIndv(3))).split(":")
								if (psire(0).apply(0) != '.' && pdam(0).apply(0) != '.'){
									if (isVar(psire(GT)) || isVar(pdam(GT))) {
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
						if (line(vcfanimals(indv))(0) != '.'){
							val curAn = line(vcfanimals(indv)).split(":")
							val refAlt = selROvAD(curAn,AD, RO, AO, GT)
							if ((isVar(curAn(GT)) || sigAD(refAlt._2)) && checkDP(curAn,DP,minDP,maxDP)){
								exFamFreq += 1
							}
						}
					}
					
				
					val debug = false				
					if (debug) {
					proBand.foreach(s => print(s"${s} "))
					 print(s"${line(1)} ${sGT} ${dGT} ${pGT} \n")
					 print(s"GT ${!valGTs.contains(proBand(GT)(0).toString + proBand(GT)(2))} ${valGTs} Pro:${proBand(GT)} PL:${checkPL(minPL, proBand)} SPL:${checkPL(minPL, sDetails)} DPL:${checkPL(minPL, dDetails)}\n")
					 print(s"DP:${checkDP(proBand, DP, minDP, maxDP)}  SDP:${checkDP(sDetails, DP, minDP, maxDP)} DDP:${checkDP(dDetails, DP, minDP, maxDP)} AN:${ances == 0} PAR:${par == 0} KID:${kids >= minKids} pop:${popFreq} \n")
					 println(s"ANS: ${ances} PAR: ${par} Kids: ${kids} Pop: ${popFreq}")
					 }

			/* -------- Denovo Check ---------- */

					/*
					var rank = ""
					altsPar match {
						case 0 => rank = "probable"
						case 1 => rank = "possible"
						case _ => rank = "unlikely"
					}
					*/
					
					
					val proRatio = selROvAD(proBand,AD, RO, AO, GT)
					val proGT = proBand(GT)
					
					val supportRatio = proRatio._2/(proRatio._2 + proRatio._1.toFloat)
					val offsupRatio = kidAlts / (kidAlts + kidRefs.toFloat)
					val popSupRatio = popALT / (popALT + popRef.toFloat)
					
					/*
					* De novo Identification logic!
					* Proband is Variant, has > min num of Alt alleles. DP is between minDP & maxDP
					* Parental DP is good, variant not present in Ancestors, Parents & is present in at least
					* minKids kids & min ratio between Ref & Alt is fine
					* If recurring allowed don't check population otherwise not present in population
					*/
					val genderCheck = if (Array("X","CHRX","CHX").contains(line(0).toUpperCase)){
						/* Check Gender from PedFile, 1 = Male, 2 = Female */
					  if (pedFile(fam._1).apply(4) == "1" && isHet(proGT)) false else true
					} else {
					  true
					}
					
					if ((ances == 0) && (kids >= minKids) && (par <= 1) && (gpAB <= minRAFreq) && ( supportRatio >= minRAFreq || offsupRatio >= minRAFreq) && ( if (reoccur && popFreq != 0) ((popSupRatio >= minRAFreq && gpdepth)) else popFreq == 0 ) && genderCheck) {
						
						var varSirePhase, varDamPhase = 0.0
						for (child <- ped._3){
						  val curAn = line(vcfanimals(child)).split(":")
							val childID = (child + "_" + fam._1)
							allChildren(childID) = phaseTracking(childID).getPhase(line(0),line(1).toInt)
							allChildren(childID) match {
								case `sire` => {sirePhase += 1; if (varChildren.contains(child)) varSirePhase += 1}
								case `dam` => {damPhase += 1 ; if (varChildren.contains(child)) varDamPhase += 1 }
								case _ =>
							}
							allChildrenState = allChildrenState + s"${child}:${if (curAn.size >= PL ) curAn(PL) else 0}:${curAn(GT)}:${allChildren(childID)}\t"		
						}	
	
						var phaseQual = ""
						val denovoHap = varSirePhase + "|" + varDamPhase + " " + sirePhase + "|" + damPhase
						var sourceHap = ""
						  
						if (varSirePhase >= 1 && varDamPhase >= 1){
							phaseQual = "Bad\t" + denovoHap
						} else {
						  if ((varSirePhase == kids && sirePhase == kids) || (varDamPhase == kids && damPhase == kids)){
						  //if ( (varSirePhase != 0 &&  varDamPhase != 0) && ((varSirePhase != 0 && varSirePhase == sirePhase) || (varDamPhase != 0 && varDamPhase == damPhase))){
								phaseQual = "Good\t" + denovoHap
							} else {
								phaseQual = "Partital\t" + denovoHap
							}
						  	sourceHap = if (varSirePhase != 0) "Pat" else "Mat"
						  	  
						  	  /*		
						  	   * 	Find GrandParent Phase for 4 Gen
						  	  */
		  	  
						  	  def grandPar(parnt: String, gSire: String, gDam: String): Unit = {
						  	  	val childID = (fam._1 + "_" + parnt)
						  	    if (phaseTracking.contains(childID)) phaseTracking(childID).getPhase(line(0),line(1).toInt) match{
						  	      case `gDam` => sourceHap = sourceHap + " grndDam"
						  	      case `gSire` => sourceHap = sourceHap + " grndSire"
						  	      case _ => //System.err.println(childID + " Phase of Parent " + phaseTracking(childID).getPhase(line(0),line(1).toInt) + " Sire: Dam:" + gSire + " " + gDam)
						  	    }
						  	  
						  	}
						  	  
						  	if (par == 1){
						  	 if (denovoParent == sire ) grandPar(sire, sgSire, sgDam)
						  	 if (denovoParent == dam) grandPar(dam, dgSire, dgDam)
						  	}
						}
						
						
						var varClass = ""
						  if ((isVar(proBand(GT)) && isVar(sDetails(GT)) && !isVar(dDetails(GT))) || (isVar(proBand(GT)) && isVar(dDetails(GT)) && !isVar(sDetails(GT)))) varClass = "gpDenovo"
						  if (!(valGTs.contains(pGT))) varClass = "denovo"
						    
						    /* Calculating Probablities*/
						    
						    if (varClass == "denovo"){
						      //Simple Case Denovo in Proband, is it Mosaic in a parent
						      var probs = Array(0,0,0,0,0,0,0)
						    }
						    
							var triNuc = if (useRef != true) "NNN" else refTable(line(0)).apply(line(1).toInt -1).toString + refTable(line(0)).apply(line(1).toInt).toString + refTable(line(0)).apply(line(1).toInt +1).toString
							var triNucAlt = if (useRef != true) "NNN" else triNuc(0) + line(4) + triNuc(2)

							var posteriorProbs = "No PL/GL"
							  if (PLexist){
							    val proPL = proBand(PL).split(",").map(_.toDouble)
							    val sirePL = sDetails(PL).split(",").map(_.toDouble)
							    val damPL = dDetails(PL).split(",").map(_.toDouble)
								  if (ped._1.size == 4){
									val sgsPL = extctPL(line(vcfanimals(pedFile(ped._2(0))(2))))
									val sgdPL = extctPL(line(vcfanimals(pedFile(ped._2(0))(3))))
									val dgsPL = extctPL(line(vcfanimals(pedFile(ped._2(0))(2))))
									val dgdPL = extctPL(line(vcfanimals(pedFile(ped._2(1))(3))))
									val proDenovo = denovoPostProb(proPL, sirePL, damPL)
									val damDenovo = denovoPostProb(damPL, dgsPL, dgdPL)
									val sireDenovo = denovoPostProb(sirePL, sgsPL, sgdPL)
								    val tmpError = denovoPostProb4gen(proPL, sirePL, damPL, sgsPL, sgdPL, dgsPL, dgdPL,ped4gen)._3
								    val mostLikelyDenovo = Array(proDenovo,sireDenovo,damDenovo).sortBy(_._2).apply(2)
								    posteriorProbs = s"${mostLikelyDenovo._2}\t${mostLikelyDenovo._1.toFloat} ${mostLikelyDenovo._2.toFloat} ${tmpError.toFloat}"
								  } else{
									  val arrayNull = Array(0.0,0.0,0.0)
									  if (ped._1.size >= 2 && (allowSire || allowDam)){
									    val tmp = if (allowSire){
									    val sgsPL = extctPL(line(vcfanimals(pedFile(ped._2(0))(2))))
									    val sgdPL = extctPL(line(vcfanimals(pedFile(ped._2(0))(3))))
									    val proDenovo = denovoPostProb(proPL, sirePL, damPL)
									    val sireDenovo = denovoPostProb(sirePL, sgsPL, sgdPL)
									    val mostLikelyDenovo = Array(proDenovo,sireDenovo).sortBy(_._2).apply(1)
									    val tmpError =  denovoPostProb4gen(proPL, sirePL, damPL, sgsPL, sgdPL, arrayNull, arrayNull, ped3_5sire)._3
									    posteriorProbs = s"${mostLikelyDenovo._2}\t${mostLikelyDenovo._1.toFloat} ${mostLikelyDenovo._2.toFloat} ${tmpError.toFloat}"  
									    } else {
									    val dgsPL = extctPL(line(vcfanimals(pedFile(ped._2(1))(2))))
									    val dgdPL = extctPL(line(vcfanimals(pedFile(ped._2(1))(3))))
									    val proDenovo = denovoPostProb(proPL, sirePL, damPL)
									    val damDenovo = denovoPostProb(damPL, dgsPL, dgdPL)
									    val mostLikelyDenovo = Array(proDenovo,damDenovo).sortBy(_._2).apply(1)
									    val tmpError = denovoPostProb4gen(proPL, sirePL, damPL, arrayNull, arrayNull, dgsPL, dgdPL, ped3_5dam)._3
									    posteriorProbs = s"${mostLikelyDenovo._2}\t${mostLikelyDenovo._1.toFloat} ${mostLikelyDenovo._2.toFloat} ${tmpError.toFloat}"
									    }
									    
									  }else {
										  val tmp = denovoPostProb(proPL, sirePL, damPL)
										  posteriorProbs = s"${tmp._2.toFloat}\t${tmp._1.toFloat} ${tmp._2.toFloat}"
									  }

								  }
							  }
				
						val pls = if (PLexist) (sDetails(PL) + "," + dDetails(PL)).split(",") else Array(".",".",".")
						var worstParent = if (pls.contains(".")) -1 else pls.map(_.toDouble).sorted.apply(2)
						
							print(s"${line(0)}\t${line(1)}\t${line(2)}\t${triNuc}>${triNucAlt}\t${line(3)}\t${line(4)}\t${line(5)}\t${fam._1}\t'${proGT}'\t${if (PLexist) proBand(PL) else -1}\t${worstParent}\t${posteriorProbs}\t${sourceHap}\t${phaseQual}\t${ances}\t${par}\t${kids}\t${desc}\t${exFamFreq}\t${popFreq}\t${popFreq.toFloat/(animalIDS.size)}\t${supportRatio}\t${offsupRatio}\t${varClass}\t" + 
								proRatio + "\t" + selROvAD(sDetails,AD, RO, AO, GT) + " " + (if (PLexist) sDetails(PL) else "0,0,0") + s"\t${sgSireAD}\t${sgDamAD}\t" + selROvAD(dDetails,AD, RO, AO, GT) + " " + (if (PLexist) dDetails(PL) else "0,0,0") + s"\t${dgSireAD}\t${dgDamAD}\t${popRef}\t${popALT} ${popSupRatio}\t${allChildrenState}\t${line(3).size}\t${line(4).size}\n")
							out_vcf.write(line.reduceLeft{(a,b) => a + "\t" + b} + "\n")

					} //eif is Denovo

				}//eisVAR
			//} // IS VAR
			} // Parent's can't make Child GT or Grandparents have no Var and one parent Does
		}//Efor fam <- trios
	}//IF Valid VCF line
}// Ewhile
	in_vcf.close
	out_vcf.close
	System.err.println("\ndone")
}//eMain

}//eObject
