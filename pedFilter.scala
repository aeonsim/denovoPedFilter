/*
*	pedFilter an advanced filtering tools for candidate de novos in population VCFs
*		Supplied with a list of filters & family information the tools will parse a
*		population VCF outputing variants which meet the provide criteria.
*/

import net.sf.samtools.util._
import java.io._
import scala.collection.mutable.HashSet
import scala.collection.mutable.HashMap

/*
* TODO
*	Performance
*	Cleaner support for missing atributes
*	Clean up and reformat code
*	Copy to git-hub
*	Version
*	Clean up CLI
*	Add max errors to CLI
*	Add flags for Errors
*/


object pedFilter{

/*
* Setup Family and Animal ID lists
*/

var family = new HashMap[String,String]

var famVar, popVar, proVar, childVar = 0

var validVars = Array("0/1","1/0","0|1","1|0","1/2","2/1","1|2","2|1","1/1","1|1","2/2","2|2")

var maxError = 0
var maxDP, maxGQ = 1000000
var minDP, minGQ, QUAL = 0
var minAD = Array(0,0)
var minPL = Array(0,0,0)
var ADerror, GQerror, DPerror = 0

var AD,DP,GQ,PL,GT = -1


def main(args: Array[String]): Unit = {
println("pedFilter input.vcf.gz family.txt GEmax0 QLmin50 DPmin4 DPmax60 GQminX ADminX,X PLminX,X,X \n GE Genotype errors allowed for low coverage")

val input = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream(args(0)))))
val inFam = new BufferedReader(new FileReader(args(2)))

/*
*	Filter Setup
*/

for (i <- args){
if (i.contains("vcf") || i.contains("txt")) {
//input = new BufferedReader(new FileReader(args(0)))
} else {
i.substring(0,5).toUpperCase match {
case "QLMIN" => QUAL = spliter(i)
case "DPMIN" => minDP = spliter(i)
case "DPMAX" => maxDP = spliter(i)
case "GQMIN" => minGQ = spliter(i)
case "GEMAX" => maxError = spliter(i)
case "ADMIN" => minAD = spliterList(i)
case "PLMIN" => minPL = spliterList(i)
case _ => println("ERROR, unknown filter field" + i);System.exit(1)
}//E Match
}//E if/else
}

val out = new BufferedWriter(new OutputStreamWriter(new BlockCompressedOutputStream(args(0) + s"GE${maxError}GQ${minGQ}DP${minDP}-${maxDP}AD${minAD(0)}.${minAD(1)}.output.vcf.gz")))

//Family.txt 0 is Parent/Grand 1 is proband 2 is child

while (inFam.ready){
var tmpLine = inFam.readLine.split("\t")
family += tmpLine(0) -> tmpLine(1)
}
inFam.close

var current = input.readLine.split("\t")
// Skips through VCF to header line
while (current(0).apply(1) == '#'){
out.write(current.reduceLeft{(a,b) => a + "\t" + b} + "\n")
current = input.readLine.split("\t")
}

//Store all Animal names
val animals = current
//out.write(animals.reduceLeft{(a,b) => a + "\t" + b} + "\n")
//val animals = current.toList.slice(9,cline.size).toArray
out.write(s"${current(0)}\t${current(1)}\t${current(2)}\t${current(3)}\t${current(4)}\t${current(5)}\t")
out.write(s"${current(6)}\t${current(7)}\t${current(8)}")
for (i <- family.keys){
if (animals.indexOf(i) != -1) out.write(s"\t${current(animals.indexOf(i))}")
}
out.write("\n")



/*
*	Start processing
*/

var pos = 0
while (input.ready){
//current = input.readLine.split("\t").par.splitAt(8)

current = input.readLine.split("\t") 
pos = current.size -1
//Apply Qual Filter
if (current(5).toDouble >= QUAL ){
//println(current(5))
//println(current.reduceLeft{(a,b) => a + "\t" + b} + "\n")
setFields(current(8))
//println(s"${GT}\t${DP}\t${GQ}\t${AD}")

while (pos > 8){
	var genotype = current(pos).split(":")
	//println("WTF\t" + animals(pos) + family.contains(animals(pos)))
	if (family.contains(animals(pos))){
	//println("Family Var")
		if (family(animals(pos)) == "0"){
			if (isVar(genotype(GT))){
				//println("fam" + s"${animals(pos)}" + genotype(GT) )
				famVar += 1
			}// Eif Variant
		} else { // Eif Parent
			if (family(animals(pos)) == "1") {
				if (isVar(genotype(GT))) {proVar += 1}
				//println("Pro" + s"${animals(pos)}" + genotype(GT) )
				}else {
				if (isVar(genotype(GT))) {childVar += 1}
				//println("Child" + s"${animals(pos)}" + genotype(GT) )
			}
		} // E Else Parent 
	} else {//E if Family
		if (passFilter(genotype)){
			popVar += 1
		}
	} // E else Family

pos -= 1
} //E While Pos
//println( s"${popVar}\t${famVar}\t${proVar}\t${childVar}") 
//println( s"error $ADerror $GQerror $DPerror")
if (popVar <= maxError && famVar == 0 && proVar == 1 && childVar >= 1){
//println("Writing")
out.write(s"${current(0)}\t${current(1)}\t${current(2)}\t${current(3)}\t${current(4)}\t${current(5)}\t")
out.write(s"${current(6)}\t${current(7)}\t${current(8)}")
for (i <- family.keys){
if (animals.indexOf(i) != -1) out.write(s"\t${current(animals.indexOf(i))}")
}
out.write("\n")
}
famVar = 0
popVar = 0
proVar = 0
childVar = 0
 ADerror = 0 
 GQerror = 0
 DPerror = 0
} // E IF QUAL 
}// E while Ready
out.close
}//E main


def spliter (field: String): Int = {
	if (field.contains("min")) field.split("min").apply(1).toInt else field.split("max").apply(1).toInt
}


def spliterList (field: String): Array[Int] = {
	if (field.contains("min")) field.split("min").apply(1).split(",").map(_.toInt) else field.split("max").apply(1).split(",").map(_.toInt)
}


def passFilter(genos: Array[String]): Boolean = {
	if (genos.size >= 2){
	//println(genos.reduceLeft{(a,b) => a + "\t" + b} + s"\t${genos.size}" +"\n")
		if (genos(GT).size == 3 && genos(GT).apply(0) ==  '0' && genos(GT).apply(2) ==  '0'){
			//Is Reference done!
			false
		} else {
			//Is variant need to filter
			if ((DP != -1 && genos(DP).toInt >= minDP && genos(DP).toInt <= maxDP) || (DP == -1)){
				if ((GQ != -1 && genos(GQ).toFloat >= minGQ )||(GQ == -1)){
					if ((AD != -1 && genos(AD)(0).toInt >= minAD(0) && genos(AD)(1).toInt >= minAD(1)) || (AD == -1)){
						true
					} else {
						ADerror += 1
						false
					} //e AD
				} else {
					GQerror += 1
					false
				}//e GQ
			} else {
				DPerror += 1
				false
			}//e DP

		} // e else
	} else false //eif size
} //e def


def setFields(geno: String): Unit = {
	AD = geno.split(":").indexOf("AD")
	DP = geno.split(":").indexOf("DP")
	GQ = geno.split(":").indexOf("GQ")
	PL = geno.split(":").indexOf("PL")
	GT = geno.split(":").indexOf("GT")
}

def isVar(geno: String): Boolean = {
	validVars.contains(geno)
}

}// E Object