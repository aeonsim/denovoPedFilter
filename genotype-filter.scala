import java.io._

object filterGenotypes {

def main (args:Array[String]): Unit = {

if (args.size == 0) {
println("scala filter file.vcf Genotype.DPmin Genotype.minGQ Genotype.minAD1,2 Genotype.maxPL-RR1,RA2,AA3\ncan supply min or Max, ie DPminDP4 DPmax20")
System.exit(1)
}

var input: BufferedReader = null
var maxDP, maxGQ = 1000000
var minDP, minGQ = 0
var minAD = Array(0,0)
var minPL = Array(0,0,0)
var AD,DP,GQ,PL = -1

def spliter (field: String): Int = {
if (field.contains("min")) field.split("min").apply(1).toInt else field.split("max").apply(1).toInt
}

def spliterList (field: String): Array[Int] = {
if (field.contains("min")) field.split("min").apply(1).split(",").map(_.toInt) else field.split("max").apply(1).split(",").map(_.toInt)
}


for (i <- args){
if (i.contains("vcf")) {
input = new BufferedReader(new FileReader(args(0)))
} else {
i.substring(0,5).toUpperCase match {
case "DPMIN" => minDP = spliter(i)
case "DPMAX" => maxDP = spliter(i)
case "GQMIN" => minGQ = spliter(i)
case "GQMAX" => maxGQ = spliter(i)
case "ADMIN" => minAD = spliterList(i)
case "PLMIN" => minPL = spliterList(i)
case _ => println("ERROR, unknown filter field");System.exit(1)
}//E Match
}//E if/else
}// E for
//println(s"maxDP $maxDP minGQ $minGQ")
var line = input.readLine.split("\t")

System.err.println(s"maxDP $maxDP minDP $minDP minGQ $minGQ maxGQ $maxGQ minAD ${minAD(1)}")

while(line(0)(0) == '#'){
println(line.reduceLeft{(a,b) => a + "\t" + b})
line = input.readLine.split("\t")
}

def setFields(geno: String): Unit = {
AD = geno.split(":").indexOf("AD")
DP = geno.split(":").indexOf("DP")
GQ = geno.split(":").indexOf("GQ")
PL = geno.split(":").indexOf("PL")
}
var DPfalse = false
if (DP == -1) DPfalse = true

while (input.ready){
print(s"${line(0)}\t${line(1)}\t${line(2)}\t${line(3)}\t${line(4)}\t${line(5)}\t${line(6)}\t${line(7)}\t${line(8)}")
setFields(line(8))
for (i <- (9 to (line.size -1 ))){
//println(s"Animal\t$i")
val geno = line(i).split(":")
if(geno.size < 4){
print("\t" + geno.reduceLeft{(a,b)=> a + ":" + b})
} else {
//System.err.println(s"${geno(DP)}\t${geno(GQ)}\t${geno(AD)}\t${geno(PL)}")
val genoAD = geno(AD).split(",")
if(
//(DPfalse ||(geno(DP) == ".") || ((geno(DP).toInt <= maxDP) && (geno(DP).toInt >= minDP))) && 
((geno(GQ) == ".") || ((geno(GQ).toInt <= maxGQ) && (geno(GQ).toInt >= minGQ))) &&
((geno(AD) == ".,." || geno(AD) == "." || geno(0) != "0/1") || ((genoAD(0).toInt >= minAD(0)) && (genoAD(1).toInt >= minAD(1)) && ((genoAD(1).toFloat / (genoAD(0).toFloat + genoAD(1).toFloat)) > 0.25)) &&
((geno(PL) == ".,.,." || geno(PL) == "." || geno(0) != "0/1") || ((geno(PL).split(",")(0).toInt >= minPL(0)) && (geno(PL).split(",")(1).toInt >= minPL(1)) && (geno(PL).split(",")(2).toInt >= minPL(2))))) 
) {
print("\t" + line(i))
} else {
print("\t" + "./.")
}
}
}// E for
print("\n")
line = input.readLine.split("\t")
}//E while
input.close
} // E Def Main

}//E Object