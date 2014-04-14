object hetInfo{

import java.io._
import net.sf.samtools.util._
import scala.collection.mutable.HashMap

def main (args: Array[String]) : Unit = {

if (args.size != 3){
println("hetInfo test.vcf.gz truth.vcf.gz SAMPLEID varcaller{gatk, fb, platy}")
System.exit(1)
}

val in_test = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream(args(0)))))
val in_truth = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream(args(1)))))

val animalID = args(2)
var truthAnimals = new HashMap[String,Array[String]]
var line = in_truth.readLine.split("\t")
var AANum, AADepth, CCNum, CCDepth, GGNum, GGDepth, TTNum, TTDepth, error = 0
var rAANum, rAADepth, rCCNum, rCCDepth, rGGNum, rGGDepth, rTTNum, rTTDepth, correct, incorrect = 0

val platform = args(3)

while (line.size < 7) line = in_truth.readLine.split("\t")

val truthLoc = line.indexOf(animalID)

while (in_truth.ready){
line = in_truth.readLine.split("\t")
truthAnimals += (line(0) + ":" + line(1)) -> line(truthLoc).split(":")
}

in_truth.close

line = in_test.readLine.split("\t")
while (line.size < 7) line = in_test.readLine.split("\t")

val testLoc = line.indexOf(animalID)

while (in_test.ready){
line = in_test.readLine.split("\t")
if (line(4).size == 1){
val details = line(testLoc).split(":")
val format = line(8).split(":")
val GT = details(format.indexOf("GT"))
var ALT, REF = 0

platform match {
case "fb" =>  ALT = details(format.indexOf("AO")).toInt ; REF = details(format.indexOf("RO")).toInt
case "gatk" => ALT = details(format.indexOf("AD")).split(",")(1).toInt ; REF = details(format.indexOf("AD")).split(",")(0).toInt
case "platy" => ALT = details(format.indexOf("NV")).toInt ; REF = details(format.indexOf("NR")).toInt - ALT
case _ => println("Invalid platform " + platform) ; System.exit(1)
}


//val DTH = details("DP")

val trueGT = truthAnimals(line(0) + ":" + line(1))(0)
if(trueGT == "0/1" || trueGT == "1/0"){
println(s"${trueGT}\t${line(3)}\t${REF}\t${line(4)}\t${ALT}\t${ALT/REF.toFloat}\t${ALT + REF}\t${GT}")
if (GT == "0/1" || GT == "1/0") correct += 1 else incorrect += 1
} else {
if (trueGT == "1/1"){
line(4) match {
case "A" => AANum += 1; AADepth += ALT;
case "C" => CCNum += 1; CCDepth += ALT;
case "G" => GGNum += 1; GGDepth += ALT;
case "T" => TTNum += 1; TTDepth += ALT;
case _ => error += 1
}

} else {
if (trueGT == "0/0"){
line(3) match {
case "A" => rAANum += 1; rAADepth += REF;
case "C" => rCCNum += 1; rCCDepth += REF;
case "G" => rGGNum += 1; rGGDepth += REF;
case "T" => rTTNum += 1; rTTDepth += REF;
case _ => error += 1
}

}
}

}
}
}
println(s"REFs\tNUM\tavgDEPTH\tALT NUM\tAlt avgDEPTH")
println(s"AA\t${rAANum}\t${rAADepth/rAANum}\t${AANum}\t${AADepth/AANum}")
println(s"CC\t${rCCNum}\t${rCCDepth/rCCNum}\t${CCNum}\t${CCDepth/CCNum}")
println(s"GG\t${rGGNum}\t${rGGDepth/rGGNum}\t${GGNum}\t${GGDepth/GGNum}")
println(s"TT\t${rTTNum}\t${rTTDepth/rTTNum}\t${TTNum}\t${TTDepth/TTNum}")
println(s"Correct Calls\t${correct}")
println(s"Incorrect Calls\t${incorrect}")
println(s"% Correct\t${correct/(incorrect + correct)}")
println(s"Errors\t${error}")
}


}