object hetInfo{

import java.io._
import net.sf.samtools.util._
import scala.collection.mutable.HashMap

def main (args: Array[String]) : Unit = {

if (args.size != 2){
println("hetInfo test.vcf.gz truth.vcf.gz")
System.exit(1)
}

val in_test = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream(args(0)))))
val in_truth = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream(args(1)))))

val animalID = "23934418"
var truthAnimals = new HashMap[String,Array[String]]
var line = in_truth.readLine.split("\t")
var AANum, AADepth, CCNum, CCDepth, GGNum, GGDepth, TTNum, TTDepth = 0

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
val details = line(testLoc).split(":")
val GT = details("GT")
val REF = details("NR").toInt
val ALT = details("NA").toInt
//val DTH = details("DP")

val trueGT = truthAnimals(line(0) + ":" + line(1))(0)

if(trueGT == "0/1" || GT == "1/0"){
println(s"${line(3)}\t${REF}\t${line(4)}\t${ALT}\t${ALT/REF.toFloat}\t${GT}")

} else {
if (trueGT == "1/1"){
line(4) match {
case "A" => AANum += 1; AADepth += ALT;
case "C" => CCNum += 1; CCDepth += ALT;
case "G" => GGNum += 1; GGDepth += ALT;
case "T" => TTNum += 1; TTDepth += ALT;
}

} else {
if (trueGT == "0/0"){
line(3) match {
case "A" => AANum += 1; AADepth += REF;
case "C" => CCNum += 1; CCDepth += REF;
case "G" => GGNum += 1; GGDepth += REF;
case "T" => TTNum += 1; TTDepth += REF;
}

}
}

}

}


}


}