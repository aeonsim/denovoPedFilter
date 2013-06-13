/*
Visualise Phase/Haplotypes
Need to know Sire, Dam, Proband, Children
Label Haplotypes S1/S2 D1/D2
Maternal/Paternal if from Beagle, Paternal/Maternal if from RTG pop
*/

import java.io._
import scala.collection.mutable.HashMap
import net.sf.samtools.util._

val input = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream("duos-phased-beagle4.vcf.gz"))))

var count = 0
var sites = 0

val pro = 9
val child = 10
var cblock = ""

var cline = input.readLine.split("\t")
while (cline(0).apply(0) == '#') cline = input.readLine.split("\t")

while (count < 1000000){

if (!(cline(pro).split(":").apply(0).contains('/') || cline(child).split(":").apply(0).contains('/'))){

var pGT = cline(pro).split(":").apply(0).split("\\|")
var cGT = cline(child).split(":").apply(0).split("\\|")

if (cblock == "") print(cline(0) + "\t" + cline(1)) 

//Assuming Pat|Mat
if (cGT(0) == pGT(1)){ //Maternal source
sites += 1
if (cblock != "M") {
print("\t" + cline(1) + "\t" + "Maternal\t" + sites + "\n")
print(cline(0) + "\t" + cline(1))
cblock = "M"
sites = 0
} //Eif "M"

} else {
if (cGT(0) == pGT(0)){
sites += 1
if (cblock != "P") {
print("\t" + cline(1) + "\t" + "Paternal\t" + sites + "\n")
print(cline(0) + "\t" + cline(1))
cblock = "P"
sites = 0
}//Eif "P"

}//Eif cGT1 = pGT1

} //E else
}
cline = input.readLine.split("\t")
count += 1
} //E While














object visPhase{

def main (args: Array[String]): Unit = {
println("scala visPhase input.vcf.gz sire dam Proband")

val sire = args(1)
val dam = args(2)
val pro = args(3)

val input = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream(args(0)))))
val out = new BufferedWriter(new FileWriter(${args(0) + ".blocks.bed"))

var current = input.readLine.split("\t")
while (current(0).apply(0) == '#' & current(0).apply(1) == '#'){
out.write(current.reduceLeft{(a,b) => a + "\t" + b} + "\n")
current = input.readLine.split("\t")
} //E While

var anPos = new HashMap[String, Int]

for ( i <- (9 to (current.size -1))){
anPos += current(i) -> i
} //E for
current = input.readLine.split("\t")

var curPhaseSire, curPhaseDam = ""

while (input.ready){

val cline = current.split("\t")
proGeno = cline.apply(anPos(pro)).split("|")

cGeno = cline.apply(3).split("|")




current = input.readLine.split("\t")
}//E While


}//E Main

}// E Object