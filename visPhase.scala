/*
Visualise Phase/Haplotypes
Need to know Sire, Dam, Proband, Children
Label Haplotypes S1/S2 D1/D2
Maternal/Paternal if from Beagle, Paternal/Maternal if from RTG pop
*/

import java.io._
import scala.collection.mutable.HashMap
import net.sf.samtools.util._

val input = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream("exTrio.50KHD.DP6-30GQ10SEQ.phased.b4.vcf.gz"))))

var count = 0
var sites = 0

val pro = 12
val child = 16
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