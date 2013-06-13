import java.io._
import net.sf.samtools.util._

object visPhase2{

def main(args: Array[String]) : Unit = {

println("VCF Parent Child")

val vcf = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream(args(0)))))
val bedout = new BufferedWriter(new FileWriter(args(0) + "." + args(2) + ".bedgraph"))

//GrandSire = 0 GrandDam = 1 Unkown = 0.5

var prevPhase = 0.5
var sw,dw,switches = 0

var cline = vcf.readLine.split("\t")

while (cline(0).apply(1) == '#')  cline = vcf.readLine.split("\t")

val animals = cline.toList.slice(9,cline.size)

val child = 9 + animals.indexOf(args(2))
val sire = 9 + animals.indexOf(args(1))

while (vcf.ready){
cline = vcf.readLine.split("\t")

if (cline(child).apply(0) == cline(sire).apply(0) && cline(child).apply(0) == cline(sire).apply(2)){
bedout.write(s"${cline(0)}\t${cline(1)}\t${cline(1)}\t${prevPhase.toString}\n")
} else {
if (cline(child).apply(0) == cline(sire).apply(0) && cline(child).apply(0) != cline(sire).apply(2)){
bedout.write(s"${cline(0)}\t${cline(1)}\t${cline(1)}\t1\n")
if (prevPhase != 1.0) {
sw += 1
dw -= 1
}
if (sw > 4){
sw = 0
dw = 0
prevPhase = 1.0
}
}
if (cline(child).apply(0) != cline(sire).apply(0) && cline(child).apply(0) == cline(sire).apply(2)){
bedout.write(s"${cline(0)}\t${cline(1)}\t${cline(1)}\t0\n")
if (prevPhase != 0.0){
sw -= 1
dw += 1
}
if (dw > 4) {
sw = 0
dw = 0
prevPhase = 0.0
}
}

}

}//Ewhile
bedout.close





} //Emain


} //Eobject