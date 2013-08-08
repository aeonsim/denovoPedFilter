import java.io._
import net.sf.samtools.util._


object visPhase2{

def main(args: Array[String]) : Unit = {

println("VCF sire dam child")

val vcf = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream(args(0)))))
val sireout = new BufferedWriter(new FileWriter(args(0) + "." + args(3) + ".sire.bedgraph"))
val damout = new BufferedWriter(new FileWriter(args(0) + "." + args(3) + ".dam.bedgraph"))

//GrandSire = 0 GrandDam = 1 Unkown = 0.5

var prevPhase,dprevPhase = 0.5
var sw,dw,switches, dsw,ddw,dswitches = 0

var cline = vcf.readLine.split("\t")

while (cline(0).apply(1) == '#')  cline = vcf.readLine.split("\t")

val animals = cline.toList.slice(9,cline.size)

val child = 9 + animals.indexOf(args(3))
val dam = 9 + animals.indexOf(args(2))
val sire = 9 + animals.indexOf(args(1))

while (vcf.ready){
cline = vcf.readLine.split("\t")

if (sire != -1 && cline(child).apply(0) == cline(sire).apply(0) && cline(child).apply(0) == cline(sire).apply(2)){
sireout.write(s"${cline(0)}\t${cline(1)}\t${cline(1)}\t${prevPhase.toString}\n")
} else {
if (sire != -1 && cline(child).apply(0) == cline(sire).apply(0) && cline(child).apply(0) != cline(sire).apply(2)){
sireout.write(s"${cline(0)}\t${cline(1)}\t${cline(1)}\t1\n")
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
if (sire != -1 && cline(child).apply(0) != cline(sire).apply(0) && cline(child).apply(0) == cline(sire).apply(2)){
sireout.write(s"${cline(0)}\t${cline(1)}\t${cline(1)}\t0\n")
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

/*
*	DAM Code
*/

if (dam != -1 && cline(child).apply(2) == cline(dam).apply(0) && cline(child).apply(2) == cline(dam).apply(2)){
damout.write(s"${cline(0)}\t${cline(1)}\t${cline(1)}\t${prevPhase.toString}\n")
} else {
if (dam != -1 && cline(child).apply(2) == cline(dam).apply(0) && cline(child).apply(2) != cline(dam).apply(2)){
damout.write(s"${cline(0)}\t${cline(1)}\t${cline(1)}\t1\n")
if (dprevPhase != 1.0) {
dsw += 1
ddw -= 1
}
if (dsw > 4){
dsw = 0
ddw = 0
dprevPhase = 1.0
}
}
if (dam != -1 && cline(child).apply(2) != cline(dam).apply(0) && cline(child).apply(2) == cline(dam).apply(2)){
damout.write(s"${cline(0)}\t${cline(1)}\t${cline(1)}\t0\n")
if (dprevPhase != 0.0){
dsw -= 1
ddw += 1
}
if (ddw > 4) {
dsw = 0
ddw = 0
dprevPhase = 0.0
}
}

}


}//Ewhile
sireout.close
damout.close

} //Emain


} //Eobject