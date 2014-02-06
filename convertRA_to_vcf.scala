import net.sf.samtools.util._
import java.io._
import scala.collection.mutable.HashMap

val input = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream("/Users/chhar0/Dropbox/PhD-Work/770k/770K-genotypes.vcf.gz"))))
val genos = new BufferedReader(new FileReader("Damona20140205_RA2.txt"))
val out = new BufferedWriter(new OutputStreamWriter(new BlockCompressedOutputStream("Damona-11k_v2.vcf.gz")))

val raGenos = new HashMap[String, Array[String]]

while (genos.ready){
var lineGenos = genos.readLine.split("\t")

raGenos += ("Chr" + lineGenos(1) + ":" + lineGenos(2)) -> lineGenos

}

def raConv(in: String): String = {
in match {
case "RR" => "0/0"
case "AA" => "1/1"
case "RA" => "0/1"
case "AR" => "0/1"
case _  => "./."
}
}


var vcfGeno = input.readLine

while (vcfGeno(1) == '#'){
out.write(vcfGeno + "\n")
vcfGeno = input.readLine
}
out.write("""##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">""" + "\n")

var clArray = vcfGeno.split("\t")

out.write(s"${clArray(0)}\t${clArray(1)}\t${clArray(2)}\t${clArray(3)}\t${clArray(4)}\t${clArray(5)}\t${clArray(6)}\t${clArray(7)}\t${clArray(8)}")
val ids = raGenos("ChrChr:Position").apply(3).split(" ")
val sizeLine = ids.size 
var counter = 0

while (counter < ids.size){
out.write("\t" + ids(counter))
counter += 1
//for (i <- ids){
//out.write( "\t" + i)
//}
}

out.write("\n")

while (input.ready){
clArray = input.readLine.split("\t")
var key = clArray(0) + ":" + clArray(1)

if (raGenos.contains(key)){
println(raGenos(key).apply(0))

var ras = raGenos(key).apply(3)
println(ras.split(' ').size + " " + sizeLine)

if (ras.split(' ').size == sizeLine || ras.split(' ').size == (sizeLine + 8)){

out.write(s"${clArray(0).toLowerCase}\t${clArray(1)}\t" + raGenos(key).apply(0) + s"\t${clArray(3)}\t${clArray(4)}\t${clArray(5)}\t${clArray(6)}\t${clArray(7)}\tGT")

val tmp = ras.split(' ')

counter = 0

while (counter < ids.size){
out.write( "\t" + raConv(tmp(counter)))
counter += 1
}

//for (i <- tmp){out.write( "\t" + raConv(i))}
out.write("\n")
//out.write(s"\t${raConv(.apply(0))}\t${raConv(ras.split(' ').apply(2))}\t${raConv(ras.split(' ').apply(1))}\n")
} else {
println("ERROR\t\t\t" + key + " " + ras)
}
} else {
println("UNKOWN KEY: " + key)
}

}

out.close
