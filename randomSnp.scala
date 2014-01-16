import java.io._

val sampleSize = 100
val input = new BufferedReader(new FileReader("770K-genotypes.vcf"))
var snps:List[String] =  Nil
val random = new java.util.Random

var cline = input.readLine.split("\t")

while (cline(0).apply(0) == '#'){
cline = input.readLine.split("\t")
}

while (input.ready){
snps = (cline(0) + "\t" + cline(1) + "\t" + cline(1)) :: snps
cline = input.readLine.split("\t")
}

val aSnps = snps.toArray


var count = 0

while (count < sampleSize){
println(aSnps(random.nextInt(aSnps.size)))
count += 1
}
