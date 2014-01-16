import net.sf.samtools.util._
import java.io._

val input = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream("beagle4-Chr1-phased-exCheckpoint.vcf.gz"))))
val out = new BufferedWriter(new OutputStreamWriter(new BlockCompressedOutputStream("out-gt.vcf.gz")))

var GTpos,GQpos = 0

def setFields(geno: String): Unit = {
//GQpos = geno.split(":").indexOf("GQ")
GTpos = geno.split(":").indexOf("GT")
}

var current = input.readLine.split("\t")
while (current(0).apply(0) == '#'){
out.write(current.reduceLeft{(a,b) => a + "\t" + b} + "\n")
current = input.readLine.split("\t")
}

while (input.ready){
setFields(current(8))
out.write(s"${current(0)}\t${current(1)}\t${current(2)}\t${current(3)}\t${current(4)}\t${current(5)}\t${current(6)}\t.\tGT:GQ")
for (i <- (9 to (current.size - 1))){
val genos = current(i).split(":")
if (genos.size >= 2){
//out.write("\t" + genos(GTpos) + ":" + genos(GQpos))
out.write("\t" + genos(GTpos))
} else {
out.write("\t.:.")
}
}
out.write("\n")
current = input.readLine.split("\t")
}

out.close
