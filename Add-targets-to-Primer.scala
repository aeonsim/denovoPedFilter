import java.io._

val input = new BufferedReader(new FileReader("primer_dna.fa"))
val out =  new BufferedWriter(new FileWriter("test.etst"))
var cprime = ""

var cline = input.readLine
while(input.ready){
if (cline(0) == '>') {
out.write(cline + "\n")
println(cline)
cprime = ""
}
cline = input.readLine
while (input.ready && cline.size != 0 && cline(0) != '>'){
//println(cprime)
cprime = cprime + cline
cline = input.readLine
}
//println(cprime.size/2)
out.write(cprime.subSequence(0,(cprime.size/2)-30) + "[" + cprime.subSequence((cprime.size/2)-30,(cprime.size/2)+31) + "]" + cprime.subSequence((cprime.size/2)+31,cprime.size) + "\n")
cprime = ""
}
//out.write(cprime.subSequence(0,(cprime.size/2)-30) + "[" + cprime.subSequence((cprime.size/2)-30,(cprime.size/2)+31) + "]" + cprime.subSequence((cprime.size/2)+31,cprime.size) + "\n")
out.close