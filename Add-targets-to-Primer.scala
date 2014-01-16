import java.io._

val input = new BufferedReader(new FileReader("generic-primers.fa.txt"))
val out =  new BufferedWriter(new FileWriter("test.etst"))
var cprime = ""

var cline = input.readLine
while(input.ready){
if (cline(0) == '>') {
out.write(cline + "\n")
//println(cline + "\n")
cprime = ""
}
cline = input.readLine
while (cline.size != 0 && cline(0) != '>'){
cprime = cprime + cline
cline = input.readLine
}
out.write(cprime.subSequence(0,(cprime.size/2)-1) + "[" + cprime.subSequence((cprime.size/2)-1,(cprime.size/2)+2) + "]" + cprime.subSequence((cprime.size/2)+2,cprime.size) + "\n")
cprime = ""
}
out.write(cprime.subSequence(0,(cprime.size/2)-1) + "[" + cprime.subSequence((cprime.size/2)-1,(cprime.size/2)+2) + "]" + cprime.subSequence((cprime.size/2)+2,cprime.size) + "\n")
out.close