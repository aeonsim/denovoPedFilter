import java.io._
import scala.collection.mutable.HashMap

val sire = new BufferedReader(new FileReader("allsirefam1.txt"))
val dam = new BufferedReader(new FileReader("alldamfam1.txt"))

val svg = new BufferedWriter(new FileWriter("test-out1.svg"))

val colour = HashMap("SM" -> "red","SF" -> "blue", "DF" -> "green", "DM" -> "orange", "D?" -> "white", "S?" -> "white")
val chrSize = HashMap("Chr1" -> 159000000, "Chr2" -> 138000000, "Chr3" -> 122000000,"Chr4" -> 121000000,"Chr5" -> 122000000,"Chr6" -> 120000000,"Chr7" -> 113000000,
"Chr8" -> 114000000,"Chr9" -> 106000000,"Chr10" -> 105000000,"Chr11" -> 108000000,"Chr12" -> 92000000,"Chr13" -> 85000000,"Chr14" -> 85000000,"Chr15" -> 86000000,
"Chr16" -> 82000000,"Chr17" -> 76000000,"Chr18" -> 67000000,"Chr19" -> 65000000,"Chr20" -> 73000000,"Chr21" -> 72000000,"Chr22" -> 62000000,"Chr23" -> 53000000,
"Chr24" -> 63000000,"Chr25" -> 43000000,"Chr26" -> 52000000,"Chr27" -> 46000000,"Chr28" -> 47000000,"Chr29" -> 52000000, "ChrX" -> 149000000)


val filter = 0

var blockStart = 0

var DamY, SireY, SireX = 30
var DamX = 110
var curPosY = SireY

var curChr = "Chr1"
var prevColour = ""

var blockEnd, blockSize, sireCount, damCount = 0

def drawShape(size:Int, X:Int, Y:Int, Type: String): String = {
s"""<rect width="50" height="${size}" x="${X}" y="${Y + blockStart}" style="fill:${colour(Type)};fill-opacity:0.8" />\n"""
}

def drawEdge(size:Int, X:Int, Y:Int): String = {
s"""<rect width="50" height="${size}" x="${X}" y="${Y}" style="fill:none;stroke:black" />\n"""
}

//SETUP File:
svg.write("<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"3500\" height=\"2000\">\n")
svg.write("""<text x="30" y="13" fill="red">Chr1</text>""")
svg.write("""<text x="30" y="28" fill="black">Paternal</text>""")
svg.write("""<text x="110" y="28" fill="black">Maternal</text>""")
svg.write(drawEdge(chrSize(curChr)*6 / 1000000,SireX,SireY))


while (sire.ready){

var cline = sire.readLine().split("\t")

if (cline(0) == curChr){ //Same Chromosome

if (cline(4).toInt > filter){ //Block is not 1 SNP
prevColour = cline(3)
blockEnd = cline(2).toInt * 6 / 1000000
blockSize = blockEnd - blockStart
if (cline(3) == prevColour){
svg.write(drawShape(blockSize,SireX,SireY,prevColour))
} else {
svg.write(drawShape(blockSize,SireX,SireY,cline(3)))
}
blockStart = blockEnd
}
} else {
blockEnd = chrSize(curChr) * 6 / 1000000
blockSize = blockEnd - blockStart
svg.write(drawShape(blockSize,SireX,SireY,prevColour))
//svg.write(drawEdge(chrSize(curChr)*6 / 1000000,SireX,SireY))
prevColour = cline(3)
blockStart = 0
blockSize = 0
blockEnd = 0
SireX = SireX + 200
curChr = cline(0)
if (cline(0) == "Chr15"){
SireY = SireY + 1000
SireX = 30
}
svg.write(s"""<text x="${SireX}" y="${SireY - 15}" fill="black">${curChr}</text>""")
svg.write(s"""<text x="${SireX}" y="${SireY -2}" fill="black">Paternal</text>""")
svg.write(s"""<text x="${SireX + 80}" y="${SireY -2}" fill="black">Maternal</text>""")
}//E else if Current
svg.write(drawEdge(chrSize(curChr)*6 / 1000000,SireX,SireY))
}//E while

blockEnd = 0
blockStart = 0
blockSize = 0
curChr = "Chr1"
prevColour = ""

svg.write("########### DAM ###########\n")

while (dam.ready){

var cline = dam.readLine().split("\t")

if (cline(0) == curChr){ //Same Chromosome

if (cline(4).toInt > filter){ //Block is not 1 SNP
prevColour = cline(3)
blockEnd = cline(2).toInt * 6 / 1000000
blockSize = blockEnd - blockStart
svg.write(drawShape(blockSize,DamX,DamY,cline(3)))
blockStart = blockEnd
}
} else {
blockEnd = chrSize(curChr) * 6 / 1000000
blockSize = blockEnd - blockStart
svg.write(drawShape(blockSize,DamX,DamY,prevColour))
//svg.write(drawEdge(chrSize(curChr)*6 / 1000000,DamX,DamY))
prevColour = cline(3)
blockStart = 0
blockSize = 0
blockEnd = 0
DamX = DamX + 200
curChr = cline(0)
if (cline(0) == "Chr15"){
DamY = DamY + 1000
DamX = 110
}
}//E else if Current
svg.write(drawEdge(chrSize(curChr)*6 / 1000000,DamX,DamY))
}//E while

svg.write("""<text x="1500" y="850" fill="black">Legend</text>""")
svg.write("""<text x="1610" y="895" fill="black">Grand Sires</text>""")
svg.write(s"""<rect width="50" height="50" x="1500" y="855" style="fill:${colour("SF")};stroke:black" />""")
svg.write(s"""<rect width="50" height="50" x="1555" y="855" style="fill:${colour("DF")};stroke:black" />""")
svg.write("""<text x="1610" y="940" fill="black">Grand Dams</text>""")
svg.write(s"""<rect width="50" height="50" x="1500" y="910" style="fill:${colour("SM")};stroke:black" />""")
svg.write(s"""<rect width="50" height="50" x="1555" y="910" style="fill:${colour("DM")};stroke:black" />""")

svg.write("</svg>")
svg.close
