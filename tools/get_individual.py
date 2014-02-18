from poptools import *
from argparser import *
ap = ArgParser(sys.argv)

outputdir = ap.getArg("--outdir") #outputdir is the folder into which a SimGenome run placed output.

generation = ap.getArg("--gen")
if generation != False:
    generation = int(generation)
else:
    print "You need to specify a generation with --gen"
    exit()

indi = ap.getArg("--id")
if indi != False:
    indi = int(indi)
else:
    print "You need to specify an individual with --id"
    exit()

p = get_individual(outputdir + "/POPS/pop.gen" + generation.__str__() + ".save.txt", indi)

print "N Genomes: 1"

for l in p[0]:
    print l
for l in p[1]:
    print l
