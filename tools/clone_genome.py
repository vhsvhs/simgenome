"""
This is a very useful script that will extract an individual from a given saved population,
clone that individual N times, and then write a new population with the clones.
"""

from poptools import *
from argparser import *
ap = ArgParser(sys.argv)

poppath = ap.getArg("--poppath")
genome = int( ap.getArg("--id") )
n = int(ap.getArg("--nclones"))

it = get_individual(poppath, genome)
its = clone_individual(it, genome, n)
#N Genomes: 40
print "N Genomes: " + n.__str__()
for l in its:
    print l