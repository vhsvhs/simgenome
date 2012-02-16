from argparser import *
from configuration import *
from genetic_algorithm import *
from population import *
from version import *
from debug_tools import *

def main():
    ap = ArgParser(sys.argv)
    
    splash()
    
    """Verify consistency of command-line arguments"""
        
    """Build a population"""
    population = Population()
    
    """Build a random fitness landscape"""
    landscape = Landscape(genome = population)
    
    ga = Genetic_Algorithm(population, landscape)
    ga.runsim()

    return [population, landscape]

def splash():
    print "======================================="
    print "."
    print ". Sim Genome"
    print "."
    print ". " + VERSIONSTRING
    print "."
    print "======================================="
    
#
# main:
#    
x = main()

    