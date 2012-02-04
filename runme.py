from argparser import *
from configuration import *
from genetic_algorithm import *
from population import *

def main():
    ap = ArgParser(sys.argv)
    
    """Verify consistency of command-line arguments"""
        
    """Build a population"""
    population = Population()
    
    """Build a random fitness landscape"""
    landscape = Landscape(genome = population)

    return [population, landscape]
    
    
    
    