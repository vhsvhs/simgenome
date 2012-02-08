from configuration import *
from landscape import *
from population import *

class Genetic_Algorithm:
    generation_counter = 0
        
    population = None
    landscape = None
        
    def __init__(self, population, landscape):
        self.population = population
        self.landscape = landscape
        
    def runsim(self):
        """This method implements the main loop of the genetic algorithm."""
        for i in range(0, MAX_GA_GENS):
            """Reproduce the population based on fitness"""
            self.population.do_reproduction()
            
            """Mutate the population"""
            self.population.do_mutations()
            
            # for debugging:
            print "\n. Calculating fitness of each genome..."
            prog = ProgressBar(0, N_GENOMES, 50, mode='dynamic', char='#')
            for genome in self.population.genomes:
                self.landscape.get_fitness( genome )
                prog.increment_progress()
            prog.finish()
          
            """Advance the generation counter"""
            self.generation_counter += 1
