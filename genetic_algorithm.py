from configuration import *
from landscape import *
from population import *

class Genetic_Algorithm:
    generation_counter = 0
        
    def runsim(self, population, landscape):
        """This method implements the main loop of the genetic algorithm."""
        for i in range(0, MAX_GA_GENS):
            """Reproduce the population based on fitness"""
            population.do_reproduction()
            
            """Mutate the population"""
            population.do_mutations()
            
            """Advance the generation counter"""
            self.generation_counter += 1
