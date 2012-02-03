from configuration import *

class Genetic_Algorithm:
    self.generation_counter = 0
    
    def __init__(self, maxg):
        self.max_generations = maxg
    
    def runsim(self, population, landscape):
        """This method implements the main loop of the genetic algorithm."""
        for i in range(0, MAX_GA_GENS):
            """Reproduce the population based on fitness"""
            population.do_reproduction()
            
            """Mutate the population"""
            population.do_mutations()
            
            """Advance the generation counter"""
            self.generation_counter += 1
