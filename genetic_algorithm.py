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
            now = datetime.now()
            
            """Reproduce the population based on fitness"""
            self.population.do_reproduction()
            
            """Mutate the population"""
            self.population.do_mutations()
            
            
            print "\n. Generation", i, ". . ."
            prog = ProgressBar(0, N_GENOMES, 50, mode='dynamic', char='#')
            
            gid_fitness = {}
            for genome in self.population.genomes:
                print "\n.Calculating fitness for individual", genome.id
                gid_fitness[ genome.id ] = self.landscape.get_fitness( genome )
                prog.increment_progress()
            prog.finish()
            
            self.print_fitness_stats(gid_fitness)
            
            """Advance the generation counter"""
            self.generation_counter += 1

            timedelta = datetime.now() - now
            print "Generation", i, "required", timedelta, "to compute."
    def print_fitness_stats(self, gid_fitness):
        """gid_fitness is a hashtable, key = genome.id, value = fitness for that genome."""
        f_vals = gid_fitness.values()
        max_f = max(f_vals)
        min_f = min(f_vals)
        mean_f = mean(f_vals)
        median_f = mean(f_vals)
        std_f = std(f_vals)
        print "Fitness... max=", max_f, "min=", min_f, "mean=", mean_f, "median=", median_f, "std=", std_f 
        
        