from configuration import *
from landscape import *
from population import *
from plots import *

class Genetic_Algorithm:
    generation_counter = 0
        
    population = None
    landscape = None    
    simid = None
    gen_gid_fitness = None
        
    def __init__(self, ap):
        id = ap.getOptionalArg("--runid")
        if id != False:
            self.simid == id 
    
    def find_max_fit_gid(self, gid_fitness):
        maxf = 0.0
        maxgid = 0
        for gid in gid_fitness:
            if gid_fitness[gid] > maxf:
                maxf = gid_fitness[gid]
                maxgid = gid
        return maxgid
    
    def runsim(self):
        """This method implements the main loop of the genetic algorithm."""
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        self.gen_gid_fitness = []
        for i in range(0, MAX_GA_GENS):
            now = datetime.now()
            
            if rank == 0:
                if i > 0:
                    """Reproduce the population based on fitness"""
                    self.population.do_reproduction()
                
                    """Mutate the population"""
                    self.population.do_mutations()
                
                prog = ProgressBar(0, N_GENOMES, 50, mode='dynamic', char='#')
            
            gid_fitness = {}
        
            """Broadcast the population to MPI slaves."""
            self.population = comm.bcast(self.population, root=0)
        
            """Calculate fitness for every individual."""
            for genome in self.population.genomes:
                """. .  .but only compute for the individuals that have been assigned to my MPI slice."""
                if is_my_item(genome.id, self.population.genomes.__len__(), rank):
                    print "Calling get_fitness for gid", genome.id, "for rank", rank
                    gid_fitness[ genome.id ] = self.landscape.get_fitness( genome, i )
                    prog.increment_progress()
            
            if rank == 0:
                prog.finish()
            
            """Gather the population from MPI slaves."""
            if rank > 0:
                my_genomes = []
                for genome in self.population.genomes:
                    if is_my_item(genome.id, self.population.genomes.__len__(), rank):
                        my_genomes.append( genome )
                comm.send([my_genomes, gid_fitness], dest=0, tag=11)  
            
            elif rank == 0:
                for slave in range(1, comm.Get_size()):
                    [their_genomes, their_gid_fitness] = comm.recv(source=slave, tag=11)
                    for d in their_genomes:
                        self.population.genomes[d.id] = d
                        gid_fitness[d.id] = their_gid_fitness[d.id]
            
            """Post-fitness. . .print stats, record data, etc."""
            self.print_fitness_stats(gid_fitness)
            self.gen_gid_fitness.append( gid_fitness )
            maxgid = self.find_max_fit_gid(gid_fitness)
            filenameseed = "PLOTS/expression.gen" + i.__str__() + ".gid" + maxgid.__str__() 
            plot_expression( self.population.genomes[maxgid], filenameseed, filenameseed, "time", "expression")
            
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
        
        