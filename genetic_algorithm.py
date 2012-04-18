from configuration import *
from landscape import *
from population import *
from plots import *
from logs import *

class Genetic_Algorithm:        
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
    
    def runsim_master(self, comm, ap):
        self.gen_gid_fitness = []
                
        """For each GA generation. . . ."""
        for i in range(0, MAX_GA_GENS):
            self.landscape.gen_counter = i
            now = datetime.now()
                                
            """gid_fitness[genome ID] = [fitness at generation i]"""
            gid_fitness = {}

            """Find my items."""
            my_items = list_my_items(self.population.list_genome_ids(), 0)
        
            """Calculate fitness for every individual."""
            for gid in self.population.genomes.keys():
                """. .  .but only compute for the individuals that have been assigned to my MPI slice."""
                if my_items.__contains__(gid):
                    #print "Genome", gid, "has", self.population.genomes[gid].genes.__len__(), "genes."
                    gid_fitness[ gid ] = self.landscape.get_fitness( self.population.genomes[gid], i , ap)
                    #prog.increment_progress()
            
            #prog.finish()
            
            """Wait for data from slaves."""
            for slave in range(1, comm.Get_size()):
                [their_genomes, their_gid_fitness] = comm.recv(source=slave, tag=11)
                for gid in their_genomes.keys():
                    self.population.genomes[gid] = their_genomes[gid]
                    gid_fitness[gid] = their_gid_fitness[gid]
            
            """Post-fitness. . .print stats, record data, etc."""
            if int(ap.getOptionalArg("--verbose")) > 2:
                pop_data = self.population.collapse()
                pop_data_pickle = pickle.dumps( pop_data )
                fout = open(ap.getArg("--runid") + "/" + POPPICKLES + "/population.gen" + i.__str__() + ".pickle", "w")
                fout.write(pop_data_pickle)
                fout.close()
                        
            self.gen_gid_fitness.append( gid_fitness )
            maxgid = self.find_max_fit_gid(gid_fitness)
            filenameseed = ap.getArg("--runid") + "/" + EXPR_PLOTS + "/expr.gen" + i.__str__() + ".gid" + maxgid.__str__() 
            plot_expression( self.population.genomes[maxgid], filenameseed, filenameseed, "time", "expression", ap)
            
            [max_f, min_f, mean_f, median_f, std_f] = self.get_fitness_stats(gid_fitness)
            
            if int(ap.getOptionalArg("--verbose")) > 2:
                timedelta = datetime.now() - now
                print "\n................................................\n" 
                print "\n. Generation", i, "required", timedelta, "to compute."
                print ""
                print "\t. population size =\t", self.population.genomes.__len__()
                print "\t. maximum fitness =\t%.3f"%max_f
                print "\t. minimum fitness =\t%.3f"%min_f
                print "\t. mean fitness =\t%.3f"%mean_f
                print "\t. median fitness =\t%.3f"%median_f
                print "\t. s.d. fitness =\t%.3f"%std_f
                #print "\n. effective popsize=", self.population.effective_popsize()
                line = "gen " + i.__str__() + "\t" + "\tmaxf=%.3f"%max_f + "\tminf= %.3f"%min_f + "\tmeanf= %.3f"%mean_f + "\tmedianf= %.3f"%median_f + "\tstdf= %.3f"%std_f 
                log_generation(ap, line)
            
            """Check for convergence on terminal conditions."""
            if mean_f == 1.0:
                print "The population arrived at an optima.  Goodbye."
                exit(1)
            
            if i < MAX_GA_GENS:                
                [min_fitness, max_fitness, sum_fitness, fitness_gid] = self.population.get_minmax_fitness(gid_fitness)
                """Mark the elite individuals..."""
                self.population.mark_elite(fitness_gid, max_fitness, min_fitness, ap)
                self.population.print_fitness(ap, gid_fitness)
                """Reproduce the population based on pre-mutation fitness"""
                self.population.do_reproduction(gid_fitness, min_fitness, max_fitness, sum_fitness, ap)
                """Mutate the population"""
                self.population.do_mutations(ap)
                """Broadcast the population to MPI slaves."""
                pop_data = self.population.collapse()
                pop_data_pickle = pickle.dumps( pop_data )
                for slave in range(1, comm.Get_size()):
                    comm.send(pop_data_pickle, dest=slave, tag=11)


    def runsim_slave(self, rank, comm, ap):
        self.gen_gid_fitness = []
        
        """For each GA generation. . . ."""
        for i in range(0, MAX_GA_GENS):        
            """gid_fitness[genome ID] = [fitness at generation i]"""
            gid_fitness = {}

            """Find my items."""
            my_items = list_my_items(self.population.list_genome_ids(), rank)
        
            """Calculate fitness for every individual."""
            for gid in my_items:
                gid_fitness[ gid ] = self.landscape.get_fitness( self.population.genomes[gid], i , ap)
                #prog.increment_progress()
                        
            """Send data to master."""
            my_genomes = {}
            for gid in my_items:
                my_genomes[gid] = self.population.genomes[gid]
            comm.send([my_genomes, gid_fitness], dest=0, tag=11)  
            
            if i < MAX_GA_GENS:
                """Get updated data from master."""
                pop_data_pickle = comm.recv(source=0, tag=11)
                pop_data = pickle.loads( pop_data_pickle ) 
                population = Population()
                population.uncollapse(pop_data)
                self.population = population

    def runsim(self, ap):
        """This is the main method of the genetic algorithm."""
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        if rank == 0:
            self.runsim_master(comm, ap)
        else:
            self.runsim_slave(rank, comm, ap)
            
 
    def get_fitness_stats(self, gid_fitness):
        """gid_fitness is a hashtable, key = genome.id, value = fitness for that genome."""
        f_vals = gid_fitness.values()
        max_f = max(f_vals)
        min_f = min(f_vals)
        mean_f = mean(f_vals)
        median_f = mean(f_vals)
        std_f = std(f_vals)
        return [max_f, min_f, mean_f, median_f, std_f]
        