from configuration import *
from landscape import *
from population import *
from plots import *
from logs import *

class Genetic_Algorithm:        
    population = None
    landscape = None    
        
    def __init__(self, ap):
        id = ap.params["runid"] 
    
    def find_max_fit_gid(self, gid_fitness):
        maxf = 0.0
        maxgid = 0
        for gid in gid_fitness:
            if gid_fitness[gid] > maxf:
                maxf = gid_fitness[gid]
                maxgid = gid
        return maxgid
    
    def masteronly_bcast(self, ap):
        """Broadcast the updated population to MPI slaves."""
        pop_data = self.population.collapse()
        pop_data_pickle = pickle.dumps( pop_data )
        for slave in range(1, comm.Get_size()):
            comm.send(pop_data_pickle, dest=slave, tag=11)
    
        if ap.params["verbosity"] >= 2:
            pop_pickle_path = ap.params["workspace"] + "/" + ap.params["runid"] + "/" + POPPICKLES + "/population.gen" + ap.params["generation"].__str__() + ".pickle"
            print "\n. Saving the population to", pop_pickle_path
            fout = open(pop_pickle_path, "w")
            fout.write(pop_data_pickle)
            fout.close()
    
    def slaveonly_recv(self):
        pop_data_pickle = comm.recv(source=0, tag=11)
        if pop_data_pickle == None:
            print "\n. Node", comm.Get_rank(), "is finishing."
            exit()
        pop_data = pickle.loads( pop_data_pickle ) 
        population = Population()
        population.uncollapse(pop_data)
        self.population = population
    
    def masteronly_whistle(self):
        """Tell all slaves that we're done.  They should exit()."""
        for slave in range(1, comm.Get_size()):
            comm.send(None, dest=slave, tag=11)
            
    
    def runsim_master(self, comm, ap):                     
        """These variables will store time measurements.  Use --perftime True to print these stats."""
        notime = 0.0
        sumtime_gen = notime
        sumtime_end_calc = notime
        sumtime_end_gather = notime
        sumtime_end_stats = notime
        sumtime_getmm = notime
        sumtime_end_evo = notime
        sumtime_end_bcast = notime
        
        if ap.params["stopconvergence"]:
            count_conv_gens = 0
        
        """For each GA generation. . . ."""
        for i in range(ap.params["generation"], ap.params["generation"]+ap.params["maxgens"]):
            ap.params["generation"] = i
            time_start_gen = datetime.utcnow()
                                                        
            gid_fitness = {}

            """Wait for data from slaves. . ."""
            for slave in range(1, comm.Get_size()):
                if ap.params["verbosity"] > 50:
                    print "\n. genetic_algorithm.py line 68 - Master is waiting for data from slave", slave
                [their_gid_fitness, their_gid_terminal_expression] = comm.recv(source=slave, tag=11)
                if slave == 1:
                    time_end_calc = datetime.utcnow()
                if ap.params["verbosity"] > 50:
                    print "\n. genetic_algorithm.py line 73 - Master received from slave", slave, ":\n. gid_fitness:", their_gid_fitness, "\n. gid_terminal_expression:", their_gid_terminal_expression
                for gid in their_gid_fitness.keys():
                    gid_fitness[gid] = their_gid_fitness[gid]
                    if ap.params["enable_epigenetics"] == True:
                        self.population.genomes[gid].gene_expr = their_gid_terminal_expression[gid]
            
            time_end_gather = datetime.utcnow()
                        
            """Get basic stats on the population's fitness distribution"""
            [max_f, min_f, mean_f, median_f, std_f] = self.get_fitness_stats(gid_fitness)

            """Print LOGS/genX.txt with the fitness of each individual at this generation."""            
            fout = open(ap.params["workspace"] + "/" + ap.params["runid"] + "/LOGS/gen" + i.__str__() + ".txt", "a")
            fout.write("ID\tfitness\n")
            gids = gid_fitness.keys()
            gids.sort()
            for gid in gids:
                mark = ""
                if gid_fitness[gid] == max_f:
                    mark = "\t*"
                fout.write(gid.__str__() + "\t%.6f"%(gid_fitness[gid]) + mark + "\n")
            fout.close()
            
            time_end_stats = datetime.utcnow()
            
            if ap.params["verbosity"] >= 1:
                timedelta = datetime.utcnow() - time_start_gen
                sumtime_gen += timedelta.total_seconds()
                print "\n................................................\n" 
                print "\n. Generation", i, "required %.3f"%timedelta.total_seconds(), "sec. to compute."
                print ". Mean generation time so far = %.3f"%(sumtime_gen/(i+1)) + " sec."
                print ""
                print "\t. Population Size =\t", self.population.genomes.__len__()
                print "\t. Maximum Fitness =\t%.3f"%max_f
                print "\t. Minimum Fitness =\t%.3f"%min_f
                print "\t. Mean Fitness =\t%.3f"%mean_f
                print "\t. Median Fitness =\t%.3f"%median_f
                print "\t. StDev of Fitness =\t%.3f"%std_f
                #
                # to-do: calculate effective pop. size
                #
                #print "\n. effective popsize=", self.population.effective_popsize()
                line = "gen " + i.__str__() + "\t" + "\tmaxf= %.3f"%max_f + "\tminf= %.3f"%min_f + "\tmeanf= %.3f"%mean_f + "\tmedianf= %.3f"%median_f + " \tstdf= %.3f"%std_f 
                log_generation(ap, line)
                                     
            time_getmm = datetime.utcnow()
            
            [min_fitness, max_fitness, sum_fitness, fitness_gid] = self.population.get_minmax_fitness(gid_fitness)
            """Mark the elite individuals..."""
            self.population.mark_elite(fitness_gid, max_fitness, min_fitness, ap)
            self.population.print_fitness(ap, gid_fitness)
            """Reproduce the population based on pre-mutation fitness"""                
            self.population.do_reproduction(gid_fitness, min_fitness, max_fitness, sum_fitness, ap)
            """Mutate the population"""
            self.population.do_mutations(ap)
            time_end_evo = datetime.utcnow()

            """Check for convergence on terminal conditions."""
            """First check user-specified conditions. . . """
            if ap.params["stopconvergence"] == True:
                if mean_f >= ap.params["fgoal"]:
                    count_conv_gens += 1
                else:
                    count_conv_gens = 0
                if count_conv_gens >= 3:
                    self.masteronly_whistle()
                    exit()
            """Next check if this is the final generation. . . """
            
            """Broadcast the updated population to MPI slaves."""
            self.masteronly_bcast(ap)
            time_end_bcast = datetime.utcnow()
        
            if ap.getOptionalArg("--perftime") == "on":
                print "\n. Performance data for generation", i, ". . ."
                ngen = i + 1
                sumtime_end_calc += (time_end_calc - time_start_gen).total_seconds()
                print "\t. mean time_end_calc %.3f sec."%(sumtime_end_calc / ngen)
                sumtime_end_gather += (time_end_gather - time_end_calc).total_seconds()
                print "\t. mean time_end_gather %.3f sec."%(sumtime_end_gather/ ngen)
                sumtime_end_stats += (time_end_stats - time_end_gather).total_seconds()
                print "\t. mean time_end_stats %.3f sec."%(sumtime_end_stats/ ngen)
                sumtime_getmm += (time_getmm - time_end_stats).total_seconds()
                print "\t. mean time_getmm %.3f sec."%(sumtime_getmm/ ngen)
                sumtime_end_evo += (time_end_evo - time_getmm).total_seconds()
                print "\t. mean time_end_evo %.3f sec."%(sumtime_end_evo/ ngen)
                sumtime_end_bcast += (time_end_bcast - time_end_evo).total_seconds()
                print "\t. mean time_end_bcast %.3f sec."%(sumtime_end_bcast/ ngen)
                                    
                print "\t. mean comm/calc ratio: %.3f"%( (sumtime_end_gather + sumtime_end_bcast) / (sumtime_end_gather + sumtime_end_bcast + sumtime_end_calc + sumtime_end_stats + sumtime_getmm + sumtime_end_evo) )

    def runsim_slave(self, rank, comm, ap):           
        """For each GA generation. . . ."""
        for i in range(ap.params["generation"], ap.params["generation"]+ap.params["maxgens"]):
            ap.params["generation"] = i
            
            """gid_fitness[genome ID] = [fitness at generation i]"""
            gid_fitness = {}
            
            """Find my items, a list of individuals in the population for which I am responsible."""
            my_items = list_my_items(self.population.list_genome_ids(), rank)
        
            """Calculate fitness for every individual."""
            for gid in my_items:
                gid_fitness[ gid ] = self.landscape.get_fitness( self.population.genomes[gid], ap)
                if ap.params["verbosity"] >= 2:
                    """and then plot the expression of all genes in this genome"""
                    title = "expr.gen" + i.__str__() + ".gid" + gid.__str__() 
                    filenameseed = ap.params["workspace"] + "/" + ap.params["runid"] + "/" + EXPR_PLOTS + "/" + title
                    write_expression_cran( self.population.genomes[gid], filenameseed, title, "time", "expression", ap)
            
            """Save the expression at the last timeslice, to be inherited by childern."""
            gid_terminal_expression = {}
            if ap.params["enable_epigenetics"] == True:
                last_timeslice = ap.params["maxtime"]
                for gid in my_items:
                    gid_terminal_expression[gid] = {}
                    for geneid in range(0, self.population.genomes[gid].genes.__len__()):
                        gid_terminal_expression[gid][geneid] = [self.population.genomes[gid].gene_expr[geneid][last_timeslice]]
                                    
            """Send data to master."""
            comm.send( [gid_fitness, gid_terminal_expression], dest=0, tag=11)  
            
            """Get updated data from master."""
            self.slaveonly_recv()

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
        