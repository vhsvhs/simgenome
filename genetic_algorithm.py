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
    
#    depricated
#    def find_max_fit_gid(self, gid_fitness):
#        maxf = 0.0
#        maxgid = 0
#        for gid in gid_fitness:
#            if gid_fitness[gid] > maxf:
#                maxf = gid_fitness[gid]
#                maxgid = gid
#        return maxgid
    
    def masteronly_bcast(self, ap):
        """Broadcast the updated population to MPI slaves."""
        pop_data = self.population.collapse()
        pop_data_pickle = pickle.dumps( pop_data )
        
        start = datetime.utcnow()
        for slave in range(1, comm.Get_size()):
            comm.send(pop_data_pickle, dest=slave, tag=11)
        ap.params["sumtime_comm"] += (datetime.utcnow() - start).total_seconds()
    
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
        ap.params["sumtime_comm"] = 0.0        
        ap.params["sumtime_gen"] = 0.0

        
        if ap.params["stopconvergence"]:
            count_conv_gens = 0
        
        """For each GA generation. . . ."""
        for i in xrange(ap.params["generation"], ap.params["generation"]+ap.params["maxgens"]):
            ap.params["generation"] = i
            time_start_gen = datetime.utcnow()
            ap.params["sumtime_generation"] = 0.0
                                                        
            gid_fitness = {} # key = genome ID, value = fitness of that genome
            slave_timing = {} # key = MPI slave ID, value = array of timing performance data

            """Wait for data from slaves. . ."""
            for slave in range(1, comm.Get_size()):
                if ap.params["verbosity"] > 50:
                    print "\n. genetic_algorithm.py line 68 - Master is waiting for data from slave", slave
                
                start = datetime.utcnow()
                [their_gid_fitness, their_gid_terminal_expression, their_timing] = comm.recv(source=slave, tag=11)
                ap.params["sumtime_comm"] += (datetime.utcnow() - start).total_seconds()
                
                slave_timing[slave] = their_timing
                
                if ap.params["verbosity"] > 50:
                    print "\n. genetic_algorithm.py line 73 - Master received from slave", slave, ":\n. gid_fitness:", their_gid_fitness, "\n. gid_terminal_expression:", their_gid_terminal_expression
                for gid in their_gid_fitness.keys():
                    gid_fitness[gid] = their_gid_fitness[gid]
                    if ap.params["enable_epigenetics"] == True:
                        self.population.genomes[gid].gene_expr = their_gid_terminal_expression[gid]
            
            time_start_fstats = datetime.utcnow()
            
            """Get basic stats on the population's fitness distribution"""
            [max_f, min_f, mean_f, median_f, std_f] = self.get_fitness_stats(gid_fitness)

            """Print LOGS/genX.txt with the fitness of each individual at this generation."""            
            fout = open(ap.params["workspace"] + "/" + ap.params["runid"] + "/" + FITNESS_DIR + "/gen" + i.__str__() + ".txt", "a")
            fout.write("ID\tfitness\n")
            gids = gid_fitness.keys()
            gids.sort()
            for gid in gids:
                mark = ""
                if gid_fitness[gid] == max_f:
                    mark = "\t*"
                fout.write(gid.__str__() + "\t%.6f"%(gid_fitness[gid]) + mark + "\n")
            fout.close()
                        
            if ap.params["verbosity"] >= 1:
                timedelta = datetime.utcnow() - time_start_gen
                ap.params["sumtime_gen"] += timedelta.total_seconds()
                print "\n................................................\n" 
                print "\n. Generation", i, "required %.3f"%timedelta.total_seconds(), "sec. to compute."
                print ". Mean generation time so far = %.3f"%(ap.params["sumtime_gen"]/(i+1)) + " sec."
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
            
            [min_fitness, max_fitness, sum_fitness, fitness_gid] = self.population.get_minmax_fitness(gid_fitness)
            """Mark the elite individuals..."""
            self.population.mark_elite(fitness_gid, max_fitness, min_fitness, ap)
            self.population.print_fitness(ap, gid_fitness)
            
            ap.params["sumtime_stats"] = (datetime.utcnow() - time_start_fstats).total_seconds()
            
            time_start_evo = datetime.utcnow()
            """Reproduce the population based on pre-mutation fitness"""                
            self.population.do_reproduction(gid_fitness, min_fitness, max_fitness, sum_fitness, ap)
            """Mutate the population"""
            self.population.do_mutations(ap)
            ap.params["sumtime_evo"] = (datetime.utcnow() - time_start_evo).total_seconds()

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
        
            if ap.getOptionalArg("--perftime") == "True":
                print "\n. Performance data for generation", i, ". . ."
                mean_slave = [0,0,0]
                for slave in slave_timing:
                    print "\t. Slave", slave, ":", slave_timing[slave]
                    mean_slave[0] += slave_timing[slave][0]
                    mean_slave[1] += slave_timing[slave][1]
                    mean_slave[2] += slave_timing[slave][2]
                n = slave_timing.keys().__len__() 
                print "\t. Mean slave:", mean_slave[0]/n, mean_slave[1]/n, mean_slave[2]/n
                
                print "\t. Total MPI communication time:", ap.params["sumtime_comm"]
            
                print "\t. Total generation time:", (datetime.utcnow() - time_start_gen).total_seconds()
                
                print "\t. Time for reproduction and mutation:", ap.params["sumtime_evo"]
                
                print "\t. Time for statistics:", ap.params["sumtime_stats"]
                
                #ngen = i + 1
                #ap.params["sumtime_end_calc"] += (time_end_calc - time_start_gen).total_seconds()
                #print "\t. mean time_end_calc %.3f sec."%(ap.params["sumtime_end_calc"] / ngen)
                #ap.params["sumtime_end_gather"] += (time_end_gather - time_end_calc).total_seconds()
                #print "\t. mean time_end_gather %.3f sec."%(ap.params["sumtime_end_gather"]/ ngen)
                #ap.params["sumtime_end_stats"] += (time_end_stats - time_end_gather).total_seconds()
                #print "\t. mean time_end_stats %.3f sec."%(ap.params["sumtime_end_stats"]/ ngen)
                #ap.params["sumtime_getmm"] += (time_getmm - time_end_stats).total_seconds()
                #print "\t. mean time_getmm %.3f sec."%(ap.params["sumtime_getmm"]/ ngen)
                #ap.params["sumtime_end_evo"] += (time_end_evo - time_getmm).total_seconds()
                #print "\t. mean time_end_evo %.3f sec."%(ap.params["sumtime_end_evo"]/ ngen)
                #ap.params["sumtime_end_bcast"] += (time_end_bcast - time_end_evo).total_seconds()
                #print "\t. mean time_end_bcast %.3f sec."%(ap.params["sumtime_end_bcast"]/ ngen)
                                    
                #print "\t. mean comm/calc ratio: %.3f"%( (ap.params["sumtime_end_gather"] + ap.params["sumtime_end_bcast"]) / (ap.params["sumtime_end_gather"] + ap.params["sumtime_end_bcast"] + ap.params["sumtime_end_calc"] + ap.params["sumtime_end_stats"] + ap.params["sumtime_getmm"]+ ap.params["sumtime_end_evo"]) )

    def runsim_slave(self, rank, comm, ap):           
        """For each GA generation. . . ."""
        for i in range(ap.params["generation"], ap.params["generation"]+ap.params["maxgens"]):
            ap.params["generation"] = i

            """We're going to time some things, for performance measurements."""
            ap.params["sumtime_ptables"] = 0.0
            ap.params["sumtime_calcprobtables"] = 0.0
            ap.params["sumtime_probexpr"] = 0.0

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
                    for geneid in xrange(0, self.population.genomes[gid].genes.__len__()):
                        gid_terminal_expression[gid][geneid] = [self.population.genomes[gid].gene_expr[geneid][last_timeslice]]
            
            #if ap.getOptionalArg("--perftime") == "True":
            timing = [ ap.params["sumtime_ptables"], ap.params["sumtime_calcprobtables"], ap.params["sumtime_probexpr"] ]
                #print "\n. Times for Node", comm.Get_rank(), "make ptables:",ap.params["sumtime_ptables"], "calc_ptables:", ap.params["sumtime_calcprobtables"], "pe:",ap.params["sumtime_probexpr"]
                        
            """SEND data to master."""
            comm.send( [gid_fitness, gid_terminal_expression, timing], dest=0, tag=11)  
            
            """RECIEVE updated data from master."""
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
        