#####################################################
#
# KO.py is derived from genetic_algorithm.py (GA)
# Whereas the GA module parallelizes over individuals,
# the KO module parallelizes over genes.
#
####################################################

from configuration import *
from landscape import *
from population import *
from plots import *
from logs import *

class KO:        
    genome = None
    landscape = None    
        
    def __init__(self, ap):
        id = ap.params["runid"] 
    
    def runko(self, ap):
        """This is the main method of the KO module"""
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        if rank == 0:
            self.runko_master(comm, ap)
        else:
            self.runko_slave(rank, comm, ap)

    def get_fitness_stats(self, gid_fitness):
        """gid_fitness is a hashtable, key = genome.id, value = fitness for that genome."""
        f_vals = gid_fitness.values()
        max_f = max(f_vals)
        min_f = min(f_vals)
        mean_f = mean(f_vals)
        median_f = mean(f_vals)
        std_f = std(f_vals)
        return [max_f, min_f, mean_f, median_f, std_f]

    def runko_master(self, comm, ap):
        rank = comm.Get_rank()
        
        print "\n. KO.py 44: Rank", rank, "enetered runko_master"
        wt_gid_fitness = {} # key = gene, value = KO'd fitness
        ko_gid_fitness = {} # key = gene, value = KO'd fitness

        for slave in range(1, comm.Get_size()):
            if ap.params["verbosity"] > 2:
                print "\n. KO.py line 40 - Master is waiting for data from slave", slave, "of", comm.Get_size()
            [their_wt_gid_fitness, their_ko_gid_fitness] = comm.recv(source=slave, tag=12)

            #if ap.params["verbosity"] > 50:
            #    print "\n. genetic_algorithm.py line 73 - Master received from slave", slave, ":\their_wt_gid_fitness:", their_wt_gid_fitness, "\their_ko_gid_fitness:", their_ko_gid_fitness
            for gid in their_wt_gid_fitness.keys():
                wt_gid_fitness[gid] = their_wt_gid_fitness[gid]
            for gid in their_ko_gid_fitness.keys():
                ko_gid_fitness[gid] = their_ko_gid_fitness[gid]
    
        """Print LOGS/koX.txt with the fitness effect of each KO'd gene for genome X."""          
        fout = open(ap.params["workspace"] + "/" + ap.params["runid"] + "/LOGS/ko.gen" + ap.params["generation"].__str__() + ".id" + ap.params["kogenome"].__str__() + ".txt", "w")
        fout.write("Gene\tWT f\tKO f\tWT-KO\n")
        gids = ko_gid_fitness.keys()
        gids.sort()
        for gid in gids:
            line = gid.__str__()
            line += "\t%.3f"%(wt_gid_fitness[gid])
            line += "\t%.3f"%(ko_gid_fitness[gid])
            line += "\t%.3f"%(wt_gid_fitness[gid]-ko_gid_fitness[gid]) 
            line += "\n"
            fout.write(line)
        fout.close()


    def runko_slave(self, rank, comm, ap):                                
        """gid_fitness[gene ID] = [fitness at generation i]"""
        wt_gid_fitness = {}
        ko_gid_fitness = {}
        
        """Find my_items, which is a list of genes that this slave is responsible for."""
        my_items = list_my_items(ap.params["trlist"], rank)
        print "\n. KO.py 82: Rank", rank, "has these items:", my_items
            
        for gid in my_items:
            """Calculate WT fitness for gene gid."""
            wt_gid_fitness[ gid ] = self.landscape.get_fitness( self.genome, ap)
            if ap.params["verbosity"] >= 2:
                """and then plot the expression of all genes in this genome"""
                title = "expr.wt.gen" + ap.params["generation"].__str__() + ".gene" + gid.__str__() 
                filenameseed = ap.params["workspace"] + "/" + ap.params["runid"] + "/" + EXPR_PLOTS + "/" + title
                write_expression_cran( self.genome, filenameseed, title, "time", "expression", ap)
            
            """Calculate KO fitness for gene gid."""
            ko_gid_fitness[ gid ] = self.landscape.get_fitness( self.genome, ap, ko=[gid])
            if ap.params["verbosity"] >= 2:
                """and then plot the expression of all genes in this genome"""
                title = "expr.ko.gen" + ap.params["generation"].__str__() + ".gene" + gid.__str__() 
                filenameseed = ap.params["workspace"] + "/" + ap.params["runid"] + "/" + EXPR_PLOTS + "/" + title
                write_expression_cran( self.genome, filenameseed, title, "time", "expression", ap)
    
        #print "\n. KO.94: slave", rank,"sending to master."
        """Send data to master."""
        comm.send( [wt_gid_fitness, ko_gid_fitness ], dest=0, tag=12)  
        
        
        