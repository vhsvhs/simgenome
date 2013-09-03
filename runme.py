from argparser import *
from configuration import *
from genetic_algorithm import *
from KO import *
from population import *
from version import *
from cli import *

def check_workspace(ap):
    if False == os.path.exists(ap.params["workspace"]):
        os.system("mkdir " + ap.params["workspace"])
    if False == os.path.exists(ap.params["workspace"] + "/" + ap.params["runid"]):
        os.system("mkdir " + ap.params["workspace"] + "/" + ap.params["runid"])
    dirs = [POPPICKLES, "LOGS", "PLOTS", EXPR_PLOTS, DBD_HISTORY, COOP_HISTORY, CONFIG_HISTORY, FITNESS_DIR]
    for d in dirs:
        if os.path.exists(ap.params["workspace"] + "/" + ap.params["runid"] + "/" + d):
            if ap.params["clearcache"]:
                """ Clear data from previous runs."""
                os.system("rm -rf " + ap.params["workspace"] + "/" + ap.params["runid"] + "/" + d)
        if False == os.path.exists(ap.params["workspace"] + "/" + ap.params["runid"] + "/" + d):
            """Build the directory"""
            os.system("mkdir " + ap.params["workspace"] + "/" + ap.params["runid"] + "/" + d)

def main():
    ap = ArgParser(sys.argv)
    read_cli(ap)

    if rank == 0:
        """Verify existence of necessary folders."""
        check_workspace(ap)

    l = []
    for s in range(0, ap.params["popsize"]):
        l.append( s )
    lmi = list_my_items(l, comm.Get_rank() )
    
    population = None
    landscape = None
    
    """Master builds a population and sends to slaves."""    
    if comm.Get_rank() == 0:
        """Build the Population..."""
        population = Population() # Creates an population data structure, but with no individuals or genomes
        if ap.params["verbosity"] > 1:
            print "\n. Building the population. . ."
        popath = ap.getOptionalArg("--pop_path")
        """If --pop_path is specified, it should be accompanied with --start_generation"""
        if popath != False:
            population.init_from_pickle(popath) # Initializes the population to be a saved population
        else:
            cli_genes = get_genes_from_file(ap)            
            population.init(ap, init_genes=cli_genes) # Initializes the population to be many clones of the user-spec. genes
        
        if ap.params["verbosity"] > 1:    
            print population.get_info()
        
        """Broadcast the updated population to MPI slaves."""
        pop_data = population.collapse()
        pop_data_pickle = pickle.dumps( pop_data )
        for slave in range(1, comm.Get_size()):
            comm.send(pop_data_pickle, dest=slave, tag=11)
    else:
        """Slaves receive the population from master."""
        pop_data_pickle = comm.recv(source=0, tag=11)
        pop_data = pickle.loads( pop_data_pickle ) 
        population = Population()
        population.uncollapse(pop_data)
        
    """Master and slaves build their own copies of the landscape, using the user specifications."""
    landscape = Landscape(ap)
    landscape.init(ap)

    if rank == 0:
        """Check for consistency with all parameters."""
        check_world_consistency(ap, population, landscape)    
    comm.Barrier()
        
    #
    # Default behavior is to run the genetic algorithm:
    #
    if ap.params["doko"] == False: # we're *not* doing a K.O. test
        ga = Genetic_Algorithm(ap)        
        ga.population = population
        ga.landscape = landscape
        
        if rank == 0 and ap.params["verbosity"] > 0:
            print "\n. Starting the Simulation. . .\n"
        comm.Barrier()
        
        ga.runsim(ap)

    #
    # If we're doing a KO test, then do that here...
    #
    if ap.params["doko"] == True:
        ko = KO(ap)
        ko.genome = population.genomes[ ap.params["kogenome"] ]
        ko.landscape = landscape
        if ap.params["verbosity"] > 1 and rank == 1:
            print "\n. Starting a knock-out experiment. . ."
            print "--> KO on each gene in individual " + ap.params["kogenome"].__str__()
            print "--> Assessing fitness at generation " + ap.params["generation"].__str__()
        ko.runko(ap)
        
    if rank == 0:
        return [population, landscape]
    else:
        return None

def splash():
    print "======================================="
    print "."
    print ". Sim Genome"
    print "."
    print ". " + VERSIONSTRING
    print "."
    print "======================================="
    
#
# main:
#    
rank = comm.Get_rank()
if rank == 0:
    splash()
mpi_check()

#
# uncomment this code to run Cprofile.
#
#if ap.getOptionalArg("--use_cprofile"):
#    prof_path = "./prof_trace." + (random.random() * 1000000).__str__() + ".cprofile"
#    cProfile.run( 'main()', prof_path )
#else:
main()