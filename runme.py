from argparser import *
from configuration import *
from genetic_algorithm import *
from population import *
from version import *
from cli import *

def check_workspace(ap):
    if False == os.path.exists(ap.params["workspace"]):
        os.system("mkdir " + ap.params["workspace"])
    if False == os.path.exists(ap.params["workspace"] + "/" + ap.params["runid"]):
        os.system("mkdir " + ap.params["workspace"] + "/" + ap.params["runid"])
    dirs = [POPPICKLES, "LOGS", "PLOTS", EXPR_PLOTS]
    for d in dirs:
        if os.path.exists(ap.params["workspace"] + "/" + ap.params["runid"] + "/" + d):
            # clear data from previous runs
            os.system("rm -rf " + ap.params["workspace"] + "/" + ap.params["runid"] + "/" + d)
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
        population = Population()
        if ap.params["verbosity"] > 1:
            print "\n. Buildng the population. . ."
        popath = ap.getOptionalArg("--pop_path")
        """If --pop_path is specified, it should be accompanied with --start_generation"""
        if popath != False:
            population.init_from_pickle(popath)
        else:
            cli_genes = get_genes_from_file(ap)
            population.init(ap, init_genes=cli_genes) # to-do: fix this! so that only one proc eventually calls gene.set_gamma (which builds random gamma distros0
        if ap.params["verbosity"] > 1:    
            print population.get_info()
        """Broadcast the updated population to MPI slaves."""
        pop_data = population.collapse()
        pop_data_pickle = pickle.dumps( pop_data )
        for slave in range(1, comm.Get_size()):
            comm.send(pop_data_pickle, dest=slave, tag=11)
    else:
        """Slaves recieve the init population from master."""
        pop_data_pickle = comm.recv(source=0, tag=11)
        pop_data = pickle.loads( pop_data_pickle ) 
        population = Population()
        population.uncollapse(pop_data)
        
        
    landscape = Landscape(ap)
    landscape.init(ap, genome = population.genomes[0])
#    """Master builds the landscape."""
#    if comm.Get_rank() == 0:
#        
#        land_data = landscape.collapse()
#        land_data_pickle = pickle.dumps( land_data )
#        for slave in range(1, comm.Get_size()):
#            comm.send(land_data_pickle, dest=slave, tag=11)
#    else:
#        land_data_pickle = comm.recv(source=0, tag=11)
#        land_data = pickle.loads( land_data_pickle ) 
#        landscape = Population()
#        landscape.uncollapse(land_data)

    
    ga = Genetic_Algorithm(ap)        
    ga.population = population
    ga.landscape = landscape
        
    if rank == 0:
        """Check for consistency with all parameters."""
        check_world_consistency(ap, population, landscape)
         
    comm.Barrier()
    
    """Setup the genetic algorithm, using the population and landscape"""
        
    """Run the simulation."""
    ga.runsim(ap)

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

#if ap.getOptionalArg("--use_cprofile"):
#   prof_path = "./prof_trace." + (random.random() * 1000000).__str__() + ".cprofile"
#    cProfile.run( 'main()', prof_path )
#else:
main()