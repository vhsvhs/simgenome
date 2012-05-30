from argparser import *
from configuration import *
from genetic_algorithm import *
from population import *
from version import *
from debug_tools import *
from cli import *

def check_workspace(ap):
    if False == os.path.exists(ap.getArg("--runid")):
        os.system("mkdir " + ap.getArg("--runid"))
    dirs = [POPPICKLES, "LOGS", "PLOTS", EXPR_PLOTS]
    for d in dirs:
        if False == os.path.exists(ap.getArg("--runid") + "/" + d):
            os.system("mkdir " + ap.getArg("--runid") + "/" + d)

def main():
    ap = ArgParser(sys.argv)
    read_cli(ap)

    l = []
    for s in range(0, ap.params["popsize"]):
        l.append( s )
    lmi = list_my_items(l, comm.Get_rank() )
    
    # for debugging:
    # print "runme.py 27, proc.", comm.Get_rank(), "-", lmi.__len__(), "items:", lmi 

    """Build a population"""    
    population = Population()
    popath = ap.getOptionalArg("--pop_path")
    """--pop_path should always be specified with --start_generation"""
    if popath != False:
        population.init_from_pickle(popath)
    else:
        cli_genes = get_genes_from_file(ap)
        population.init(ap, init_genes=cli_genes)
                      
    """Build a fitness landscape"""
    landscape = Landscape(ap)
    cli_timepatterns = get_timepatterns_from_file(ap)
    landscape.init(ap, genome = population.genomes[0], tp=cli_timepatterns)
        
    if rank == 0:
        """Verify existence of necessary folders."""
        check_workspace(ap)

        """Check for consistency with all parameters."""
        check_world_consistency(ap, population, landscape)
                 
    comm.Barrier()
    
    """Setup the genetic algorithm, using the population and landscape"""
    ga = Genetic_Algorithm(ap)
    ga.population = population
    ga.landscape = landscape
        
    """Run the simulation."""
    ga.runsim(ap)

    if rank == 0:
        return [ga.gen_gid_fitness, population, landscape]
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
x = main()