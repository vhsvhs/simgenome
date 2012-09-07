from argparser import *
from configuration import *
from genetic_algorithm import *
from population import *
from version import *
from cli import *

def check_workspace(ap):
    if False == os.path.exists(ap.getArg("--workspace")):
        os.system("mkdir " + ap.getArg("--workspace"))
    if False == os.path.exists(ap.getArg("--workspace") + "/" + ap.getArg("--runid")):
        os.system("mkdir " + ap.getArg("--workspace") + "/" + ap.getArg("--runid"))
    dirs = [POPPICKLES, "LOGS", "PLOTS", EXPR_PLOTS]
    for d in dirs:
        if False == os.path.exists(ap.getArg("--workspace") + "/" + ap.getArg("--runid") + "/" + d):
            os.system("mkdir " + ap.getArg("--workspace") + "/" + ap.getArg("--runid") + "/" + d)
    # clear expression histories from previous runs
    os.system("rm -rf " + ap.getArg("--workspace") + "/" + ap.getArg("--runid") + "/" + EXPR_PLOTS + "/*")

def main():
    ap = ArgParser(sys.argv)
    read_cli(ap)

    l = []
    for s in range(0, ap.params["popsize"]):
        l.append( s )
    lmi = list_my_items(l, comm.Get_rank() )
    
    """Build a population"""    
    population = Population()
    if comm.Get_rank() == 0 and ap.params["verbosity"] > 1:
        print "\n. Buildng the population. . ."
    popath = ap.getOptionalArg("--pop_path")
    """If --pop_path is specified, it should be accompanied with --start_generation"""
    if popath != False:
        population.init_from_pickle(popath)
    else:
        cli_genes = get_genes_from_file(ap)
        population.init(ap, init_genes=cli_genes)
    if comm.Get_rank() == 0 and ap.params["verbosity"] > 1:    
        print population.get_info()
                
                      
    """Build a fitness landscape"""
    landscape = Landscape(ap)
    landscape.init(ap, genome = population.genomes[0])
        
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