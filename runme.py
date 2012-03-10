from argparser import *
from configuration import *
from genetic_algorithm import *
from population import *
from version import *
from debug_tools import *

def main():
    ap = ArgParser(sys.argv)
    if rank == 0:
        """Verify consistency of command-line arguments"""
            
        """Build a population"""
        population = Population()
        population.init()
        
        """Build a random fitness landscape"""
        landscape = Landscape()
        landscape.init(genome = population.genomes[0])
            
        """Broadcast the population and landscape to slaves"""
        pop_data = population.collapse()
        pop_data_pickle = pickle.dumps( pop_data )
        land_data = landscape.collapse()
        land_data_pickle = pickle.dumps( land_data )
        for slave in range(1, comm.Get_size()):
            comm.send([pop_data_pickle, land_data_pickle], dest=slave, tag=11)
        #print "debug test:", rank, ":", population.genomes[0].genes[0].urs
    else:
        """Recieve the population and landscape from master"""
        [pop_data_pickle, land_data_pickle] = comm.recv(source=0, tag=11)
        pop_data = pickle.loads( pop_data_pickle ) 
        population = Population()
        population.uncollapse(pop_data)
        land_data = pickle.loads( land_data_pickle ) 
        landscape = Landscape()
        landscape.uncollapse(land_data)

        #print "debug test:", rank, ":", population.genomes[0].genes[0].urs
    
    """Setup the genetic algorithm, using the population and landscape"""
    ga = Genetic_Algorithm(ap)
    ga.population = population
    ga.landscape = landscape
        
    """Run the simulation."""
    ga.runsim()

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
print "\n. MPI process", rank, "of", comm.Get_size(), "is alive."
comm.Barrier()
x = main()