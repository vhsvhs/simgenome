from configuration import *
from genome import *

class Population:
    genomes = []
    
    def __init__(self):
        print "\n. Building a new population..."
        init_genes = None
        for i in range(0, N_GENOMES):
            """Fill the population with N_GENOMES copies of the seed genome."""
            genome = Genome( i, init_genes )
            if init_genes == None:
                init_genes = genome.genes
            self.genomes.append( genome )
    
    def contains_genome(self, id):
        """Does the population contain a genome with ID == id?"""
        for genome in self.genomes:
            if genome.id == id:
                return True
        return False
    
    def generate_unique_genomeid(self):
        """Returns a unique integer ID for a new genome"""
        randid = random.randint(0,1000000)
        while (True == self.contains_genome(randid) ):
            randid = random.randint(0,1000000)
        return randid
    
    def do_mutations(self):
        """Introduce mutations into (potentially) all genomes in the population"""
        n_point_mutations = None
        n_deletions = None
        n_duplications = None
        
    def do_reproduction(self):
        """reproduces a proportion of genomes based on their fitness"""
        pass
    

    
    