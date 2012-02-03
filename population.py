from configuration import *

class Population:
    self.genomes = []
    
    def __init__(self):
        """Build a seed genome."""
        genome = Genome()
        for i in range(0, N_GENOMES):
            self.genomes.append( Genome() )
    
    def do_mutations(self):
        """Introduce mutations into (potentially) all genomes in the population"""
        n_point_mutations = None
        n_deletions = None
        n_duplications = None
        
    def do_reproduction(self):
        """reproduces a proportion of genomes based on their fitness"""
        pass
    

    
    