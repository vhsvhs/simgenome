from configuration import *
from gene import *

class Genome:
    id = None
    
    """genes 0 through N_TR-1 are TR genes, genes N_TR through N_REPORTER are reporter genes."""
    genes = []
    
    def __init__(self, id):
        self.id = id
        """Add N_TR number of transcription factor genes"""
        for i in range(0, N_TR):
            self.genes.append( Gene(self.generate_unique_geneid(), has_dbd=True) )
        """Add N_REPORTER number of transcription factor genes"""
        for i in range(0, N_REPORTER):
            self.genes.append( Gene(self.generate_unique_geneid(), has_dbd=False) )
    
    def contains_gene(self, id):
        """Does the genome contain a gene with ID = id?"""
        for gene in self.genes:
            if gene.id == id:
                return True
        return False
    
    def generate_unique_geneid(self):
        """Returns a unique integer ID for a new gene"""
        randid = random.randint(0,1000000)
        while (True == self.contains_gene(randid) ):
            randid = random.randint(0,1000000)
        return randid
    
    def get_expression_levels(self):