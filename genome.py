from configuration import *
from gene import *

class Genome:
    id = None
    
    """genes 0 through N_TR-1 are TR genes, genes N_TR through N_REPORTER+N_TR-1 are reporter genes."""
    genes = []
    
    def __init__(self, id, init_genes = None):
        p = ProgressBar()
        
        self.id = id
        if init_genes == None:
            """Add N_TR number of transcription factor genes"""
            for i in range(0, N_TR):
                self.genes.append( Gene(i, has_dbd=True) )
                p.render(i, 'step %s' % i)
            """Add N_REPORTER number of transcription factor genes"""
            for i in range(0, N_REPORTER):
                self.genes.append( Gene(N_TR + i, has_dbd=False) )
                p.render(i, 'step %s' % i)
        else:
            self.genes = list(init_genes)
    
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
        pass