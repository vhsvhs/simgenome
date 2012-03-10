from configuration import *
from gene import *

class Genome:
    id = None
    """genes 0 through N_TR-1 are TR genes, genes N_TR through N_REPORTER+N_TR-1 are reporter genes."""
    genes = []
    
    """Data computed by the last call to Landscape.get_fitness."""
    gene_expr = {}     
    
    def __init__(self, id):
        self.id = id
    
    def init(self, init_genes = None):
        """init_genes will be used for self.genes, unless it's None"""
        if init_genes == None:
            """Add N_TR number of transcription factor genes"""
            for i in range(0, N_TR):
                repressor = False
                if i%2:
                    repressor = True
                self.genes.append( Gene(i, has_dbd=True, repressor=repressor) )
            """Add N_REPORTER number of transcription factor genes"""
            for i in range(0, N_REPORTER):
                self.genes.append( Gene(N_TR + i, has_dbd=False) )
        else:
            self.genes = list(init_genes)
    
    def uncollapse(self, data):
        gids = data[0].keys()
        gids.sort()
        for gid in gids:
            this_gene = Gene(data[0][gid][0], data[0][gid][1], data[0][gid][2], data[0][gid][3], data[0][gid][4])
            self.genes.append( this_gene )
        self.gene_expr = data[1]
    
    def collapse(self):
        data = {}
        for g in self.genes:
            data[g.id] = g.collapse()
        return [data, self.gene_expr]
    
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
    
    def count_cis_seq_len(self):
        count = 0
        for g in self.genes:
            count += g.urs.__len__()
        return count