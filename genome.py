from configuration import *

class Genome:
    self.id = None
    self.genes = []
    
    def __init__(self):
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
        """Returns a unique integer ID for a gene"""
        randid = random.randint(0,1000000)
        while (True == self.contains_gene(randid) ):
            randid = random.randint(0,1000000)
        return randid