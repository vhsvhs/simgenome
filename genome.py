from configuration import *
from gene import *

class Genome:
    id = None
    """genes 0 through ap.params["numtr"]-1 are TR genes, genes ap.params["numtr"] through N_REPORTER+N_TR-1 are reporter genes."""
    genes = []
    
    """Data computed by the last call to Landscape.get_fitness."""
    gene_expr = {}     
    
    def __init__(self, id):
        self.id = id
        self.genes = []
        self.gene_expr = {}
        self.tfcoop = None # tfcoop gets transformed to create gamma
        self.gamma = None # 3-D Numpy array, gamma[gene id, gene id, distance] = cooperative/competitive factor.  1.0 = no effect on binding, <1.0 = decreases binding, >1.0 = strengthens binding.
        self.is_elite = False
    
    def init(self, ap, init_genes=None, init_expression=None):
        if init_genes == None:
            """Add ap.params["numtr"] number of transcription factor genes"""
            for i in ap.params["rangetrs"]:
                repressor = False
                if i%2:
                    repressor = True
                self.genes.append( Gene(i, ap.params["init_urs_len"], has_dbd=True, repressor=repressor) )
            """Add N_REPORTER number of transcription factor genes"""
            for i in range(0, ap.params["numreporter"]):
                self.genes.append( Gene(ap.params["numtr"] + i, ap.params["init_urs_len"], has_dbd=False) )
        else:
            self.genes = init_genes
        
        if ap.params["enable_epigenetics"] == True and init_expression != None:
            self.gene_expr = init_expression
        
        # init the TF coop matrix:
        if ap.params["coopinit"] == "random":
            self.tfcoop = numpy.random.gamma(2.0,10.0, (ap.params["numtr"], ap.params["numtr"])) - 4.0
            # to-do: ensure that values in tfcoop range form -1 to +infinity
        else:
            self.tfcoop = zeros( (ap.params["numtr"], ap.params["numtr"]), dtype=float)
        self.set_gamma(ap)
        
    def uncollapse(self, data):
        gids = data[0].keys()
        gids.sort()
        for gid in gids:
            this_gene = Gene(data[0][gid][0], data[0][gid][1], data[0][gid][2], data[0][gid][3], data[0][gid][4], data[0][gid][5])
            self.genes.append( this_gene )
        self.gene_expr = data[1]
        self.tfcoop = data[2]
        self.gamma = data[3]
    
    def collapse(self):
        data = {}
        for g in self.genes:
            data[g.id] = g.collapse()
        return [data, self.gene_expr, self.tfcoop, self.gamma]
    
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

    def coopfunc(self, g, d):
        """Calculates the degree of binding cooperativity between two TFs binding distance d apart.
        g is positive for synergistic interactions, and negative for antagonistic interactions.
        g ranges from -1 to +inifinity. g of zero means no effect on binding."""
        return 1 + g * math.exp( (-1)*(d**2)/V_RATE_OF_COOP_DECAY );

    def set_gamma(self, ap):
        """Prec-calculates the cooperative/competitive binding interactions between all TFs."""
        self.gamma = zeros( (ap.params["numtr"], ap.params["numtr"], ap.params["maxgd"]), dtype=float)
        for i in ap.params["rangetrs"]:
            for j in ap.params["rangetrs"]:
                for d in ap.params["rangegd"]:
                    #print i, j, d
                    self.gamma[i,j,d] = self.coopfunc( self.tfcoop[i,j], d)
        print "tfcoop, genome", self.id 
        print self.tfcoop
        print "gamma, genome", self.id
        print self.gamma