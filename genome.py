from configuration import *
from gene import *

GLOBAL_GEN_COUNTER = 0
GLOBAL_T_COUNTER = 0

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
        self.is_elite = False
    
    def init(self, ap):
        """init_genes will be used for self.genes, unless it's None"""
        init_genes = self.get_genes_from_file(ap)
        if init_genes == None:
            #print "Genome", self.id, "has", self.genes.__len__(), "genes."
            """Add ap.params["numtr"] number of transcription factor genes"""
            for i in range(0, ap.params["numtr"]):
                repressor = False
                if i%2:
                    repressor = True
                self.genes.append( Gene(i, ap.params["init_urs_len"], has_dbd=True, repressor=repressor) )
            """Add N_REPORTER number of transcription factor genes"""
            for i in range(0, N_REPORTER):
                self.genes.append( Gene(ap.params["numtr"] + i, has_dbd=False) )
            #print "Genome", self.id, "has", self.genes.__len__(), "genes."
        else:
            self.genes = init_genes
    
    """Returns either a list of genes read from a file, OR returns None if the user
    did not specify to use genes from a file."""
    def get_genes_from_file(self, ap):
        #print "Getting genes from file..."
        this_pwm = None
        if ap.getOptionalArg("--pwmpath"):
            pwmpath = ap.getOptionalArg("--pwmpath")
            this_pwm = PWM()
            this_pwm.read_from_file( pwmpath )
        if ap.getOptionalArg("--genepath"): 
            ret_genes = []
            genepath = ap.getOptionalArg("--genepath")
            fin = open(genepath, "r")
            #print "\n. Reading genes from", genepath
            for l in fin.readlines():
                if l.startswith("#"):
                    continue
                else:
                    tokens = l.split()
                    if tokens.__len__() >= 4:
                        this_id = int(tokens[0])
                        this_has_dbd = int(tokens[1])
                        this_repressor = int(tokens[2])
                        if this_repressor == 0:
                            this_repressor = False
                        elif this_repressor == 1:
                            this_repressor = True
                        this_urs = tokens[3]
                        this_gene = Gene(this_id, this_urs.__len__(), urs=this_urs, has_dbd=this_has_dbd, repressor=this_repressor, pwm=this_pwm) 
                        ret_genes.append(this_gene)
                        #print "\n.", this_id, this_has_dbd, this_repressor, this_urs
            fin.close()
            count_tf = 0
            for gid in ret_genes:
                if gid.has_dbd:
                    count_tf += 1
            ap.params["numtr"] = count_tf
            return ret_genes
        else:
            return None
        
    
    def uncollapse(self, data):
        gids = data[0].keys()
        gids.sort()
        for gid in gids:
            this_gene = Gene(data[0][gid][0], data[0][gid][1], data[0][gid][2], data[0][gid][3], data[0][gid][4], data[0][gid][5])
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