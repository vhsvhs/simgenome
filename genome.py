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
        self.is_elite = False
    
    def init(self, ap, init_genes=None, init_expression=None):
        if init_genes == None:
            """Add ap.params["numtr"] number of transcription factor genes"""
            for i in ap.params["trlist"]:
                repressor = False
                if i%2:
                    repressor = True
                self.genes.append( Gene(i, ap.params["init_urs_len"], has_dbd=True, repressor=repressor, ap=ap) )
            """Add N_REPORTER number of transcription factor genes"""
            for i in range(0, ap.params["numreporter"]):
                self.genes.append( Gene(ap.params["numtr"] + i, ap.params["init_urs_len"], has_dbd=False, ap=ap) )
        else:
            self.genes = init_genes
        
        if ap.params["enable_epigenetics"] == True and init_expression != None:
            self.gene_expr = init_expression
        
        # init the TF coop matrix:
        # 10/4 - moved to gene.py
#        if ap.params["coopinit"] == "random":
#            self.tfcoop = numpy.random.gamma(2.0,10.0, (ap.params["numtr"], ap.params["numtr"])) - 4.0
#            # to-do: ensure that values in tfcoop range form -1 to +infinity
#        else:
#            self.tfcoop = zeros( (ap.params["numtr"], ap.params["numtr"]), dtype=float)
#        self.set_gamma(ap)
        
        self.print_dbds(ap)
        # print the gamma matrix
        self.print_gamma(ap)

        
    def uncollapse(self, data):
        gids = data[0].keys()
        gids.sort()
        for gid in gids:
            this_gene = Gene(data[0][gid][0], data[0][gid][1], urs=data[0][gid][2], has_dbd=data[0][gid][3], repressor=data[0][gid][4], pwm=data[0][gid][5], tfcoop=data[0][gid][6], gamma=data[0][gid][7], name=data[0][gid][8])
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

    
    def count_cis_seq_len(self):
        count = 0
        for g in self.genes:
            count += g.urs.__len__()
        return count

    
    def print_gamma(self, ap):
        foutpath = ap.params["workspace"] + "/" + ap.params["runid"] + "/" + EXPR_PLOTS + "/coop.gen" + ap.params["generation"].__str__() + ".gid" + self.id.__str__() + ".txt"
        lines = []
        lines.append("dist.\tTFi\tTFj\tgamma")
        for gene in self.genes:
            if gene.has_dbd:
                for j in ap.params["trlist"]:                    
                    for d in ap.params["rangegd"]:
                        this_line = d.__str__() + "\t" + gene.id.__str__() + "\t" + j.__str__() 
                        this_line += "\t" + gene.gamma[j,d].__str__()
                        lines.append(this_line)
        fout = open(foutpath, "w")
        for l in lines:
            fout.write(l + "\n")
        fout.close()

    def print_dbds(self, ap):
        foutpath = ap.params["workspace"] + "/" + ap.params["runid"] + "/" + EXPR_PLOTS + "/dbds.gen" + ap.params["generation"].__str__() + ".gid" + self.id.__str__() + ".txt"
        lines = []
        for gene in self.genes:
            if gene.has_dbd:
                lines.append("\nGene " + gene.id.__str__() + ":\n")
                lines.append( gene.pwm.__str__())
        fout = open(foutpath, "w")
        for l in lines:
            fout.write(l + "\n")
        fout.close()

        