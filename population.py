from configuration import *
from genome import *

class Population:
    genomes = []
    
    def __init__(self):
        pass
    
    def init(self):
        print "\n. Building a new population..."
        init_genes = None
        for i in range(0, N_GENOMES):
            """Fill the population with N_GENOMES copies of the seed genome."""
            genome = Genome(i)
            genome.init( init_genes )
            if init_genes == None:
                init_genes = genome.genes
            self.genomes.append( genome )
    
    def uncollapse(self, data):
        for gid in data:
            this_genome = Genome(gid)
            this_genome.uncollapse(data[gid])
            self.genomes.append( this_genome )
    
    def collapse(self):
        data = {}
        for g in self.genomes:
            data[g.id] = g.collapse()
        return data
        
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
        for g in self.genomes:
            n_point_mutations = int(g.count_cis_seq_len() * MU)
            n_deletions = None
            n_duplications = None
            
            #print "\n. Introducing", n_point_mutations, "point mutations to individual", g.id
            
            for i in range(0, n_point_mutations):
                # pick a random gene
                rand_gene = random.randint(0, g.genes.__len__()-1)
                # pick a random site
                rand_site = random.randint(0, g.genes[rand_gene].urs.__len__()-1)
                # mutate!
                curr_state = g.genes[rand_gene].urs[rand_site]
                ALPHABET.remove(curr_state)
                new_state = random.choice( ALPHABET )
                ALPHABET.append(curr_state)
                new_urs = ""
                for j in range(0, g.genes[rand_gene].urs.__len__()):
                    if j == rand_site:
                        new_urs += new_state.__str__()
                    else:
                        new_urs += g.genes[rand_gene].urs[j]
                self.genomes[g.id].genes[rand_gene].urs = new_urs
        
    def do_reproduction(self):
        """reproduces a proportion of genomes based on their fitness"""
        pass
    
    def compare_two_genomes(self, idx, idy):
        for j in range(0, self.genomes[idx].genes.__len__()):
            if self.genomes[idx].genes[j].urs.__contains__( self.genomes[idy].genes[j].urs ):
                print "Genomes", idx, "and", idy, "differ in URS on gene", j
    

    
    