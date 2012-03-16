from configuration import *
from genome import *

class Population:
    genomes = {}
    
    def __init__(self):
        pass
    
    def init(self, ap):
        print "\n. Building a new population..."
        for i in range(0, N_GENOMES):
            """Fill the population with N_GENOMES copies of the seed genome."""
            genome = Genome(i)
            print "\n+ Genome ", i
            genome.init( ap )
            self.genomes[ i ] = genome
    
    def uncollapse(self, data):
        for gid in data:
            this_genome = Genome(gid)
            this_genome.uncollapse(data[gid])
            self.genomes[gid] = this_genome
    
    def collapse(self):
        data = {}
        for gid in self.genomes.keys():
            data[gid] = self.genomes[gid].collapse()
        return data
    
    def list_genome_ids(self):
        l = self.genomes.keys()
        l.sort()
        return l
            
    def contains_genome(self, id):
        """Does the population contain a genome with ID == id?"""
        for gids in self.genomes.keys():
            if genome.id == gid:
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
        for gid in self.genomes.keys():
            n_point_mutations = int(self.genomes[gid].count_cis_seq_len() * MU)
            n_deletions = None
            n_duplications = None
            
            #print "\n. Introducing", n_point_mutations, "point mutations to individual", g.id
            
            for i in range(0, n_point_mutations):
                # pick a random gene
                rand_gene = random.randint(0, self.genomes[gid].genes.__len__()-1)
                # pick a random site
                rand_site = random.randint(0, self.genomes[gid].genes[rand_gene].urs.__len__()-1)
                # mutate!
                curr_state = self.genomes[gid].genes[rand_gene].urs[rand_site]
                ALPHABET.remove(curr_state)
                new_state = random.choice( ALPHABET )
                ALPHABET.append(curr_state)
                new_urs = ""
                for j in range(0, self.genomes[gid].genes[rand_gene].urs.__len__()):
                    if j == rand_site:
                        new_urs += new_state.__str__()
                    else:
                        new_urs += self.genomes[gid].genes[rand_gene].urs[j]
                self.genomes[gid].genes[rand_gene].urs = new_urs
        
    def do_reproduction(self):
        """reproduces a proportion of genomes based on their fitness"""
        pass
    
    def compare_two_genomes(self, idx, idy):
        for j in range(0, self.genomes[idx].genes.__len__()):
            if self.genomes[idx].genes[j].urs.__contains__( self.genomes[idy].genes[j].urs ):
                print "Genomes", idx, "and", idy, "differ in URS on gene", j
    

    
    