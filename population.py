from configuration import *
from genome import *

class Population:
    genomes = {}
    
    def __init__(self):
        self.genomes = {}
    
    def init(self, ap, init_genes=None):
        if init_genes == None:
            if comm.Get_rank() == 0:
                print "\n. Building a random population..."
        else:
            if comm.Get_rank() == 0:
                print "\n. Creating a population, specified by", ap.getOptionalArg("--urspath"), "and",ap.getOptionalArg("--pwmpath") 
        for i in range(0, ap.params["popsize"]):
            """Fill the population with ap.params["popsize"] copies of the seed genome."""
            #print "+ Genome ", i
            self.genomes[ i ] = Genome(i)
            self.genomes[ i ].init( ap, init_genes=init_genes)
    
    def init_from_pickle(self, picklepath):
        fin = open(picklepath, "r")
        pop_data = pickle.load( fin )
        fin.close()
        self.uncollapse(pop_data) 
    
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
    
    def effective_popsize(self):
        return self.genomes.__len__()
    
    def mark_elite(self, fitness_gid, max_fitness, min_fitness, ap):
        if max_fitness == min_fitness:
            return
        """Mark the elite individuals."""
        count_elite = 0
        max_elite = self.genomes.__len__() * ap.params["eliteproportion"]
        fkeys = fitness_gid.keys()
        fkeys.sort()
        fkeys.reverse()
        for fitness in fkeys:
            for gid in fitness_gid[fitness]:
                if count_elite < max_elite:
                    self.genomes[gid].is_elite = True
                    count_elite += 1
                else:
                    self.genomes[gid].is_elite = False
    
    def print_fitness(self, ap, gid_fitness):
        if int(ap.getOptionalArg("--verbose")) > 2:
            print "\n................................................\n"
            for gid in gid_fitness:
                marka = ""
                if self.genomes[gid].is_elite:
                    marka = "(elite)"
                print "\t. ID", gid, "f= %.3f"%gid_fitness[gid], marka
            print "\n................................................"
    
    def do_mutations(self, ap):
        if int(ap.getOptionalArg("--verbose")) > 2:
            print "\n. The population is mutating. . .\n"
            
        """Introduce mutations into (potentially) all genomes in the population"""
        for gid in self.genomes:
            mu = ap.params["mu"] 
            if self.genomes[gid].is_elite == True:
                mu = ap.params["elitemu"] 
                
            #print "debug population.py 92 - ID", gid, "mu", mu
            
            if mu > 0.0:
                """How many mutations should we make?"""
                n_point_mutations = int(self.genomes[gid].count_cis_seq_len() * mu)
                n_deletions = None
                n_duplications = None
                if int(ap.getOptionalArg("--verbose")) >= 2:                
                    print "\t.", n_point_mutations, "cis mutations to individual", gid
                """URS mutations...."""
                for i in range(0, n_point_mutations):
                    """Pick a random gene"""
                    rand_gene = random.randint(0, self.genomes[gid].genes.__len__()-1)
                    """Pick a random site"""
                    rand_site = random.randint(0, self.genomes[gid].genes[rand_gene].urs.__len__()-1)
                    """Mutate!"""
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
                """PWM mutations..."""
                rand_roll = random.random()
                if rand_roll < mu:
                    rand_tr_id = random.randint(0, ap.params["numtr"]-1)
                    if int(ap.getOptionalArg("--verbose")) >= 2:
                        print "\t+ mutating PWM ", rand_tr_id, "in individual", gid
                    self.genomes[gid].genes[rand_tr_id].pwm.mutate(ap)
        if int(ap.getOptionalArg("--verbose")) > 2:
            print "\n"
                
                    
    def fitness_cdf_sampler(self, min, max, sum, gid_fitness):
        if sum == 0:
            """In this cases, all individuals have the same fitness. So pick one randomly"""
            x = random.sample(gid_fitness, 1)[0]
            return x
        else:
            randp = random.random() * sum
            gids = gid_fitness.keys()
            gids.sort()
            gids.reverse()
            randpsum = 0.0
            for i in gids:
                randpsum += gid_fitness[i] - min
                if randpsum > randp:
                    return i


    def get_minmax_fitness(self, gid_fitness):
        """First, find max and min fitness."""
        min_fitness = 1.0
        max_fitness = 0.0
        reverse_hash = {} # key = fitness, value = array of 1 or more IDs.
        for gid in gid_fitness:
            this_fitness = gid_fitness[gid]
            if reverse_hash.__contains__(this_fitness):
                reverse_hash[this_fitness].append( gid )
            else:
                reverse_hash[this_fitness] = [gid]     
            if gid_fitness[gid] < min_fitness:
                min_fitness = gid_fitness[gid]
            elif gid_fitness[gid] > max_fitness:
                max_fitness = gid_fitness[gid]
        sum_fitness = 0.0
        for gid in gid_fitness:
            sum_fitness += (gid_fitness[gid] - min_fitness)
        return [min_fitness, max_fitness, sum_fitness, reverse_hash]
        
    def do_reproduction(self, gid_fitness, min_fitness, max_fitness, sum_fitness, ap):        
        if int(ap.getOptionalArg("--verbose")) > 2:
            print "\n. The population is reproducing, selectively based on fitness. . .\n"
        
        new_genomes = {}
        """For now we assume the new population has the same size as its parent population."""
        
        for child_gid in gid_fitness:
            new_genomes[child_gid] = Genome(child_gid)                         
            if self.genomes[child_gid].is_elite == True:
                if int(ap.getOptionalArg("--verbose")) > 2:
                    print "\t+ new child", child_gid, "=", child_gid, "cloned."
                new_genomes[child_gid].is_elite = True
                new_genomes[child_gid].init(ap, init_genes = self.genomes[child_gid].genes)
            else:
                """Select the parents from the fitness CDF."""
                parent1 = self.fitness_cdf_sampler(min_fitness, max_fitness, sum_fitness, gid_fitness) 
                parent2 = self.fitness_cdf_sampler(min_fitness, max_fitness, sum_fitness, gid_fitness)
                if int(ap.getOptionalArg("--verbose")) > 2:
                    print "\t+ new child", child_gid, "=", parent1, "X", parent2
                new_genomes[child_gid].is_elite = False
                new_genes = []
                for geneid in range(0, self.genomes[parent1].genes.__len__()):
                    """Flip a coin; heads the child gets parent1's copy of this gene, tails the child gets parent2's copy."""
                    if random.randint(0,1):
                        parentgene = self.genomes[parent1].genes[geneid]
                    else:
                        parentgene = self.genomes[parent2].genes[geneid]
                    copypwm = None
                    if parentgene.has_dbd:
                        copypwm = PWM(copyfrom=parentgene.pwm)
                    gene_copy = Gene(geneid, parentgene.urs.__len__(), urs=parentgene.urs, has_dbd=parentgene.has_dbd, repressor=parentgene.is_repressor,pwm=copypwm)
                    new_genes.append( gene_copy )
                new_genomes[child_gid].init(ap, init_genes=new_genes)
            
        """Save the new children genomes."""
        self.genomes = new_genomes
            
    def compare_two_genomes(self, idx, idy):
        for j in range(0, self.genomes[idx].genes.__len__()):
            if self.genomes[idx].genes[j].urs.__contains__( self.genomes[idy].genes[j].urs ):
                print "Genomes", idx, "and", idy, "differ in URS on gene", j
    