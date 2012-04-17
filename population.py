from configuration import *
from genome import *

class Population:
    genomes = {}
    
    def __init__(self):
        self.genomes = {}
    
    def init(self, ap):
        print "\n. Building a new population..."
        for i in range(0, ap.params["popsize"]):
            """Fill the population with ap.params["popsize"] copies of the seed genome."""
            self.genomes[ i ] = Genome(i)
            print "+ Genome ", i
            self.genomes[ i ].init( ap )
    
    def uncollapse(self, data):
        for gid in data:
            this_genome = Genome(gid)
            this_genome.uncollapse(data[gid]    )
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
        max_elite = self.genomes.__len__() * 1.0 / ELITE_PROPORTION
        fkeys = fitness_gid.keys()
        fkeys.sort()
        for fitness in fkeys:
            for gid in fitness_gid[fitness]:
                if count_elite < max_elite:
                    self.genomes[gid].is_elite = True
    
    def do_mutations(self, ap):
        if int(ap.getOptionalArg("--verbose")) > 2:
            print "\n. The population is mutating. . .\n"
            

        """Introduce mutations into (potentially) all genomes in the population"""
        for gid in self.genomes.keys():
            mu = ap.params["mu"] 
            if self.genomes[gid].is_elite:
                mu = ap.params["elitemu"] 
            
            n_point_mutations = int(self.genomes[gid].count_cis_seq_len() * mu)
            n_deletions = None
            n_duplications = None
            
            if int(ap.getOptionalArg("--verbose")) > 2:                
                print "\t", n_point_mutations, "point mutations to individual", gid
            
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
                #print new_urs
        
            rand_roll = random.random()
            if rand_roll < mu:
                rand_tr_id = random.randint(0, ap.params["numtr"]-1)
                print "mutating PWM ", rand_tr_id
                print "old PWM=", self.genomes[gid].genes[rand_tr_id].pwm.P.__str__()
                self.genomes[gid].genes[rand_tr_id].pwm.mutate(ap)
        if int(ap.getOptionalArg("--verbose")) > 2:
            print "\n"
                
                    
    def fitness_cdf_sampler(self, min, max, sum, gid_fitness):
        if sum == 0:
            x = random.sample(gid_fitness, 1)[0]
            #print "CDF sample: random", x, "from", max
            return x
        else:
            randp = random.random() * sum
            gids = gid_fitness.keys()
            gids.sort()
            randpsum = 0.0
            for i in gids:
                randpsum += gid_fitness[i] - min
                if randpsum > randp:
                    print "CDF sample", i, "from", max
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
            print "\n. The population is reproducing selectively based on fitness. . .\n"
            for gid in gid_fitness:
                marka = ""
                if self.genomes[gid].is_elite:
                    marka = "+"
                print "\tID", gid, "f= %.3f"%gid_fitness[gid], marka
            #print gid_fitness
        
        new_genomes = {}
        for child_gid in gid_fitness:
            """Here the new population has the same size as its parent population."""
            new_genomes[child_gid] = Genome(child_gid) 
            
            """Select two parents, by selecting from the fitness CDF."""
            parent1 = self.fitness_cdf_sampler(min_fitness, max_fitness, sum_fitness, gid_fitness)
            parent2 = self.fitness_cdf_sampler(min_fitness, max_fitness, sum_fitness, gid_fitness)
            if int(ap.getOptionalArg("--verbose")) > 2:
                print "\tchild", child_gid, ":", parent1, "X", parent2
                
            """Cross the parents"""
            for geneid in range(0, self.genomes[parent1].genes.__len__()):
                """Flip a coin; heads the child gets parent1's copy of this gene, tails the child gets parent2's copy."""
                gene_copy = None
                if random.randint(0,1):
                    gene_copy = self.genomes[parent1].genes[geneid]
                else:
                    gene_copy = self.genomes[parent2].genes[geneid]
                new_genomes[child_gid].genes.append( gene_copy )
            
        """Save the new children genomes."""
        self.genomes = new_genomes
            
    def compare_two_genomes(self, idx, idy):
        for j in range(0, self.genomes[idx].genes.__len__()):
            if self.genomes[idx].genes[j].urs.__contains__( self.genomes[idy].genes[j].urs ):
                print "Genomes", idx, "and", idy, "differ in URS on gene", j
    

    
    