from configuration import *
from genome import *

class Population:
    genomes = {}
    
    def __init__(self):
        self.genomes = {}
    
    def init(self, ap, init_genes=None): 
        for i in range(0, ap.params["popsize"]):
            """Fill the population with ap.params["popsize"] copies of the seed genome."""
            self.genomes[ i ] = Genome(i)
            self.genomes[ i ].init( ap, init_genes=init_genes)
    
    def init_from_pickle(self, picklepath):
        if comm.Get_rank() == 0:
            print "\n. Loading the pickled population at", picklepath
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
        
    def get_info(self):
        """Returns a multi-line string with basic information about the population."""
        ret = "--> Population size: "
        ret += (self.genomes.__len__()).__str__()
        ret += "\n--> Genes per individual: "
        ret += (self.genomes[0].genes.__len__()).__str__()
        return ret
        
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
        if ap.params["verbosity"] > 2:
            print "\n................................................\n"
            for gid in gid_fitness:
                marka = ""
                if self.genomes[gid].is_elite:
                    marka = "(elite)"
                print "\t. ID", gid, "f= %.3f"%gid_fitness[gid], marka
            print "\n................................................"
    
    def do_mutations(self, ap):
        if ap.params["verbosity"] > 2:
            print "\n. The population is mutating. . .\n"

        """Print LOGS/mu.genX.txt with the number of mutations made to each individual at this generation."""            
        fout = open(ap.params["workspace"] + "/" + ap.params["runid"] + "/LOGS/mu.gen" + ap.params["generation"].__str__() + ".txt", "a")
        fout.write("ID\tN_cis\tN_urslen\tN_PWM\tN_p2p\n")


        """First, we build some relevant probability distributions."""
        from scipy import stats
        ncis_norm = stats.norm(loc=ap.params["cismu"], scale=0.0)
        ncisindel_norm = stats.norm(loc=ap.params["urslenmu"], scale=NORM_SD)
        ndbd_norm = stats.norm(loc=ap.params["dbdmu"], scale=NORM_SD)
        ndbdindel_norm = stats.norm(loc=ap.params["pwmlenmu"], scale=NORM_SD)
        np2p_norm = stats.norm(loc=ap.params["p2pmu"], scale=NORM_SD)
    
        """Introduce mutations into (potentially) all genomes in the population"""
        for gid in self.genomes:
            mu = ap.params["mu"] 
            if self.genomes[gid].is_elite == True:
                mu = ap.params["elitemu"] 
            
            n_point_mutations = 0
            n_cis_indels = 0
            n_dbd_mutations = 0
            n_dbd_indels = 0
            n_p2p_changes = 0
            
            if mu > 0: # if --mu 1, basically.  --mu enables/disables all mutations.
                """cis nt mutations"""
                n_point_mutations = int(self.genomes[gid].count_cis_seq_len() * ncis_norm.rvs() )
                if n_point_mutations < 0:
                    n_point_mutations = 0
                if ap.params["verbosity"] >= 2:                
                    print "\t> ID", gid, ":", n_point_mutations, "random cis sequence substitutions."
                for i in range(0, n_point_mutations):
                    rand_gene = random.randint(0, self.genomes[gid].genes.__len__()-1)
                    self.genomes[gid].genes[rand_gene].mutate_urs()
                
                """cis indels"""
                n_cis_indels = int(self.genomes[gid].genes.__len__() * ncis_norm.rvs() )
                for i in range(0, n_cis_indels):
                    rand_gene = random.randint(0, self.genomes[gid].genes.__len__()-1)
                    if ap.params["verbosity"] >= 2:
                        print "\t> ID", gid, ": mutating length of regulatory region in gene", rand_gene
                    self.genomes[gid].genes[rand_gene].mutate_urs_len()
                
                """DBD mutations..."""
                n_dbd_mutations = int(ap.params["numtr"] * ndbd_norm.rvs() )
                for  i in range(n_dbd_mutations):
                    rand_tr_id = random.randint(0, ap.params["numtr"]-1)
                    if ap.params["verbosity"] >= 2:
                        print "\t> ID", gid, ": mutating DNA-binding specificity of TF", rand_tr_id
                    self.genomes[gid].genes[rand_tr_id].pwm.mutate(ap)
                
                """DBD indels"""
                n_dbd_indels = int(ap.params["numtr"] * ndbdindel_norm.rvs() )
                for i in range(n_dbd_mutations):
                    rand_tr_id = random.randint(0, ap.params["numtr"]-1)
                    did = self.genomes[gid].genes[rand_tr_id].pwm.mutate_len(ap)               
                    if did == True:
                        if ap.params["verbosity"] >= 2:
                            print "\t> ID", gid, ": mutating length of PWM for TF", rand_tr_id
                            
                
                """Gamma mutations...."""
                n_p2p_changes = int(ap.params["numtr"] * np2p_norm.rvs() )
                for i in range(n_p2p_changes):
                    rand_tr_id = random.randint(0, ap.params["numtr"]-1)
                    if ap.params["verbosity"] >= 2:
                        print "\t> ID", gid, ": mutating cofactor interactions for TR ", rand_tr_id
                    self.genomes[gid].genes[rand_tr_id].mutate_gamma(ap)
            
            fout.write(gid.__str__() + ("\t%d"%n_point_mutations).__str__() + ("\t%d"%n_cis_indels).__str__() + ("\t%d"%n_dbd_mutations).__str__() + ("\t%d"%n_dbd_indels).__str__() + ("\t%d"%n_p2p_changes).__str__() + "\n" )
        
        """Close the LOGS/mu.genX.txt"""
        fout.close()

        if ap.params["verbosity"] > 2:
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
        if ap.params["verbosity"] > 2:
            print "\n. The population is reproducing, selectively based on fitness. . .\n"
        
        foutpath = ap.params["workspace"] + "/" + ap.params["runid"] + "/LOGS/mating.gen" + ap.params["generation"].__str__() + ".txt"
        lines = []

        new_genomes = {}
        """For now we assume the new population has the same size as its parent population."""
        
        for child_gid in gid_fitness:
            new_genomes[child_gid] = Genome(child_gid)                         
            if self.genomes[child_gid].is_elite == True:
                if ap.params["verbosity"] > 2:
                    l = "\t> The elite child #" + child_gid.__str__() + " is cloned."# + child_gid.__str__() + " cloned."
                    print l
                    lines.append(l)
                new_genomes[child_gid].is_elite = True
                new_genomes[child_gid].init(ap, init_genes = self.genomes[child_gid].genes, init_expression=self.genomes[child_gid].gene_expr)
            else:
                """Select the parents from the fitness CDF."""
                parent1 = self.fitness_cdf_sampler(min_fitness, max_fitness, sum_fitness, gid_fitness) 
                
                do_sex = False
                if random.random() < ap.params["sexual_ratio"]:
                    do_sex = True
                    
                if do_sex:
                    """sexual reproduction with two parents."""
                    parent2 = self.fitness_cdf_sampler(min_fitness, max_fitness, sum_fitness, gid_fitness)
                else:
                    """clonal reproduction with one parent."""
                    parent2 = parent1
                if ap.params["verbosity"] > 2:
                    if do_sex:
                        sextoken = "cross"
                    else:
                        sextoken = "clone"
                    l = "\t> child " + sextoken + " " + child_gid.__str__() + " = " + parent1.__str__() + " X " + parent2.__str__()
                    print l
                    lines.append(l)
                new_genomes[child_gid].is_elite = False
                new_genes = []
                for geneid in range(0, self.genomes[parent1].genes.__len__()):
                    """Flip a coin; heads the child gets parent1's copy of this gene, tails the child gets parent2's copy."""
                    if random.randint(0,1):
                        parentgene = self.genomes[parent1].genes[geneid]
                    else:
                        parentgene = self.genomes[parent2].genes[geneid]
                    copypwm = None
                    copygamma = None
                    if parentgene.has_dbd:
                        copypwm = PWM(copyfrom=parentgene.pwm)
                        copygamma = numpy.copy( parentgene.gamma )
                        copytfcoop = numpy.copy( parentgene.tfcoop )
                        #print "debug population.py 215 parent gamma:\n", parentgene.gamma, "\n child gamma:\n", copygamma
                    gene_copy = Gene(geneid, parentgene.urs.__len__(), urs=parentgene.urs, has_dbd=parentgene.has_dbd, repressor=parentgene.is_repressor,pwm=copypwm,gamma=copygamma,tfcoop=copytfcoop)
                    new_genes.append( gene_copy )
                new_genomes[child_gid].init(ap, init_genes=new_genes, init_expression=self.genomes[parent1].gene_expr)
        """Save the new children genomes."""
        self.genomes = new_genomes
            
        fout = open(foutpath, "w")
        fout.write("\n. Generation " + ap.params["generation"].__str__() + "\n")
        for l in lines:
            fout.write(l + "\n")
        fout.close()
            
    def compare_two_genomes(self, idx, idy):
        for j in range(0, self.genomes[idx].genes.__len__()):
            if self.genomes[idx].genes[j].urs.__contains__( self.genomes[idy].genes[j].urs ):
                print "Genomes", idx, "and", idy, "differ in URS on gene", j
    