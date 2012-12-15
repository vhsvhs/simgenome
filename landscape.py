from cli import *
from rules import *

class Landscape:
    """The Landscape class holds an array of Rule_Collection objects.  Each Rule_Collection can express
    different overlapping expression patterns.  For example, one Rule_Collection can define the expression
    pattern needed for cell lifecycle, while a second Rule_Collection can define the expression
    patterns for a stress response."""
    rulecollections = []
    r = None
    t_counter = 0
    
    def __init__(self, ap):
        self.rulecollections = {} # array of Rule_Collection objects
        self.r = None
        self.t_counter = ap.params["generation"]
    
    def init(self, ap, genome=None, tp=None):
        tp = get_input_rules_from_file(ap)
        
        """Verify that we found input patterns and fitness rules."""
        if tp == None:
            if comm.Get_rank() == 0:
                print "\n. Hmmm, I didn't find any fitness rules in the file", ap.getOptionalArg("--patternpath")
            exit()
        
        if comm.Get_rank() == 0:
            print "\n. Building the fitness landscape described in", ap.getOptionalArg("--patternpath")
        self.rulecollections = tp
        
        if comm.Get_rank() == 0:
            for t in self.rulecollections:
                print t
        
    def uncollapse(self, data, ap):
        for d in data[0]:
            this_rulecollection = Rule_Collection(d)
            this_rulecollection.uncollapse( data[d] )
            self.rulecollections.append( this_rulecollection )

    
    def collapse(self):
        data_rules = {}
        for t in self.rulecollections:
            data_rules[ t.rid ] = t.collapse()
        return data_rules

    def build_random_fitness_rules(self, genome):
        """Creates a random fitness goal, using genome as the seed."""    
        """First, select a basal TR to activate the goal."""
        randid = random.randint(0, ap.params["numtr"]-1)
        
        """Next, build a random pattern.
        For now, the random Landscape contains only one time pattern."""
        rand_rulecollection = Rule_Collection(randid)
        rand_rulecollection.init_basic_test()
        self.rulecollections.append( rand_rulecollection )
    

    
    def fitness_helper(self, gene_expr_level, rid):
        """Returns the fitness of an individual, given their expression levels for all genes."""
        expr_error = 0.0
        max_expr_error = 0.0
        for rule in self.rulecollections[rid].rules:
            obs_expr = gene_expr_level[rule.reporter_id][rule.timepoint]
            if False == rule.rule_type(obs_expr, rule.expression_level):
                expr_error += abs(obs_expr - rule.expression_level)
            max_expr_error += 1.0
        expr_error = expr_error / max_expr_error # Normalize the expression error by the max. possible error
        return math.exp(FITNESS_SCALAR * expr_error) # FITNESS_SCALAR is defined in configuration.py
                    
    
    def get_fitness(self, genome, ap, ko=[]):
        """Calculates the fitness of the given genome, over all time patterns in the fitness landscape.
        Any gene IDs specified in the list 'ko' will be set to MINIMUM_ACTIVITY_LEVEL at each timeslice."""
                
        if ap.params["verbosity"] >= 99:
            notime = 0.0
            timestart = datetime.utcnow()
                
        # First, build the r vector...
        self.r = []
        self.maxr = 0
        for x in ap.params["rangetrs"]:
            self.r.append( genome.genes[x].pwm.P.__len__() )
            if self.r[x] > self.maxr:
                self.maxr = self.r[x]
        # Deal with the no-occupancy case
        self.r.append(1)
                
        rid_fitness = {}
        
        # for each rule pattern:
        for rid in self.rulecollections:
            """Next, initialize expression levels, using either epigenetically inherited levels,
            or, as default, set expression to zero.
            If epigenetics is enabled, then gene expression has already been established at this point."""
            if ap.params["enable_epigenetics"] == False:
                for gene in genome.genes:
                    genome.gene_expr[gene.id] = [MINIMUM_ACTIVITY_LEVEL]
            
            """The expression levels for each TF will be stored in this hashtable."""
            tf_expr_level = {} # key = gene id, value = array of expression values for each timeslice
                    
            for timeslice in range(0, ap.params["maxtime"]):
                self.t_counter = timeslice
    
                """Manually set expression values for genes defined in the input rules."""   
                if timeslice in self.rulecollections[rid].inputs:
                    for i in self.rulecollections[rid].inputs[timeslice]: # for each input rule i...
                        this_gene = i[0]
                        this_expr_level = i[1]
                        genome.gene_expr[this_gene][timeslice] = this_expr_level
                
                """Knock-out any genes that have been specified to be KO'd."""
                """NOTE: if there is a rule for the KO'd gene, the KO wins over the input rule."""
                for gid in ko:
                    genome.gene_expr[gid][timeslice] = MINIMUM_ACTIVITY_LEVEL
    
                """Print a special line, but only if it's the 0th timeslice."""
                if timeslice == 0:
                    if ap.params["verbosity"] > 5:
                        """Print a report to the screen."""
                        for gene in genome.genes:
                            if ap.params["verbosity"] > 5:
                                marka = "      "
                                if gene.has_dbd and gene.is_repressor == False:
                                    marka = "[act.]"
                                elif gene.has_dbd and gene.is_repressor:
                                    marka = "[rep.]"
                                print "r:", rid, "gen.:", ap.params["generation"], "\tt 0", "\tID", genome.id, "\tgene", gene.id, marka, "\tpe: n/a", "\texpr: %.3f"%genome.gene_expr[ gene.id ][timeslice]
                        print ""
    
                
                """ Update all TF expression levels."""         
                for tf in ap.params["rangetrs"]:
                    tf_expr_level[ genome.genes[tf].id ] = genome.gene_expr[ genome.genes[tf].id ][timeslice]
                                           
                for gene in genome.genes:                
                    if gene.id in ko:
                        pe = 0.0
                    else:
                        """ Calculate the delta G of binding on the cis-region for every gene.
                        pe ranges from -0.5 (repression) to 0.5 (strong activation)"""
                        pe = self.get_expr_modifier(genome, gene, tf_expr_level, ap)                
                                        
                    if pe > 0.0:
                        expr_modifier = ap.params["growth_rate"] * pe
                    elif pe < 0.0:
                        expr_modifier = ap.params["decay_rate"] * pe
                    else:
                        expr_modifier = 0.0 # this will result in no change to expression level
                    """New expression level equals the old expression level, plus the modification (which may be positive or negative)"""
                    new_expr_level = genome.gene_expr[gene.id][timeslice] + expr_modifier;
                    genome.gene_expr[gene.id].append( new_expr_level )
                    
                    """ Prevent out-of-bounds expression levels."""
                    if genome.gene_expr[gene.id][timeslice+1] > MAXIMUM_ACTIVITY_LEVEL:
                        genome.gene_expr[gene.id][timeslice+1] = MAXIMUM_ACTIVITY_LEVEL
                    if genome.gene_expr[gene.id][timeslice+1] < MINIMUM_ACTIVITY_LEVEL:
                        genome.gene_expr[gene.id][timeslice+1] = MINIMUM_ACTIVITY_LEVEL
                
                    """ Print a report to stdout."""
                    if ap.params["verbosity"] > 5:                    
                        expr_delta = 0.0
                        expr_delta = genome.gene_expr[ gene.id ][timeslice+1] - genome.gene_expr[ gene.id ][timeslice]
                        marka = "      "
                        if gene.has_dbd and gene.is_repressor == False:
                            marka = "[act.]"
                        elif gene.has_dbd and gene.is_repressor:
                            marka = "[rep.]"
                        print "r:", rid, "gen.:", ap.params["generation"], "\tt", timeslice+1, "\tID", genome.id, "\tgene", gene.id, marka, "\tpe: %.3f"%pe, "\texpr: %.3f"%genome.gene_expr[ gene.id ][timeslice+1], "\td: %.3f"%expr_delta 
    
                #
                # to-do: if gene expression has not changed from the last timeslice
                #    then we've reached equilibrium, so stop cycling through time slices.
                #
                if ap.params["verbosity"] > 5:
                    print ""

            rid_fitness[rid] = self.fitness_helper(genome.gene_expr, rid)
        
        """ Calculate fitness over all rule conditions."""
        fitness = 0.0
        countrid = 0
        for rid in rid_fitness:
            countrid += 1.0
            fitness += rid_fitness[rid]
        fitness /= countrid
        
        if ap.params["verbosity"] > 5:
            print "\n. At generation", ap.params["generation"], ", individual", genome.id, "has fitness %.5f"%fitness, "\n"
        if ap.params["verbosity"] >= 99:
            timeend = datetime.utcnow()
            print ". Generation", ap.params["generation"], "for individual", genome.id, "took %.3f"%(timeend-timestart).total_seconds(), "seconds.\n"
        return fitness
    
    def get_expr_modifier(self, genome, gene, tf_expr_levels, ap):
        """returns a floating-point value, the expression level of gene, given the TF expression levels"""
        pe = []        
        ptables = ProbTable( ap.params["numtr"], ap.params["maxgd"], gene.urs.__len__() )
        #print ptables
        ptables = self.calc_prob_tables(genome, gene, tf_expr_levels, ptables, ap)          
        pe = self.prob_expr(genome, ptables, gene, tf_expr_levels, ap)        
        pe = pe - 0.5 # This will make pe range from -0.5 to +0.5.  
        return pe
    
    def calc_prob_tables(self, genome, gene, rel_tf_expr, ret, ap):
        """returns a ProbTable object, named ret."""
        if ap.params["verbosity"] > 100:
            print "\n\n. CALC_PROB_TABLES gene", gene.id
        L = gene.urs.__len__()        
        for x in range(0, L): # foreach site in gene's upstream region
            sum_cpr = 0.0
            for i in ap.params["rangetrs"]:   # foreach transcription factor
                pwm_tmp = genome.genes[i].pwm.specificity( x, gene.urs )
                #print "TF", i, "binds", gene.urs, "at site", x, "with %.3f"%pwm_tmp, "bits."
                sum_cpt = 0.0
                for j in ap.params["rangetrs+"]: # +1 to also consider the empty case
                    for d in ap.params["rangegd"]:                        
                        if (L-x < self.r[i]): 
                            """ CASE 1: TF i's PWM is too wide to start binding at site x"""
                            ret.cpa[i,j,d,x] = 0 # P @ x = 0
                            sum_cpt += ret.cpa[i,j,d,x]
                            sum_cpr += ret.cpa[i,j,d,x]
                            if ap.params["verbosity"] > 100:
                                print "case 1:", i, j, d, x, ret.cpa[i,j,d,x]
                            continue
                        if (j < ap.params["numtr"]):  # TF j is real, not the empty slot.
                            if (d < MIN_TF_SEPARATION):
                                # If d == 0, then this case should never be reached.
                                """ CASE 4: the distance between TFs i and j is too small."""
                                # we forbid TFs to bind this close together
                                ret.cpa[i,j,d,x] = 0 # then P @ x = 0
                                if ap.params["verbosity"] > 100:
                                    print "case 4:", i, j, d, x, ret.cpa[i,j,d,x]
                                continue
                            elif (L - self.r[i] - d - self.r[j] < 0):
                                """CASE 3: TF i and j cannot both fit on the sequence."""
                                ret.cpa[i,j,d,x] = 0 # then P @ x = 0
                                if ap.params["verbosity"] > 100:
                                    print "case 3:", i, j, d, x, ret.cpa[i,j,d,x]
                                continue
                            else:
                                """Case 5: TF i and TF j can fit with distance d between them"""
                                ret.cpa[i,j,d,x] = rel_tf_expr[i] * pwm_tmp #* (genome.genes[i].gamma[j, d] + genome.genes[j].gamma[i,d]) 
                                sum_cpt += ret.cpa[i,j,d,x]
                                sum_cpr += ret.cpa[i,j,d,x]
                                if ap.params["verbosity"] > 100:
                                    print "case 5:", i, j, d, x, ret.cpa[i,j,d,x]
                                continue
                        elif (j == ap.params["numtr"]): # j is the empty slot
                            """CASE 6: there is no j:"""
                            ret.cpa[i,j,d,x] = rel_tf_expr[i] * pwm_tmp;
                            sum_cpt += ret.cpa[i,j,d,x]
                            sum_cpr += ret.cpa[i,j,d,x]
                            if ap.params["verbosity"] > 100:
                                print "case 6:", i, j, d, x, ret.cpa[i,j,d,x]
                            continue
                ret.cpt[i,x] = sum_cpt
            ret.cpr[x] = sum_cpr
        #if gene.id == 2:
        #    print "\n\nPTABLE FOR GENE 2:"
        #    print ret
        return ret
        
    def sample_cdf(self, site, ptables, ap):
        """Picks a random configuration starting at site, drawing IID from ptables.
        ptables.cpr must contain the cummulative marginal distributions.
        This method returns [i,j,d], where TF is binding to site, followed by distance
        d, and then TF j binds."""
        
        # Use Numpy's cumsum method to calculate a cummulative probability distribution for all binding events at site.
        # c will be a flattened version of ptables.cpa, for the subarray specific to 'site'.
#        c = cumsum(ptables.cpa[:,:,:,site])
#        #print "c=", c
#        # Draw a random value from the CDF. . .
#        randp = random.random() * ptables.cpr[site]
#        # Find where that random value lives in ptables.cpa. . .
#        x = searchsorted(c, randp)
#        #print "randp=", randp
#        #print "x=", x
#        reti = x/(ap.params["numtr"]+1)/MAX_GD
#        retj = reti/MAX_GD
#        retd = x%MAX_GD
#        #print "new answer:", [reti, retj, retd]
#        return [reti, retj, retd]
        
    
        if ptables.cpr[site] == 0.0:
            return None
#            return [ap.params["numtr"], ap.params["numtr"], 1]
        
        
        #
        # This is the old way....
        #  
        randp = random.random() * ptables.cpr[site]
        sump = 0.0
        reti = 0
        retj = 0
        retd = 0
        
        # for debugging: as a test, just pick a random cell from the appropriate cpa
        # range, rather than picking an IID cell from the cpr-based range.
        #i = random.randint(0, ap.params["numtr"]-1)
        #j = random.randint(0, ap.params["numtr"])
        #d = random.randint(0, MAX_GD-1)
        
        
        for i in ap.params["rangetrs"]:
            reti = i
            for j in ap.params["rangetrs+"]:
                retj = j
                for d in ap.params["rangegd"]:
                    retd = d
                    sump += ptables.cpa[i,j,d,site]
                    #this next line is useful for debugging, but it totally explodes the runtime:
                    #print "landscape 340, site", site, "tf", i, "tf", j, "d", d, "sump", sump, "randp", randp
                    if sump > randp:
                        break
                if sump > randp:
                    break
            if sump > randp:
                break
        #print "old answer:", [reti, retj, retd]
        return [reti, retj, retd]
        
    
    def prob_expr(self, genome, ptables, gene, tf_expr_levels, ap, lam=None):        
        """Returns a floating-point value corresponding to the expression level of gene,
        given the ProbTable ptables"""        
        #print "prob_expr, gene", gene.id, gene.urs
        
        pe_sum = 0
        min_r = min( self.r )
        configurations = {}
        """configurations: key = site, value = array of arrays, [i,j,d] samples"""
        
        for sample in range(0, ap.params["iid_samples"]):
            
            """"1. Build a configuration this_config, by sampling cells from ptables.cpa"""          
            this_config = {} # key = site, value = the TF bound starting at this site.
            site = 0
            while (site < gene.urs.__len__()):
                if ptables.cpr[site] == 0.0:
                    """Skip sites where nothing is binding...."""
                    site += 1
                else:
                    """Sample a configuration..."""        
                    [i, j, d] = self.sample_cdf(site, ptables, ap)                                  
                    this_config[site] = i

                    
                    """Save the configuration..."""
                    if False == configurations.__contains__(site):
                        configurations[site] = []
                    configurations[site].append( [i,j,d] )

                    """Advance the site counter..."""
                    if i < ap.params["numtr"]:
                        site += genome.genes[i].pwm.P.__len__()
                    if site < gene.urs.__len__():
                        site += d
                    if j < ap.params["numtr"]:
                        if (site + genome.genes[j].pwm.P.__len__()) < gene.urs.__len__():
                            this_config[site] = j
                            site += genome.genes[j].pwm.P.__len__()
                
            """"2. Calculate the binding energy of the configuration."""
            sum_lambda_act = 0.0
            sum_lambda_rep = 0.0
            for site in this_config:
                tf = this_config[site]
                if tf < ap.params["numtr"]: # there will be no binding energy for the empty configuration:
                    """Get the strength of TF binding at this site..."""
                    tf_specificity = genome.genes[tf].pwm.specificity(site, gene.urs)
                    #print "tf", tf, "binds gene", gene.id, "site", site, "tf_specificity", tf_specificity
                    if genome.genes[tf].is_repressor:
                        sum_lambda_rep += tf_specificity
                    else:
                        sum_lambda_act += tf_specificity
            
            """ 3. the probability of this configuration equals:
            # this is what Kevin did:
            #pe_sum = (1/(1+math.exp(-1*sum_lambda)))
            # but this is what I'm doing instead:
            # This incorporates the hill equation"""
            

            k_act = sum_lambda_act
            k_rep = sum_lambda_rep
            this_pe = (1/( 1+math.exp(-1*ap.params["pe_scalar"]*(k_act-k_rep) ) ))
            pe_sum += this_pe
            #print "debug landscape.py 417 - gene", gene.id, "k_act", k_act, "k_rep", k_rep, "this_pe", this_pe
        
        if ap.params["verbosity"] > 3:
            """Print what we've sampled..."""
            self.print_configuration(configurations, genome, gene, ap)
        return (pe_sum / ap.params["iid_samples"])
    
    
    def print_configuration(self, configs, genome, gene, ap):        
        
        #if gene.id == 2:
        #    for key in configs:
        #        for c in configs[key]:
        #            print key, c, configs[key].__len__()
        
        """configs is a hashtable of configurations...
            configs[site] = array of triples [gene i, gene j, distance]"""
        foutpath = ap.params["workspace"] + "/" + ap.params["runid"] + "/" + EXPR_PLOTS + "/config.gen" + ap.params["generation"].__str__() + ".gid" + genome.id.__str__() + ".txt"
        try:
            fout = open( foutpath , "a")
        except IOError:
            print "\n. I had a problem opening", foutpath
            exit()

        if self.t_counter == 0 and gene.id == 0:        
            fout.write("Notes\n")
            fout.write("[+] : activator\n")
            fout.write("[-] : repressor\n")
            fout.write("p   : proportion of IID samples\n")
            fout.write("e   : energy bound by TF in this configuration\n\n")
        
        #fout.write("#\n# For example, \"site 3 :  1 [-] 0.672 0.120\" indicates that site 3 is bound by repressor 1\n# in 67.2 percent of IID samples with activation energy = 0.120.\n# Activation energies < 0.5 indicate repression, = 0.5 indicate no activation,\n# > 0.5 indicates activation.#\n#\n")
        fout.write(". TIME " + self.t_counter.__str__() + " GENE " + gene.id.__str__() + "\t" + gene.urs + "\n")
            
        """configs is array of arrays, [site, tf_i, tf_j, distance between i and j]"""
        sites = configs.keys()
        sites.sort()
        tf_count = {} # key = site, value = hash, with key = tf, value = proportion of samples
        for site in sites:
            if site not in tf_count:
                tf_count[site] = {}
            
    
        # transform configs into tf_count
        for site in sites:
            for c in configs[site]:
                i = c[0]
                j = c[1]
                d = c[2]
                if i not in tf_count[site]:
                    tf_count[site][i] = 0.0
                tf_count[site][i] += 1.0 / ap.params["iid_samples"]
                
                if j < ap.params["numtr"]:
                    adj_site = site + self.r[i] + d
                    if adj_site + self.r[j] < gene.urs.__len__():                
                        if adj_site not in sites:
                            tf_count[adj_site] = {}
                        if False == tf_count[adj_site].__contains__(j):
                            tf_count[adj_site][j] = 0.0
                        tf_count[adj_site][j] += 1.0 / ap.params["iid_samples"]

        # print tf_count
        insites = tf_count.keys()
        insites.sort()
        for site in insites:
            line = None
            for tf in tf_count[site]:
                if tf < ap.params["numtr"]:
                    if line == None:
                        line = "site " + site.__str__() + " :"
                    pe = genome.genes[tf].pwm.specificity(site, gene.urs)
                    if genome.genes[tf].is_repressor:
                        tf_type = "[-]"
                    else:
                        tf_type = "[+]"
                    line += " \t" + tf.__str__() + " " + tf_type + " p=%.3f"%tf_count[site][tf] + " e=%.3f"%pe
            if line != None:
                fout.write(line + "\n")
        fout.write("\n")
        fout.close()
        
        
        