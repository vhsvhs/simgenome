from cli import *
from rules import *

class Landscape:
    """The Landscape class holds an array of Rule_Collection objects.  Each Rule_Collection can express
    different overlapping expression patterns.  For example, one Rule_Collection can define the expression
    pattern needed for cell lifecycle, while a second Rule_Collection can define the expression
    patterns for a stress response."""
    rulecollection = []
    inputpatterns = {} # key = timepoint, value = array of tuples (gene, expr. value)
    gamma = None
    r = None
    t_counter = 0
    
    def __init__(self, ap):
        self.rulecollection = [] # array of Rule_Collection objects
        self.inputpatterns = {}  # hashtable, key = timepoint, value = array of tuples (gene, expression level).
        self.gamma = None
        self.r = None
        self.t_counter = ap.params["generation"]
    
    def init(self, ap, genome=None, tp=None):
        [tp,ip] = get_input_rules_from_file(ap)
        
        """Verify that we found input patterns and fitness rules."""
        if tp == None:
            if comm.Get_rank() == 0:
                print "\n. Hmmm, I didn't find any fitness rules in the file", ap.getOptionalArg("--patternpath")
            exit()
        if ip == None:
            if comm.Get_rank() == 0:
                print "\n. Hmmm, I didn't find any input patterns in the file", ap.getOptionalArg("--patternpath")
            exit()
        
        if comm.Get_rank() == 0:
            print "\n. Building the fitness landscape described in", ap.getOptionalArg("--patternpath")
        self.rulecollection = tp
        self.inputpatterns = ip
        self.set_gamma(ap)
        
        if comm.Get_rank() == 0:
            for t in self.rulecollection:
                print t
        
    def uncollapse(self, data, ap):
        for d in data[0]:
            this_rulecollection = Rule_Collection(d)
            this_rulecollection.uncollapse( data[d] )
            self.rulecollection.append( this_rulecollection )
            self.set_gamma(ap)
        for t in data[1]:
            self.inputpatterns[t] = data[1][t]
    
    def collapse(self):
        data_rules = {}
        for t in self.rulecollection:
            data_rules[ t.basal_gene_id ] = t.collapse()
        
        data_inputs = {}
        for t in self.inputpatterns:
            data_inputs[ t ] = self.inputpatterns[t]
        return [data_rules, data_inputs]

    def build_random_fitness_rules(self, genome):
        """Creates a random fitness goal, using genome as the seed."""    
        """First, select a basal TR to activate the goal."""
        randid = random.randint(0, ap.params["numtr"]-1)
        
        """Next, build a random pattern.
        For now, the random Landscape contains only one time pattern."""
        rand_rulecollection = Rule_Collection(randid)
        rand_rulecollection.init_basic_test()
        self.rulecollection.append( rand_rulecollection )
    
    def coopfunc(self, g, d):
        """Calculates the degree of binding cooperativity between two TFs binding distance d apart.
        g is positive for synergistic interactions, and negative for antagonistic interactions.
        g ranges from -1 to +inifinity."""
        return 1 + g * math.exp( (-1)*(d**2)/V_RATE_OF_COOP_DECAY );

    def set_gamma(self, ap):
        """Prec-calculates the cooperative/competitive binding interactions between all TFs."""
        # g is the binding coop term (see the coopfunc).
        # for now, g is all zeroes, so there is no cooperative binding.
        #
        # to-do: grab gamma values from the command-line, or radomly initialize them
        # from an a priori distribution.
        g = zeros( (ap.params["numtr"], ap.params["numtr"]), dtype=float)
        self.gamma = zeros( (ap.params["numtr"], ap.params["numtr"], ap.params["maxgd"]), dtype=float)
        for i in ap.params["rangetrs"]:
            for j in ap.params["rangetrs"]:
                for d in ap.params["rangegd"]:
                    print i, j, d
                    self.gamma[i,j,d] = self.coopfunc( g[i,j], d)
    
    def fitness_helper(self, gene_expr_level):
        expr_error = 0.0
        max_expr_error = 0.0
        for pattern in self.rulecollection:
            for rule in pattern.rules:
                obs_expr = gene_expr_level[rule.reporter_id][rule.timepoint]
                if False == rule.rule_type(obs_expr, rule.expression_level):
                    expr_error += abs(obs_expr - rule.expression_level)
                max_expr_error += 1.0
        return math.exp(FITNESS_SCALAR * expr_error)
                    
    
    def get_fitness(self, genome, ap):
        """Calculates the fitness of the given genome, over all time patterns in the fitness landscape.
        Returns a floating-point value."""
                
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
        
        # Next, initialize expression levels, using either epigenetically inherited levels,
        # or, as default, set expression to zero.
        # If epigenetics is enabled, then gene expression has already been established at this point.
        if ap.params["enable_epigenetics"] == False or genome.gene_expr.__len__() < 1:
            for gene in genome.genes:
                genome.gene_expr[gene.id] = [MINIMUM_ACTIVITY_LEVEL]
        
        # The expression levels for each TF will be stored in this hashtable.
        tf_expr_level = {} # key = gene id, value = array of expression values for each timeslice
                
        for timeslice in range(0, ap.params["maxtime"]):
            self.t_counter = timeslice

            # manually set expression values for genes defined in the input rules.   
            if timeslice in self.inputpatterns:
                for i in self.inputpatterns[timeslice]:
                    this_gene = i[0]
                    this_expr_level = i[1]
                    genome.gene_expr[this_gene][timeslice] = this_expr_level

            # print a special line, but only if it's the 0th timeslice.
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
                            print "gen.", ap.params["generation"], "\tt 0", "\tID", genome.id, "\tgene", gene.id, marka, "\tact: n/a", "\texpr: %.3f"%genome.gene_expr[ gene.id ][timeslice]
                    print ""

            
            # Update TF expression levels...           
            for tf in ap.params["rangetrs"]:
                tf_expr_level[ genome.genes[tf].id ] = genome.gene_expr[ genome.genes[tf].id ][timeslice]
                                       
            for gene in genome.genes:                
                # Calculate the delta G of binding on the cis-region for every gene.
                pe = self.get_expression(genome, gene, tf_expr_level, ap)                
                expr_modifier = 2.0 * pe
                new_expr_level = genome.gene_expr[gene.id][timeslice] * expr_modifier;
                genome.gene_expr[gene.id].append( new_expr_level )
                
                # Prevent out-of-bounds expression levels.
                if genome.gene_expr[gene.id][timeslice+1] > MAXIMUM_ACTIVITY_LEVEL:
                    genome.gene_expr[gene.id][timeslice+1] = MAXIMUM_ACTIVITY_LEVEL
                if genome.gene_expr[gene.id][timeslice+1] < MINIMUM_ACTIVITY_LEVEL:
                    genome.gene_expr[gene.id][timeslice+1] = MINIMUM_ACTIVITY_LEVEL
            
                # Print a report to stdout.
                if ap.params["verbosity"] > 5:                    
                    expr_delta = 0.0
                    expr_delta = genome.gene_expr[ gene.id ][timeslice+1] - genome.gene_expr[ gene.id ][timeslice]
                    marka = "      "
                    if gene.has_dbd and gene.is_repressor == False:
                        marka = "[act.]"
                    elif gene.has_dbd and gene.is_repressor:
                        marka = "[rep.]"
                    print "gen.", ap.params["generation"], "\tt", timeslice+1, "\tID", genome.id, "\tgene", gene.id, marka, "\tact: %.3f"%pe, "\texpr: %.3f"%genome.gene_expr[ gene.id ][timeslice+1], "\td: %.3f"%expr_delta 

            #
            # to-do: if gene expression has not changed from the last timeslice
            #    then we've reached equilibrium, so stop cycling through time slices.
            #
            if ap.params["verbosity"] > 5:
                print ""

        fitness = self.fitness_helper(genome.gene_expr)
        if ap.params["verbosity"] > 5:
            print "\n. At generation", ap.params["generation"], ", individual", genome.id, "has fitness %.5f"%fitness, "\n"
        if ap.params["verbosity"] >= 99:
            timeend = datetime.utcnow()
            print ". Generation", ap.params["generation"], "for individual", genome.id, "took %.3f"%(timeend-timestart).total_seconds(), "seconds.\n"
        
        return fitness
    
    def get_expression(self, genome, gene, tf_expr_levels, ap):
        """returns a floating-point value, the expression level of gene, given the TF expression levels"""
        pe = []        
        ptables = ProbTable( ap.params["numtr"], ap.params["maxgd"], gene.urs.__len__() )
        ptables = self.calc_prob_tables(genome, gene, tf_expr_levels, ptables, ap)  
        #print ptables
        #exit()
        return self.prob_expr(genome, ptables, gene, tf_expr_levels, ap)          
    
    def calc_prob_tables(self, genome, gene, rel_tf_expr, ret, ap):
        """returns a ProbTable object, named ret."""
        if ap.params["verbosity"] > 100:
            print "\n\n. CALC_PROB_TABLES gene", gene.id
        L = gene.urs.__len__()        
        for x in range(0, L): # foreach site in gene's upstream region
            sum_cpr = 0.0
            for i in ap.params["rangetrs"]:   # foreach transcription factor
                pwm_tmp = genome.genes[i].pwm.prob_binding( x, gene.urs )
                #print "TF", i, "binds", gene.urs, "at site", x, "with %.3f"%pwm_tmp, "bits."
                sum_cpt = 0.0
                for j in ap.params["rangetrs+"]: # +1 to also consider the empty case
                    #print i, j, x
                    for d in ap.params["rangegd"]:                        
                        if (L-x < self.r[i]): 
                            """ CASE 1: TF i's PWM is too wide to start binding at site x"""
                            ret.cpa[i,j,d,x] = 0 # P @ x = 0
                            sum_cpt += ret.cpa[i,j,d,x]
                            sum_cpr += ret.cpa[i,j,d,x]
                            if ap.params["verbosity"] > 100:
                                print "case 1:", i, j, d, x, ret.cpa[i,j,d,x]
                            continue
                        # The implied else case is that TF i can fit at site x.
                        # Victor: I think we should first do cases 3, 4, and 5
#                        elif (L-x >= self.r[i] and (j == i or j == ap.params["numtr"]) and d == 0): # // if TF i can bind here
#                            """CASE 2: TF i can start binding at site x
#                            and TF i is identical to j."""
#                            ret.cpa[i,j,d,x] = rel_tf_expr[i] * pwm_tmp # basic case, no competition or cooperation
#                            sum_cpt += ret.cpa[i,j,d,x]
#                            sum_cpr += ret.cpa[i,j,d,x]
#                            if ap.params["verbosity"] > 100:
#                                print "case 2:", i, j, d, x, ret.cpa[i,j,d,x]
#                            continue
                        if (j < ap.params["numtr"]):  # TF j is real, not the empty slot.
                            if (d < MIN_TF_SEPARATION):
                                # If d == 0, then this case should never be reached.
                                """ CASE 4: the distance between TFs i and j is too small."""
                                # we forbid TFs to bind this close together
                                ret.cpa[i,j,d,x] = 0 # then P @ x = 0
                                sum_cpt += ret.cpa[i,j,d,x]
                                sum_cpr += ret.cpa[i,j,d,x]
                                if ap.params["verbosity"] > 100:
                                    print "case 4:", i, j, d, x, ret.cpa[i,j,d,x]
                                continue
                            elif (L - self.r[i] - d - self.r[j] < 0):
                                """CASE 3: TF i and j cannot both fit on the sequence."""
                                ret.cpa[i,j,d,x] = 0 # then P @ x = 0
                                sum_cpt += ret.cpa[i,j,d,x]
                                sum_cpr += ret.cpa[i,j,d,x]
                                if ap.params["verbosity"] > 100:
                                    print "case 3:", i, j, d, x, ret.cpa[i,j,d,x]
                                continue
                            else:
                                ret.cpa[i,j,d,x] = rel_tf_expr[i] * pwm_tmp * self.gamma[j, i, d] 
                                sum_cpt += ret.cpa[i,j,d,x]
                                sum_cpr += ret.cpa[i,j,d,x]
                                if ap.params["verbosity"] > 100:
                                    print "case 5:", i, j, d, x, ret.cpa[i,j,d,x]
                                continue
#                            elif (i != j):
#                                """CASE 5: TFs i and j are different AND they can both fit on the URS..."""
#                                ret.cpa[i,j,d,x] = rel_tf_expr[i] * pwm_tmp * self.gamma[j, i, d]
#                                sum_cpt += ret.cpa[i,j,d,x]
#                                sum_cpr += ret.cpa[i,j,d,x]
#                                if ap.params["verbosity"] > 100:
#                                    print "case 5:", i, j, d, x, ret.cpa[i,j,d,x]
#                                continue
#                            elif (i == j):
#                                ret.cpa[i,j,d,x] = rel_tf_expr[i] * pwm_tmp # basic case, no competition or cooperation
#                                sum_cpt += ret.cpa[i,j,d,x]
#                                sum_cpr += ret.cpa[i,j,d,x]
#                                if ap.params["verbosity"] > 100:
#                                    print "case 2:", i, j, d, x, ret.cpa[i,j,d,x]
#                                continue
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
                    #print "site", site, "tf", i, "tf", j, "d", d, "sump", sump, "randp", randp
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
            
            """"1. build a configuration c_k, by sampling cells from ptables.cpa"""          
            this_config = {} # key = site, value = the TF bound starting at this site.
            site = 0
            while (site < gene.urs.__len__()):
                [i, j, d] = self.sample_cdf(site, ptables, ap)
                
                #
                # Useful output for debugging:
                #
                #print "\n. Sampled i=", i, "j=", j, "d=", d
                #if i < ap.params["numtr"]:
                #    print "i:",i, genome.genes[i].pwm.P
                #if j < ap.params["numtr"]:
                #    print "j:",j,genome.genes[j].pwm.P
                #print "URS:", gene.urs
                #print i, j, d
                #
                #
                
                if False == configurations.__contains__(site):
                    configurations[site] = []
                configurations[site].append( [i,j,d] )
                this_config[site] = i
                if i < ap.params["numtr"]:
                    site += genome.genes[i].pwm.P.__len__()
                else:
                    site += 1 # for the non-occuped case
                site += d
                
                if site < gene.urs.__len__(): # don't create a configuration beyond the length of the URS
                    if False == configurations.__contains__(site):
                        configurations[site] = []
                    this_config[site] = j
                    if j < ap.params["numtr"]:
                        site += genome.genes[j].pwm.P.__len__()
                    else:
                        site += 1
                
            """"2. Calculate the binding energy of the configuration."""
            sum_lambda_act = 0.0
            sum_lambda_rep = 0.0
            for site in this_config:
                tf = this_config[site]
                if tf != ap.params["numtr"]: # there will be no binding energy for the empty configuration:
                    """Get the strength of TF binding at this site..."""
                    tf_specificity = genome.genes[tf].pwm.prob_binding(site, gene.urs)
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
            #print "gene", gene.id, "k_act", k_act, "k_rep", k_rep, "this_pe", this_pe
        if ap.params["verbosity"] > 3:
            self.print_configuration(configurations, genome, gene, ap)
        return (pe_sum / ap.params["iid_samples"])
    
    
    def print_configuration(self, configs, genome, gene, ap):        
        """configs is a hashtable of configurations...
            configs[site] = array of triples [gene i, gene j, distance]"""
        foutpath = ap.getArg("--workspace") + "/" + ap.getArg("--runid") + "/" + EXPR_PLOTS + "/config.gen" + ap.params["generation"].__str__() + ".gid" + genome.id.__str__() + ".txt"
        fout = open( foutpath , "a")
        fout.write(". TIME " + self.t_counter.__str__() + " GENE " + gene.id.__str__() + "\t" + gene.urs + "\n")
            
        """configs is array of arrays, [site, tf_i, tf_j, distance between i and j]"""
        sites = configs.keys()
        sites.sort()
        for site in sites:
            tf_count = {}
            for c in configs[site]:
                if False == tf_count.__contains__( c[0] ):
                    tf_count[ c[0] ] = 0
                tf_count[ c[0] ] += 1.0 / ap.params["iid_samples"]

            line = "site " + site.__str__() + " :"
            for tf in tf_count:
                pe = 0.0
                if tf < ap.params["numtr"]:
                    pe = genome.genes[tf].pwm.prob_binding(site, gene.urs)
                tf_type = "[+]"
                if genome.genes[tf].is_repressor:
                    tf_type = "[-]"
                line += " \t" + tf.__str__() + " " + tf_type + " %.3f"%tf_count[tf] + " %.3f"%pe
            fout.write(line + "\n")
        fout.write("\n")
        fout.close()

    
      