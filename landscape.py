from configuration import *
from genome import *
from probtable import *

"""RULE_TYPES defines a set of boolean comparison criteria for comparing
the observed expression levels of reporter genes to their optimal levels.
These operator functions will be applied as follows:
operator.f( observed expression, optimal expression)
where 'f' is the operator function."""
RULE_TYPES = [operator.ge, operator.eq, operator.le]

class Fitness_Rule:
    timepoint = None
    expression_level = None
    reporter_id = None
    rule_type = None # points to a boolean operator function.

    def __init__(self, t, e, id, r):
        self.timepoint = t
        self.expression_level = e
        self.reporter_id = id
        self.rule_type = r
    
    def __str__(self):
        l = "\t+ At time "
        l += self.timepoint.__str__()
        l += " expression of gene "
        l += self.reporter_id.__str__()
        l += " must be "
        if self.rule_type == operator.ge:
            l += ">"
        if self.rule_type == operator.le:
            l += "<"
        if self.rule_type == operator.eq:
            l += "="
        l += " "
        l += self.expression_level.__str__()
        l += ".\n"
        return l
    
    def collapse(self):
        return [self.timepoint, self.expression_level, self.reporter_id, self.rule_type]

class Time_Pattern:
    """Each instance of the class Time_Pattern defines a unique set of temporal expression patterns
    for a subset of reporter genes in each genome.  See the parent class, Landscape, for how to
    calculate the fitness of a given Time_Pattern.  This class can be initialized with an a priori
    desired expression pattern, or the pattern can be randomly initialized.
    
    The optimal pattern is defined using the following information, for example:
        reporter x, above or equal to 50% expression @ timepoint A, max expression at timepoint B, below 50% expression @ timepoint C.
        reporter y, 0% expression @ timepoint A, 0% expression at timepoint C, above or equal to 50% at timepoint D."""

    basal_gene_id = None
    rules = []
    
    def __init__(self, basal_gene_id):
        self.basal_gene_id = basal_gene_id
        self.rules = []
        
    def collapse(self):
        data = []
        for r in self.rules:
            #print "debug 72, rule"
            data.append( r.collapse() )
        return data
    
    def uncollapse(self, data):
        for rule in data:
            #print "debug 78, rule"
            this_rule = Fitness_Rule(rule[0], rule[1], rule[2], rule[3])
            self.rules.append( this_rule )
    
    def init_random(self):
        last_time = 0
        for i in range(0, ap.params["ngoals"]):
            """First, expression > 0.5"""
            last_time += random.randint(1, GOAL_PULSE_MAX)
            if last_time > ap.params["maxtime"]:
                last_time = ap.params["maxtime"]
            e = 0.5
            rand_ri = random.randint(ap.params["numtr"], N_REPORTER+ap.params["numtr"]-1)
            r = RULE_TYPES[0] 
            self.rules.append( Fitness_Rule(last_time, e, rand_ri, r) )
            """Next, expression at maximum"""     
            last_time += random.randint(1, GOAL_PULSE_MAX)
            e = 1.0
            r = RULE_TYPES[1]
            if last_time > ap.params["maxtime"]:
                last_time = ap.params["maxtime"]
            self.rules.append( Fitness_Rule(last_time, e, rand_ri, r) )
            """Finally, expression < 0.5"""
            last_time += random.randint(1, GOAL_PULSE_MAX)
            e = 0.5
            r = RULE_TYPES[2]
            if last_time > ap.params["maxtime"]:
                last_time = ap.params["maxtime"]
            self.rules.append( Fitness_Rule(last_time, e, rand_ri, r) )
    
    def init_basic_test(self):
        """This method is for debugging/testing only."""
        time = ap.params["maxtime"]-1
        e = 0.5
        rand_ri = random.randint(ap.params["numtr"], N_REPORTER+ap.params["numtr"]-1) 
        r = RULE_TYPES[0] 
        self.rules.append( Fitness_Rule(time, e, rand_ri, r) )
    
    def __str__(self):
        ret = ""
        ret += "\n+ Rule: activation of basal gene " + self.basal_gene_id.__str__() + " should cause the following:\n"
        for r in self.rules:
            ret += r.__str__()
        return ret
    
class Landscape:
    """An array of Time_Patterns"""
    timepatterns = []
    gamma = None
    r = None
    t_counter = 0
    gen_counter = 0
    
    def __init__(self):
        self.timepatterns = []
        self.gamma = None
        self.r = None
    
    def init(self, ap, genome=None, tp=None):
        if tp == None:
            self.init_random(genome)
        else:
            self.timepatterns = tp
        self.set_gamma(ap)
        
        for t in self.timepatterns:
            print t
        
    def uncollapse(self, data, ap):
        #print "debug 142"
        for d in data:
            #print "debug 143, time pattern for basal gene id", d
            this_timepattern = Time_Pattern(d)
            this_timepattern.uncollapse( data[d] )
            self.timepatterns.append( this_timepattern )
            self.set_gamma(ap)
    
    def collapse(self):
        data = {}
        for t in self.timepatterns:
            #print "debug 153, timepattern for basal gene id", t.basal_gene_id
            data[ t.basal_gene_id ] = t.collapse()
        return data

    def init_random(self, genome):
        """Creates a random fitness goal, using genome as the seed."""    
        """First, select a basal TR to activate the goal."""
        randid = random.randint(0, ap.params["numtr"]-1)
        #print "basal id = ", randid
        
        """Next, build a random pattern.
        For now, the random Landscape contains only one time pattern."""
        rand_timepattern = Time_Pattern(randid)
        rand_timepattern.init_basic_test()
        self.timepatterns.append( rand_timepattern )
    
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
        self.gamma = zeros( (ap.params["numtr"], ap.params["numtr"], MAX_GD), dtype=float)
        for i in range (0, ap.params["numtr"]):
            for j in range(0, ap.params["numtr"]):
                for d in range(0, MAX_GD):
                    self.gamma[i,j,d] = self.coopfunc( g[i,j], d)
    
    def fitness_helper(self, gene_expr_level):
        expr_error = 0.0
        max_expr_error = 0.0
        for pattern in self.timepatterns:
            for rule in pattern.rules:
                obs_expr = gene_expr_level[rule.reporter_id][rule.timepoint]
                if False == rule.rule_type(obs_expr, rule.expression_level):
                    expr_error += abs(obs_expr - rule.expression_level)
                max_expr_error += 1.0
        #print expr_error
        #print max_expr_error
        return math.exp(FITNESS_SCALAR * expr_error)
        #return (max_expr_error - expr_error) / max_expr_error
                    
    
    def get_fitness(self, genome, generation, ap):
        """Calculates the fitness of the given genome, over all time patterns in the fitness landscape.
        Returns a floating-point value."""
        
        """First, build the r vector..."""
        self.r = []
        self.maxr = 0
        for x in range(0, ap.params["numtr"]):
            self.r.append( genome.genes[x].pwm.P.__len__() )
            if self.r[x] > self.maxr:
                self.maxr = self.r[x]
        """Deal with the no-occupancy case"""
        self.r.append(1)
        
        """ gene_expr[gene ID][timeslice] = expression level of gene at timeslice"""
        genome.gene_expr = {}
        for gene in genome.genes:
            genome.gene_expr[gene.id] = [MINIMUM_ACTIVITY_LEVEL] # all genes begin with zero
            
        """ tf_expr_level[gene ID] = current expression level"""
        """ tf_expr_level is used as short-term variable inside the loop "for timeslice..." """
        list_of_basal_gids = []
        tf_expr_level = {}
        for timepattern in self.timepatterns:
            genome.gene_expr[timepattern.basal_gene_id] = [MAXIMUM_ACTIVITY_LEVEL]
            list_of_basal_gids.append( timepattern.basal_gene_id )
                
        for timeslice in range(0, ap.params["maxtime"]):
            self.t_counter = timeslice
            # 1a. Ensure that basal genes remain activated
            for timepattern in self.timepatterns:
                if genome.gene_expr[timepattern.basal_gene_id][timeslice] < MAXIMUM_ACTIVITY_LEVEL:
                    genome.gene_expr[timepattern.basal_gene_id][timeslice] = MAXIMUM_ACTIVITY_LEVEL
            
            # for debugging:
            #if timeslice > 10:
            #    for timepattern in self.timepatterns:
            #        gene_expr[timepattern.basal_gene_id].append(MINIMUM_ACTIVITY_LEVEL)               
            
            
            """ 1b. Update TF expression levels..."""            
            for tf in range(0, ap.params["numtr"]):
                tf_expr_level[ genome.genes[tf].id ] = genome.gene_expr[ genome.genes[tf].id ][timeslice]
                #print "landscape 270", tf, genome.genes[tf].id
                                       
            for gene in genome.genes:
                if False == list_of_basal_gids.__contains__(gene.id):
                    """Calculate the delta G of binding on the cis-region for every gene."""
                    pe = self.get_expression(genome, gene, tf_expr_level, ap)
                    
                    #expr_modifier = 20.0**pe - 1.0 
                    expr_modifier = 2.0 * pe
                    #print "pe=", pe, "growth_mod=", expr_modifier
                    #if expr_modifier > ap.params["growth_rate"]:
                    #    expr_modifier = ap.params["growth_rate"]
                    new_expr_level = genome.gene_expr[gene.id][timeslice] * expr_modifier;
                    genome.gene_expr[gene.id].append( new_expr_level )
                else:
                    pe = 1.0
                    genome.gene_expr[gene.id].append( genome.gene_expr[gene.id][timeslice] )
                
                # if the delta G is high enough, then transcribe the gene.
                #if pe >= ACTIVATION_THRESHOLD:

                #else: # otherwise, decay the expression level of the gene.
                #    new_expr_level = genome.gene_expr[gene.id][timeslice] / ap.params["decay_rate"];
                #    genome.gene_expr[gene.id].append( new_expr_level )
            
                """ Prevent out of bounds expression. """
                if genome.gene_expr[gene.id][timeslice+1] > MAXIMUM_ACTIVITY_LEVEL:
                    genome.gene_expr[gene.id][timeslice+1] = MAXIMUM_ACTIVITY_LEVEL
                if genome.gene_expr[gene.id][timeslice+1] < MINIMUM_ACTIVITY_LEVEL:
                    genome.gene_expr[gene.id][timeslice+1] = MINIMUM_ACTIVITY_LEVEL
            
                # print report to screen
                if int(ap.getOptionalArg("--verbose")) > 2:
                    expr_delta = 0.0
                    #if timeslice > 1:
                    expr_delta = genome.gene_expr[ gene.id ][timeslice+1] - genome.gene_expr[ gene.id ][timeslice]
                    marka = "*"
                    if gene.has_dbd == False:
                        marka = ""
                    print "gen.", generation, "\tt", timeslice+1, "\tID", genome.id, "\tgene", gene.id, marka, "\tact: %.3f"%pe, "\texpr: %.3f"%genome.gene_expr[ gene.id ][timeslice+1], "\td: %.3f"%expr_delta 

            #
            # to-do: if gene expression has not changed from the last timeslice
            #    then we've reached equilibrium, so stop cycling through time slices.
            #
            print ""

        #
        # to-do: assess fitness of genome, using gene_expr
        #
        fitness = self.fitness_helper(genome.gene_expr)
        print "\n. At gen.", generation, ", individual", genome.id, "has fitness %.5f"%fitness, "\n"
        return fitness
    
    def get_expression(self, genome, gene, tf_expr_levels, ap):
        """returns a floating-point value, the expression level of gene, given the TF expression levels"""
        pe = []        
        ptables = ProbTable( ap.params["numtr"], MAX_GD, gene.urs.__len__() )
        ptables = self.calc_prob_tables(genome, gene, tf_expr_levels, ptables, ap)   
        #print "PTABLE for gene", gene.id
        #print ptables   
        return self.prob_expr(genome, ptables, gene, tf_expr_levels, ap)          
    
    def calc_prob_tables(self, genome, gene, rel_tf_expr, ret, ap):
        """returns a ProbTable object. ret is the ProbTable that should be returned."""
        if int(ap.getOptionalArg("--verbose")) > 100:
            print "\n\n. CALC_PROB_TABLES gene", gene.id
        L = gene.urs.__len__()        
        for x in range(0, L): # foreach site in gene's upstream region
            sum_cpr = 0.0
            for i in range(0, ap.params["numtr"]):   # foreach transcription factor
                pwm_tmp = genome.genes[i].pwm.prob_binding( x, gene.urs )
                #print "TF", i, "binds", gene.urs, "at site", x, "with %.3f"%pwm_tmp, "bits."
                sum_cpt = 0.0
                for j in range(0, ap.params["numtr"]+1): # +1 to also consider the empty case
                    #print i, j, x
                    for d in range(0, MAX_GD):                        
                        if (L-x < self.r[i]): 
                            """ CASE 1: TF i's PWM is too wide to start binding at site x"""
                            ret.cpa[i,j,d,x] = 0 # P @ x = 0
                            sum_cpt += ret.cpa[i,j,d,x]
                            sum_cpr += ret.cpa[i,j,d,x]
                            if int(ap.getOptionalArg("--verbose")) > 100:
                                print "case 1:", i, j, d, x, ret.cpa[i,j,d,x]
                            continue
                        elif (L-x >= self.r[i] and (j == i or j == ap.params["numtr"]) and d == 0): # // if TF i can bind here
                            """CASE 2: TF i can start binding at site x
                            and TF i is identical to j."""
                            ret.cpa[i,j,d,x] = rel_tf_expr[i] * pwm_tmp # basic case, no competition or cooperation
                            sum_cpt += ret.cpa[i,j,d,x]
                            sum_cpr += ret.cpa[i,j,d,x]
                            if int(ap.getOptionalArg("--verbose")) > 100:
                                print "case 2:", i, j, d, x, ret.cpa[i,j,d,x]
                            continue
                        elif (j < ap.params["numtr"]): # else, cooperative/competitive binding...
                            if (d < MIN_TF_SEPARATION):
                                """ CASE 4: the distance between TFs i and j is too small."""
                                # we forbid TFs to bind this close together
                                ret.cpa[i,j,d,x] = 0
                                sum_cpt += ret.cpa[i,j,d,x]
                                sum_cpr += ret.cpa[i,j,d,x]
                                if int(ap.getOptionalArg("--verbose")) > 100:
                                    print "case 4:", i, j, d, x, ret.cpa[i,j,d,x]
                                continue
                            elif (x+1 - self.r[i] - d - self.r[j] < 0):
                                """CASE 3: TF i and j cannot both fit on the sequence."""
                                ret.cpa[i,j,d,x] = 0 # then P @ x = 0
                                sum_cpt += ret.cpa[i,j,d,x]
                                sum_cpr += ret.cpa[i,j,d,x]
                                if int(ap.getOptionalArg("--verbose")) > 100:
                                    print "case 3:", i, j, d, x, ret.cpa[i,j,d,x]
                                continue
                            """CASE 5: TFs i and j are different AND they can both fit on the URS..."""
                            ret.cpa[i,j,d,x] = rel_tf_expr[i] * pwm_tmp * self.gamma[j, i, d]
                            sum_cpt += ret.cpa[i,j,d,x]
                            sum_cpr += ret.cpa[i,j,d,x]
                            if int(ap.getOptionalArg("--verbose")) > 100:
                                print "case 5:", i, j, d, x, ret.cpa[i,j,d,x]
                            continue
                        elif (j == ap.params["numtr"]) or (i == j):
                            """CASE 6: there is no j:"""
                            ret.cpa[i,j,d,x] = rel_tf_expr[i] * pwm_tmp;
                            sum_cpt += ret.cpa[i,j,d,x]
                            sum_cpr += ret.cpa[i,j,d,x]
                            if int(ap.getOptionalArg("--verbose")) > 100:
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
        
        randp = random.uniform(0.0, ptables.cpr[site])
        #print "sample_cdf for site", site, "cpr sum =", ptables.cpr[site], "rand=", randp
        # now determine which TF and d value randp corresponds to.
        sump = 0.0
        reti = 0
        retj = 0
        retd = 0
        
        # debugging: as a test, just pick a random cell from the appropriate cpa
        # range, rather than picking an IID cell from the cpr-based range.
        #i = random.randint(0, ap.params["numtr"]-1)
        #j = random.randint(0, ap.params["numtr"])
        #d = random.randint(0, MAX_GD-1)
        
        
        for i in range(0, ap.params["numtr"]):
            reti = i
            for j in range(0, ap.params["numtr"]+1):
                retj = j
                for d in range(0, MAX_GD):
                    retd = d
                    sump += ptables.cpa[i,j,d,site]
                    #print "site", site, "tf", i, "tf", j, "d", d, "sump", sump, "randp", randp
                    if sump > randp:
                        break
                if sump > randp:
                    break
            if sump > randp:
                break
        return [reti, retj, retd]
        
        
    #def i_to_k(self, i):
    #    """Convert information bits to Kd"""
    #    return 1/(math.exp(INFO_ALPHA*i))
    
    def prob_expr(self, genome, ptables, gene, tf_expr_levels, ap, lam=None):
        """Returns a floating-point value corresponding to the expression level of gene,
        given the ProbTable ptables"""        
        #print "prob_expr, gene", gene.id, gene.urs
        
        pe_sum = 0
        min_r = min( self.r )
        configurations = {}
        """configurations: key = site, value = array of arrays, [i,j,d] samples"""
        for sample in range(0, IID_SAMPLES):
            
            """"1. build a configuration c_k, by sampling cells from ptables.cpa"""          
            this_config = {} # key = site, value = the TF bound starting at this site.
            site = 0
            while (site < gene.urs.__len__()):
                [i, j, d] = self.sample_cdf(site, ptables, ap)
                #print i, j, d
                if False == configurations.__contains__(site):
                    configurations[site] = []
                configurations[site].append( [i,j,d] )
                this_config[site] = i
                #print "site", site
                if i < ap.params["numtr"]-1:
                    site += genome.genes[i].pwm.P.__len__()
                else:
                    site += 1 # for the non-occuped case
                site += d
                
                if False == configurations.__contains__(site):
                    configurations[site] = []
                this_config[site] = j
                #print "site", site
                if j < ap.params["numtr"]-1:
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
            this_pe = (1/( 1+math.exp(-1*PE_SCALAR*(k_act-k_rep) ) ))
            pe_sum += this_pe
            #print "gene", gene.id, "k_act", k_act, "k_rep", k_rep, "this_pe", this_pe
        if ap.getOptionalArg("--verbose") > 2:
            self.print_configuration(configurations, genome, gene, ap)
        return (pe_sum / IID_SAMPLES)
    
    
    def print_configuration(self, configs, genome, gene, ap):
        foutpath = ap.getArg("--runid") + "/" + EXPR_PLOTS + "/gen" + self.gen_counter.__str__() + ".gid" + genome.id.__str__() + ".txt"
        fout = open( foutpath , "a")
        fout.write(". TIME " + self.t_counter.__str__() + " GENE " + gene.id.__str__() + "\t" + gene.urs + "\n")
        
        #fout.write("\n. URS binding at generation " + self.gen_counter.__str__() + " timeslice " + self.t_counter.__str__() + "\n" )        
        """configs is array of arrays, [site, tf_i, tf_j, distance between i and j]"""
        sites = configs.keys()
        sites.sort()
        for site in sites:
            tf_count = {}
            for c in configs[site]:
                if False == tf_count.__contains__( c[0] ):
                    tf_count[ c[0] ] = 0
                tf_count[ c[0] ] += 1.0 / configs[site].__len__()

            line = "site " + site.__str__()
            for tf in tf_count:
                pe = 0.0
                if tf < ap.params["numtr"]:
                    pe = genome.genes[tf].pwm.prob_binding(site, gene.urs)
                line += " \t" + tf.__str__() + " (%.3f"%tf_count[tf] + " %.3f)"%pe
            fout.write(line + "\n")
            #print line
        fout.write("\n")
        fout.close()

    
      