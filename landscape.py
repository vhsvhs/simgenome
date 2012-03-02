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
    
    def __init__(self, basal_gene_id, rules = None):
        print "\n. Building a new fitness landscape..."
        self.basal_gene_id = basal_gene_id
        if rules != None:
            self.rules = rules
        else:
            self.init_random()
    
    def init_random(self):
        last_time = 0
        for i in range(0, GOAL_COMPLEXITY):
            """First, expression > 0.5"""
            last_time += random.randint(1, GOAL_PULSE_MAX)
            if last_time > MAX_TIME:
                last_time = MAX_TIME
            e = 0.5
            rand_ri = random.randint(N_TR, N_REPORTER+N_TR-1)
            r = RULE_TYPES[0] 
            self.rules.append( Fitness_Rule(last_time, e, rand_ri, r) )
            """Next, expression at maximum"""     
            last_time += random.randint(1, GOAL_PULSE_MAX)
            e = 1.0
            r = RULE_TYPES[1]
            if last_time > MAX_TIME:
                last_time = MAX_TIME
            self.rules.append( Fitness_Rule(last_time, e, rand_ri, r) )
            """Finally, expression < 0.5"""
            last_time += random.randint(1, GOAL_PULSE_MAX)
            e = 0.5
            r = RULE_TYPES[2]
            if last_time > MAX_TIME:
                last_time = MAX_TIME
            self.rules.append( Fitness_Rule(last_time, e, rand_ri, r) )
    
    def __str__(self):
        ret = ""
        ret += "Time pattern for basal gene " + self.basal_gene_id.__str__() + " " + self.rules.__len__().__str__() + " rules."
        return ret
    
class Landscape:
    """An array of Time_Patterns"""
    timepatterns = []
    gamma = None
    r = None
    
    def __init__(self, scapepath = None, genome = None):
        """scapepath is a filepath."""
        if genome != None:
            self.init_random(genome)
        elif scapefile != None:
            self.init_spec(scapefile)
        self.set_gamma()

    def init_random(self, genome):
        """Creates a random fitness goal, using genome as the seed."""    
        """First, select a basal TR to activate the goal."""
        randid = random.randint(0, N_TR-1)
        #print "basal id = ", randid
        
        """Next, build a random pattern.
        For now, the random Landscape contains only one time pattern."""
        self.timepatterns.append( Time_Pattern(randid) )
    
    def init_spec(self, scapepath):
        """The landscape specification file should be formatted as follows:
            time expression_level reporter_id rule
            on each line."""
        fin = open(scapepath, "r")
        for l in fin.readlines():
            tokens = l.split()
            if tokens.__len__() < 4:
                print "Hmmm, an error occurred. A line in your landscape specification file is formatted incorrectly."
                print "This line:\n", l, "\n"
                exit()
        fin.close()

    def coopfunc(self, g, d):
        """Calculates the degree of binding cooperativity between two TFs binding distance d apart.
        g is positive for synergistic interactions, and negative for antagonistic interactions.
        g ranges from -1 to +inifinity."""
        return 1 + g * math.exp( (-1)*(d**2)/V_RATE_OF_COOP_DECAY );

    def set_gamma(self):
        """Prec-calculates the cooperative/competitive binding interactions between all TFs."""
        # g is the binding coop term (see the coopfunc).
        # for now, g is all zeroes, so there is no cooperative binding.
        #
        # to-do: grab gamma values from the command-line, or radomly initialize them
        # from an a priori distribution.
        g = zeros( (N_TR, N_TR), dtype=float)
        self.gamma = zeros( (N_TR, N_TR, MAX_GD), dtype=float)
        for i in range (0, N_TR):
            for j in range(0, N_TR):
                for d in range(0, MAX_GD):
                    self.gamma[i,j,d] = self.coopfunc( g[i,j], d)
    
    def fitness_helper(self, gene_expr_level):
        fitsum = 0.0
        for pattern in self.timepatterns:
            for rule in pattern.rules:
                obs_expr = gene_expr_level[rule.reporter_id][rule.timepoint]
                if rule.rule_type(obs_expr, rule.expression_level):
                    fitsum += 1.0
                else:
                    fitsum += 1.0 / abs(obs_expr - rule.expression_level)
        return fitsum
                    
    
    def get_fitness(self, genome):
        """Calculates the fitness of the given genome, over all time patterns in the fitness landscape.  Returns a floating-point value."""
        # Build the r vector...
        self.r = []
        self.maxr = 0
        for x in range(0, N_TR):
            self.r.append( genome.genes[x].pwm.P.__len__() )
            if self.r[x] > self.maxr:
                self.maxr = self.r[x]
        
            
        # gene_expr[gene ID][timeslice] = expression level of gene at timeslice
        gene_expr = {}
        last_gene_expr = {} 
        for gene in genome.genes:
            gene_expr[gene.id] = [MINIMUM_ACTIVITY_LEVEL] # all genes begin with zero
            last_gene_expr[gene.id] = [MINIMUM_ACTIVITY_LEVEL] # all genes begin with zero
            
        # tf_expr_level[gene ID] = current expression level
        # tf_expr_level is used as short-term variable inside the loop "for timeslice..."
        tf_expr_level = {}
        for timepattern in self.timepatterns:
            gene_expr[timepattern.basal_gene_id] = [MAXIMUM_ACTIVITY_LEVEL]
                
        for timeslice in range(1, MAX_TIME):

            # 1a. Ensure that basal genes remain activated
            for timepattern in self.timepatterns:
                if gene_expr[timepattern.basal_gene_id][timeslice-1] < MAXIMUM_ACTIVITY_LEVEL:
                    gene_expr[timepattern.basal_gene_id][timeslice-1] = MAXIMUM_ACTIVITY_LEVEL
            
            # for debugging:
            #if timeslice > 10:
            #    for timepattern in self.timepatterns:
            #        gene_expr[timepattern.basal_gene_id].append(MINIMUM_ACTIVITY_LEVEL)               
            
            
            # 1b. Update TF expression levels...            
            for tf in range(0, N_TR):
                tf_expr_level[ genome.genes[tf].id ] = gene_expr[ genome.genes[tf].id ][timeslice-1]
                

            print "\n+++++++++++++++\nTIME", timeslice
            
            #
            # to-do: this would be a good place for simple parallelization
            #
            
            for gene in genome.genes:
                # calculate the delta G of binding on the cis-region for every gene.
                pe = self.get_expression(genome, gene, tf_expr_level)
                
                # if the delta G is high enough, then transcribe the gene.
                if pe >= ACTIVATION_THRESHOLD:
                    new_expr_level = gene_expr[gene.id][timeslice-1] * GROWTH_FACTOR;
                    gene_expr[gene.id].append( new_expr_level )
                else: # otherwise, decay the expression level of the gene.
                    new_expr_level = gene_expr[gene.id][timeslice-1] / DECAY_FACTOR;
                    gene_expr[gene.id].append( new_expr_level )
            
                # sanity check:
                if gene_expr[gene.id][timeslice] > MAXIMUM_ACTIVITY_LEVEL:
                    gene_expr[gene.id][timeslice] = MAXIMUM_ACTIVITY_LEVEL
                if gene_expr[gene.id][timeslice] < MINIMUM_ACTIVITY_LEVEL:
                    gene_expr[gene.id][timeslice] = MINIMUM_ACTIVITY_LEVEL
            
            # debugging:
            #for gene in genome.genes:
                print "gene", gene.id, "expr= %.3f"%gene_expr[ gene.id ][timeslice], "pe=%.3f"%pe
            #    if gene.has_dbd:
            #        print gene.pwm
            

            #
            # to-do: if gene expression has not changed from the last timeslice
            #    then we've reached equilibrium, so stop cycling through time slices.
            #


        #
        # to-do: assess fitness of genome, using gene_expr
        #
        return self.fitness_helper(gene_expr)
    
    def get_expression(self, genome, gene, tf_expr_levels):
        """returns a floating-point value, the expression level of gene, given the TF expression levels"""
        pe = []        
        ptables = ProbTable( N_TR, MAX_GD, gene.urs.__len__() )
        ptables = self.calc_prob_tables(genome, gene, tf_expr_levels, ptables)        
        return self.prob_expr(genome, ptables, gene, tf_expr_levels)          
    
    def calc_prob_tables(self, genome, gene, rel_tf_expr, ret):
        """returns a ProbTable object. ret is the ProbTable that should be returned."""
        L = gene.urs.__len__()        
        for x in range(0, L): # foreach site in gene's upstream region
            sum_cpr = 0.0
            for i in range(0, N_TR):   # foreach transcription factor
                pwm_tmp = genome.genes[i].pwm.prob_binding( x, gene.urs )
                #print "TF", i, "binds", gene.urs, "at site", x, "with %.3f"%pwm_tmp, "units."
                sum_cpt = 0.0
                for j in range(0, N_TR+1): # +1 to also consider the empty case
                    for d in range(0, MAX_GD):                        
                        # CASE 1: TF i's PWM is too wide to start binding at site x
                        if (L-x < self.r[i]): 
                            ret.cpa[i,j,d,x] = 0 # P @ x = 0
                            sum_cpt += ret.cpa[i,j,d,x]
                            sum_cpr += ret.cpa[i,j,d,x]
                            #print "case 1:", i, j, d, x, ret.cpa[i,j,d,x], self.r[i]
                            continue
                        # CASE 2: TF i can start binding at site x
                        #    and TF i is identical to j.
                        elif (L-x >= self.r[i] and (j == i or j == N_TR) and d == 0): # // if TF i can bind here
                            ret.cpa[i,j,d,x] = rel_tf_expr[i] * pwm_tmp # basic case, no competition or cooperation
                            sum_cpt += ret.cpa[i,j,d,x]
                            sum_cpr += ret.cpa[i,j,d,x]
                            #print "case 2:", i, j, d, x, ret.cpa[i,j,d,x]
                            continue
                        elif (j < N_TR): # else, cooperative/competitive binding...
                            # CASE 4: the distance between TFs i and j is too small.
                            if (d < MIN_TF_SEPARATION):
                                # we forbid TFs to bind this close together
                                ret.cpa[i,j,d,x] = 0
                                sum_cpt += ret.cpa[i,j,d,x]
                                sum_cpr += ret.cpa[i,j,d,x]
                                #print "case 4:", i, j, d, x, ret.cpa[i,j,d,x]
                                continue
                            # CASE 3: TF i and j cannot both fit on the sequence
                            elif (x+1 - self.r[i] - d - self.r[j-1] < 0):
                                ret.cpa[i,j,d,x] = 0 # then P @ x = 0
                                sum_cpt += ret.cpa[i,j,d,x]
                                sum_cpr += ret.cpa[i,j,d,x]
                                #print "case 3:", i, j, d, x, ret.cpa[i,j,d,x]
                                continue
                            # CASE 5: TFs i and j are different AND they can both fit on the URS...
                            ret.cpa[i,j,d,x] = rel_tf_expr[i] * pwm_tmp * self.gamma[j, i, d]
                            sum_cpt += ret.cpa[i,j,d,x]
                            sum_cpr += ret.cpa[i,j,d,x]
                            #print "case 5:", i, j, d, x, ret.cpa[i,j,d,x]
                            continue
                        # CASE 6: there is no j:
                        #elif (j == N_TR) or (i == j):
                        #    ret.cpa[i,j,d,x] = rel_tf_expr[i] * pwm_tmp;
                        #   sum_cpt += ret.cpa[i,j,d,x]
                        #   sum_cpr += ret.cpa[i,j,d,x]
                        #   #print "case 6:", i, j, d, x, ret.cpa[i,j,d,x]
                        #   continue
                ret.cpt[i,x] = sum_cpt
            ret.cpr[x] = sum_cpr
        return ret
    
    def sample_cdf(self, site, ptables):
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
        #i = random.randint(0, N_TR-1)
        #j = random.randint(0, N_TR)
        #d = random.randint(0, MAX_GD-1)
        
        
        for i in range(0, N_TR):
            reti = i
            for j in range(0, N_TR+1):
                retj = j
                for d in range(0, MAX_GD):
                    retd = d
                    sump += ptables.cpa[i,j,d,site]
                    #print "tf", i, "tf", j, "d", d, "sump", sump
                    if sump > randp:
                        break
                if sump > randp:
                    break
            if sump > randp:
                break
        return [reti, retj, retd]
        
        
    def i_to_k(self, i):
        """Convert information bits to Kd"""
        return 1/(math.exp(INFO_ALPHA*i))
    
    def prob_expr(self, genome, ptables, gene, tf_expr_levels, print_configs=False, lam=None):
        """Returns a floating-point value corresponding to the expression level of gene,
        given the ProbTable ptables"""        
        #print "prob_expr, gene", gene.id, gene.urs
        
        pe_tmp = 0
        min_r = min( self.r )
        configurations = {} # key = site, value = array of arrays, [i,j,d] samples
        for sample in range(0, IID_SAMPLES):
            
            # 1. build a configuration c_k, by sampling cells from ptables.cpa            
            this_config = {} # key = site, value = the TF bound starting at this site.
            site = 0
            while (site < gene.urs.__len__()):
                [i, j, d] = self.sample_cdf(site, ptables)
                if False == configurations.__contains__(site):
                    configurations[site] = []
                configurations[site].append( [i,j,d] )
                if i < N_TR:
                    this_config[site] = i
                    site += genome.genes[i].pwm.P.__len__()
                site += d
                if j < N_TR:
                    this_config[site] = j
                    site += genome.genes[j].pwm.P.__len__()
                
            # 2. Calculate the binding energy of the configuration.
            sum_lambda_act = 0.0
            sum_lambda_rep = 0.0
            for site in this_config:
                tf = this_config[site]
                
                # Get the strength of TF binding at this site
                tf_specificity = genome.genes[tf].pwm.prob_binding(site, gene.urs)
                #print "site", site, "tf_specificity", tf_specificity
                if genome.genes[tf].is_repressor:
                    sum_lambda_rep += tf_specificity
                else:
                    sum_lambda_act += tf_specificity
            
            # 3. the probability of this configuration equals:
            # this is what Kevin did:
            #pe_tmp = (1/(1+math.exp(-1*sum_lambda)))
            # but this is what I'm doing instead:
            # This incorporates the hill equation
            
            #k_act = self.i_to_k(sum_lambda_act)
            #k_rep = self.i_to_k(sum_lambda_rep)
            k_act = sum_lambda_act
            k_rep = sum_lambda_rep
            this_pe = (1/( 1+50*math.exp(-1*DELTA_G_SCALAR*(k_act-k_rep) ) ))
            pe_tmp += this_pe
            #print "k_act", k_act, "k_rep", k_rep, "pe_tmp", pe_tmp
        #self.print_configuration(configurations, genome, gene)
        return (pe_tmp / IID_SAMPLES)
    
    
    def print_configuration(self, configs, genome, gene):
        """configs is array of arrays, [site, tf_i, tf_j, distance between i and j]"""
        sites = configs.keys()
        sites.sort()
        for site in sites:
            tf_count = {}
            for c in configs[site]:
                if False == tf_count.__contains__( c[0] ):
                    tf_count[ c[0] ] = 0
                tf_count[ c[0] ] += 1.0 / configs[site].__len__()

            line = "site: " + site.__str__()
            for tf in tf_count:
                line += " \ttf: " + tf.__str__() + " %.3f"%tf_count[tf]
            print line

    
      