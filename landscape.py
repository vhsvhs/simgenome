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
        if rules != None:
            self.basal_gene_id = basal_gene_id
            self.rules = rules
        else:
            self.init_random()
    
    def init_random(self):
        last_time = 0
        for i in range(0, GOAL_COMPLEXITY):
            """First, expression > 0.5"""
            last_time += random.randint(1, GOAL_PULSE_MAX)
            e = 0.5
            rand_ri = random.randint(N_TR, N_REPORTER+N_TR-1)
            r = RULE_TYPES[0] 
            self.rules.append( Fitness_Rule(last_time, e, rand_ri, r) )
            """Next, expression at maximum"""     
            last_time += random.randint(1, GOAL_PULSE_MAX)
            e = 1.0
            r = RULE_TYPES[1]
            self.rules.append( Fitness_Rule(last_time, e, rand_ri, r) )
            """Finally, expression < 0.5"""
            last_time += random.randint(1, GOAL_PULSE_MAX)
            e = 0.5
            r = RULE_TYPES[2]
            self.rules.append( Fitness_Rule(last_time, e, rand_ri, r) )
    
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
        """Precalculates the cooperative/competitive binding interactions between all TFs."""
        # g is the binding coop term (see the coopfunc).
        # for now, g is all zeroes, so there is no cooperative binding.
        g = zeroes( (N_TR, N_TR), dtype=float)
        self.gamma = zeros( (N_TR, N_TR, MAX_GD), dtype=float)
        for i in range (0, N_TR):
            for j in range(0, N_TR):
                for d in range(0, MAX_GD):
                    self.gamma[i,j,d] = self.coopfunc( g[i,j], d)
    
    def get_fitness(self, genome):
        """Calculates the fitness of the given genome, over all time patterns in the fitness landscape.  Returns a floating-point value."""
        # Build the r vector...
        self.r = []
        self.maxr = 0
        for x in range(0, N_TR):
            self.r.append( genome.genes[x].pwm.P.__len__() )
            if self.r[x] > self.maxr:
                self.maxr = self.r[x]
        
        # Build the initial TF expression levels:
        # All TFs start with zero expression, except for those genes
        # that are activated by time patterns (i.e., they are the basal
        # signal-induced TF).  Basal TFs start with full expression.
        tf_expr_level = {}
        for tf in range(0, N_TR):
            tf_expr_level[ genome.genes[tf].id ] = 0.0
        for tp in self.timepatterns:
            tf_expr_level[tp.basal_gene_id] = 1.0

        # Get expression for each timeslice...
        # gene_expr[gene ID][timeslice] = expression level of gene at timeslice
        gene_expr = {} 
        for gene in genome.genes:
            gene_expr[gene.id] = []
        
        for timeslice in range(0, MAX_TIME):
            #
            # to-do: this would be a good place for simple parallelization
            #
            for gene in genome.genes:
                pe = self.get_expression(genome, gene, tf_expr_level)
                gene_expr[gene.id].append( pe )
            
            # Update TF expression levels for the next time iteration...
            for tf in range(0, N_TR):
                tf_expr_level[ genome.genes[tf].id ] = gene_expr[ genome.genes[tf].id ][timeslice]
            
            # debugging:
            print "debug landscape.py 143"
            for gene in genome.genes:
                print gene.id, gene_expr[ gene.id ][timeslice]

        #
        # to-do: assess fitness of genome, using gene_expr
        #
        
        # debugging:
        return 1.0
    
    def get_expression(self, genome, gene, tf_expr_level):
        """returns a floating-point value, the expression level of gene, given the TF expression levels"""
        pe = []
        ptables = ProbTable( N_TR, MAX_GD, gene.urs.__len__() )
        ptables = self.calc_prob_tables(genome, gene, tf_expr_level, ptables)
        return self.prob_expr(genome, ptables, gene)          
    
    def calc_prob_tables(self, genome, gene, rel_tf_expr, ret):
        """returns a ProbTable object. ret is the ProbTable that should be returned."""
        L = gene.urs.__len__()

        for x in range(0, L): # foreach site in gene's upstream region
            #print "debug calc_prob_tables, site", x
            for i in range(0, N_TR):   # foreach transcription factor
                
                # pwm_tmp: probability of TF i binding to the sequence with right edge at position x ?
                pwm_tmp = genome.genes[i].pwm.prob_binding( x, gene.urs )
                
                for j in range(0, N_TR+1): # +1 to also consider the empty case
                    for d in range(0, MAX_GD):                        
                        # CASE: TF i's PWM is too wide to start binding at site x
                        if (x+1 - self.r[i] < 0): 
                            ret.cpa[i,j,d,x] = 0 # P @ x = 0
                            continue
                        # CASE: TF i can start binding at site x
                        #    and TF i is identical to j.
                        if (x+1 - self.r[i] == 0 and j == 0 and d == 0): # // if TF i can bind here
                            ret.cpa[i,j,d,x] = rel_tf_expr[i] * pwm_tmp # basic case, no competition or cooperation
                            continue
                        if (j < N_TR): 
                            # CASE: TF i and j cannot both fit on the sequence
                            if (x+1 - self.r[i] - d - self.r[j-1] < 0):
                                ret.cpa[i,j,d,x] = 0 # then P @ x = 0
                                continue
                            # CASE: the distance between TFs i and j is too small.
                            if (d < MIN_TF_SEPARATION):
                                # we forbid TFs to bind this close together
                                ret.cpa[i,j,d,x] = 0
                                continue
                            # CASE: TFs i and j are different AND they can both fit on the URS...
                            ret.cpa[i,j,d,x] = rel_tf_expr[i] * pwm_tmp * self.gamma[j, i, d]
                            continue
                        # CASE: there is no j:
                        if j == N_TR:
                            ret.cpa[i,j,d,x] = rel_tf_expr[i] * pwm_tmp;
                            continue
    return ret
    
    #
    # continue here!
    #
    # also calculate cpt and cpr, but not cpm
    #

    
    def sample_cdf(self, site, ptables):
        """Picks a random configuration starting at site, drawing IID from ptables.
        ptables.cpr must contain the cummulative marginal distributions.
        This method returns [i,j,d], where TF is binding to site, followed by distance
        d, and then TF j binds."""
        
        randp = random.uniform(0.0, ptables.cpr[site])
        # now determine which TF and d value randp corresponds to.
        sump = 0.0
        i = 0
        j = 0
        d = 0
        for i in range(0, N_TR):
            for j in range(0, N_TR+1):
                for d in range(0, MAX_GD):
                    sump += ptables.cpa[i,j,d,site]
                    if sump > randp:
                        break
        # at this point, we've chosen TF i binding to site, followed by distance d, and then TF j
        return [i, j, d]
        
    
    def prob_expr(self, genome, ptables, gene, print_configs=False, lam=None,):
        """Returns a floating-point value corresponding to the expression level of gene,
        given the ProbTable ptables"""        
        # lam[i] is the maximum expression level of gene i
        if lam == None:
            lam = []
            for tf in range(0, N_TR + N_REPORTER):
                lam.append(1.0)
        
        L = gene.urs.__len__()
        pe_tmp = 0
        min_r = min( self.r )

        
        for sample in range(0, MAX_GA_GENS):
            # 1. build a configuration c_k, by sampling cells from ptables.cpa            
            this_config = {} # key = TF id, value = site at which TF binds.  Not all TFs are necessarily in this_config.
            site = 0
            while (site < L):
                [i, j, d] = sample_cdf()
                this_config[i] = site
                site += genome.genes[i].pwm.P.__len__()
                site += d
                this_config[j] = site
                site += genome.genes[j].pwm.P.__len__()
            
            #
            # 2. Then calculate P(E|c_k)
            #  . Sum the lambda values for all the TFs in the configuration.
            sum_lambda = 0.0
            for tf in this_config:
                sum_lamba += lam[tf]
            
            
            # 3. the probability of this configuration equals:
            pe_tmp = (1/(1+math.exp(-1*sum_lambda)))
            #
            #
        return pe_tmp/MAX_GA_GENS
    
    
      