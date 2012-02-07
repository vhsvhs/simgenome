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
        print "Building a new fitness landscape..."
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
    r = None
    
    def __init__(self, scapepath = None, genome = None):
        """scapepath is a filepath."""
        if genome != None:
            self.init_random(genome)
        elif scapefile != None:
            self.init_spec(scapefile)
            
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
    
    def get_fitness(self, genome):
        """Calculates the fitness of the given genome.  Returns a floating-point value."""
        
        print "Calculating fitness of genome", genome.id
        
        # build the r vector
        self.r = []
        self.maxr = 0
        for x in range(0, N_TR):
            self.r.append( genome.genes[x].pwm.P.__len__() )
            if self.r[x] > self.maxr:
                self.maxr = self.r[x]
        
        gene_expr = {}
        for gene in genome.genes:
             pe = self.get_expression(genome, gene)
             #gene_expr[ gene.id ] = pe
        
        #
        # to-do: asses fitness...
        #
        
        
        # debugging:
        return pe
    
    def get_expression(self, genome, gene):
        """ Returns an array of expression values."""
        
        #print "Calculating expression for gene", gene.id, " in genome", genome.id
        
        #if (pe.size() != conditions.size())
        #throw SimError("num. of desired expression levels must match num. cond.");
        # loop through the conditions first pre-computing tables then est. expr. */
        #ptables = ProbTable( N_TR, MAX_GD, gene.urs.__len__() )
        #for i in range(0, conditions.size() ):
            #ptables = self.calc_prob_tables(genome, conditions[i], alpha, gamma, gene, ptables)
            #pe.append( self.prob_expr(genome, ptables, gene, i) )
        #return pe
        return
    
    def calc_prob_tables(self, genome, rel_tf_expr, alpha, gm, gene, ret):
        """returns a ProbTable object"""
        L = gene.urs.__len__()

        for x in range(0, L): # foreach site in gene's upstream region
            for i in range(0, N_TR):   # foreach transcription factor
                pwm_tmp = p_pwm(i, x, gene.urs)
                for j in range(0, N_TR+1):
                    for d in range(0, MAX_GD):
                        if (x+1 - self.r[i] < 0):
                            ret.cpa.mem[ ret.cpa.pos([i,j,d,x]) ] = 0 # then probability = 0
                            continue
                        if (x+1 - self.r[i] == 0 and j == 0 and d == 0): # // if TF i can bind here
                            ret.cpa.mem[ ret.cpa.pos([i,j,d,x]) ] = alpha[i] * rel_tf_expr[i] * pwm_tmp
                            continue
                        if (j > 0):
                            if (x+1 - self.r[i] - d - self.r[j-1] < 0):
                                ret.cpa.mem[ ret.cpa.pos([i,j,d,x]) ] = 0
                                continue
                            if (d < minimum_tf_separation):
                                # we forbid TFs this close together */
                                ret.cpa.mem[ ret.cpa.pos( [i,j,d,x] ) ] = 0 
                                continue
                            ret.cpa.mem[ ret.cpa.pos( [i,j,d,x] ) ] = ret.cpt.mem[ ret.cpt.pos([j-1,x-r[i]-d] ) ]  * gm.mem[ gm.pos([j-1, i, d]) ] * alpha[i] * (rel_tf_expr)[i] * pwm_tmp
                            continue
                        else:  # /* previous TF is D or further away or no other TF */
                            if (d > 0): # /* we don't use this slot */
                                ret.cpa.mem[ ret.cpa.pos([i,0,d,x]) ]= 0 
                                continue
                            if (x+1 - self.r[i] - MAX_GD <= 0): # /* no other TF */
                                ret.cpa.mem[ ret.cpa.pos([i,0,0,x]) ] = alpha[i] * (rel_tf_expr)[i] * pwm_tmp;
                                continue
                            ret.cpa.mem[ ret.cpa.pos([i,0,0,x])] = (1+ret.cpr.mem[ ret.cpr.pos([ x - self.r[i]-MAX_GD ]) ]) * alpha[i] * (rel_tf_expr)[i] * pwm_tmp;
                            continue
                ans1 = ret.cpa.mem[x*N_TR*(N_TR+1)*MAX_GD+i : (N_TR+1)*MAX_GD : N_TR] # // collect elements from this slice of cpa
                ret.cpt.mem[ ret.cpt.pos([i,x])] = sum(ans1) # // ... and sum those elements into cpt
            if (x == 0): 
                ans2 = ret.cpt.mem[0 : N_TR : 1]
                ret.cpt.mem[0] = sum(ans2)
            else: 
                ans3 = ret.cpt.mem[x*N_TR : N_TR : 1]
                ret.cpr.mem[ x ] = ret.cpr.mem[ x-1 ] + sum(ans3)
    
       # /* make the margin tables */
        for k1 in range(0, N_TR):
            for k2 in range(0, L):
                csum = 0
                ret.cpm.mem[ ret.cpm.pos([0, k1, k2]) ] = 0
                for k3 in range(0, (int)((N_TR+1)*MAX_GD) ): # // foreach slice in cpa
                    csum += ret.cpa.mem[k2*N_TR*(N_TR+1)*MAX_GD+k1+N_TR*k3]
                    ret.cpm.mem[((N_TR+1)*MAX_GD+1)*(k2*N_TR+k1)+1+k3] = csum
    
       # /* now make s_array cumulative, we need the first two values to be 0, 1
       #  * and shifting the result of the cumulative sum over by 2 */
        tmp1 = 1 
        tmp2 = ret.cpa.mem[0];
        tmp3 = ret.cpa.mem[1];
        ret.cpa.posi[0] = 0;
        for i in range(2, reduce(mul, ret.cpa.dims) ):
            ret.cpa.mem[i-1] = tmp1
            tmp1 += tmp2
            tmp2 = tmp3
            tmp3 = ret.cpa.mem[i]
        ret.cpa.mem[i-1] = tmp1
        tmp1 += tmp2
        ret.cpa.mem[i] = tmp1
        tmp1 += tmp3
        ret.cpa.mem[i+1] = tmp1
    
       # /* now we add one to s_row to account for the empty configuration */
        for i in range(0, ret.cpr.mem.__len__()):
            ret.cpr.mem[i] += 1
    
        return ret
    
    def prob_expr(self, genome, L, M, D, r, lam, ptables, logistic_parameter, print_configs):
        """ L = """ 
        #int up, min_r, cell, left_tf, right_tf, d, tbl_size, pos;
        #int i=steps;
        pe_tmp = 0
        min_r = min( self.r )
        mmp1d = M * (M+1) * D
        mmp1 = M * (M+1)
        mp1d = (M+1) * D
        
        for sample in range(0, MAX_GA_GENS):
            sum_lam = lam[0]
            left_tf = 0
            up = L
            d = 0
            while (up >= min_r):
                if (left_tf == 0):
                    tbl_size = mmp1d * up + 1 # size of ptables.cpa + 1 (for empty configuration)
                            #it's called like this: sample_cdf(int tbl_size, double tbl_sum, valarray<double> &tbl, int offset)
                    cell = sample_cdf(tbl_size, ptables.cpr.mem[up-1], ptables.cpa.mem, 0) # to-do: implement sample_cdf
                    if (cell == 0): break
                    right_tf = ((cell-1) % M)+1
                    pos = ((cell-1) / mmp1d)+1
                    left_tf = ((cell-1) % mmp1) / M
                    d = ((cell-1) % mmp1d) / mmp1
                else: 
                    cell = sample_cdf(mp1d, ptables.cpt.mem[ ptables.cpt.pos([left_tf-1, up-1] ) ], ptables.cpm.mem, ((up-1)*M + left_tf-1) * (mp1d+1))
                    right_tf = left_tf
                    pos = up
                    left_tf = cell % (M+1)
                    d = cell / (M+1)
                if (left_tf == 0): 
                    up = pos - self.r[right_tf-1] - D
                else: 
                    up = pos - self.r[right_tf-1] - d
                sum_lam += lam[right_tf]
                if (print_configs):
                    print right_tf, ":", pos, ", "
            sum_lam = sum_lam/L*1000
            if (print_configs):
                print "\n"
            pe_tmp += (1/(1+exp((-1)*sum_lam/logistic_parameter)))
        return pe_tmp/MAX_GA_GENS
          
          