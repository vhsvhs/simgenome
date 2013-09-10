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
    pe_sum = 0.0   
    lock = threading.Lock()
    temp_genome = None
 
    def __init__(self, ap):
        self.rulecollections = {} # array of Rule_Collection objects
        self.r = None
        self.t_counter = ap.params["generation"]
    
    def init(self, ap, tp=None):
        tp = get_input_rules_from_file(ap)
        
        """Verify that we found input patterns and fitness rules."""
        if tp == None:
            if comm.Get_rank() == 0:
                print "\n. Hmmm, I didn't find any fitness rules in the file", ap.getOptionalArg("--patternpath")
            exit()
        
        if comm.Get_rank() == 0:
            if ap.params["verbosity"] > 0:
                print "\n. Building the fitness landscape described in", ap.getOptionalArg("--patternpath")
        self.rulecollections = tp
        
        #if comm.Get_rank() == 0:
        #    for t in self.rulecollections:
        #        print t
        
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
                expr_error += abs(obs_expr - rule.expression_level) * rule.multiplier
            max_expr_error += 1.0 * rule.multiplier
        expr_error = expr_error / max_expr_error # Normalize the expression error by the max. possible error
        return math.exp(FITNESS_SCALAR * expr_error) # FITNESS_SCALAR is defined in configuration.py
        
        
    def get_fitness_thread(self, gene, ap, tf_expr_level, ko, timeslice):
        if gene.id in ko:
            pe = 0.0
        else:
            """ Calculate the delta G of binding on the cis-region for every gene.
                See the notes in the code for the function "get_expr_modifier".
                pe ranges from -0.5 (repression) to 0.5 (strong activation)"""
            pe = self.get_expr_modifier(self.temp_genome, gene, tf_expr_level, ap)                
                            
        if pe > 0.0:
            expr_modifier = ap.params["growth_rate"] * pe
        elif pe < 0.0:
            expr_modifier = ap.params["decay_rate"] * pe
        else:
            expr_modifier = 0.0 # this will result in no change to expression level
        
        """New expression level equals the old expression level, plus the modification (which may be positive or negative)"""
        #print self.temp_genome
        self.lock.acquire()
        new_expr_level = self.temp_genome.gene_expr[gene.id][timeslice] + expr_modifier;
        self.temp_genome.gene_expr[gene.id].append( new_expr_level )
        
        """ Prevent out-of-bounds expression levels."""
        if self.temp_genome.gene_expr[gene.id][timeslice+1] > MAXIMUM_ACTIVITY_LEVEL:
            self.temp_genome.gene_expr[gene.id][timeslice+1] = MAXIMUM_ACTIVITY_LEVEL
        if self.temp_genome.gene_expr[gene.id][timeslice+1] < MINIMUM_ACTIVITY_LEVEL:
            self.temp_genome.gene_expr[gene.id][timeslice+1] = MINIMUM_ACTIVITY_LEVEL
        self.lock.release()
                    
    def get_fitness(self, genome, ap, ko=[]):

        """Calculates the fitness of the given genome, over all time patterns in the fitness landscape.
        Any gene IDs specified in the list 'ko' will be set to MINIMUM_ACTIVITY_LEVEL at each timeslice."""
                
        if ap.params["verbosity"] >= 99:
            notime = 0.0
            timestart = datetime.utcnow()
                
        # First, build the r vector...
        self.r = []
        self.maxr = 0
        for x in ap.params["trlist"]:
            #print "landscape.py 86", x, genome.genes[x], genome.genes[x].pwm
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
                    
            for timeslice in xrange(0, ap.params["maxtime"]):
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
                                print "r:", rid, "\tgen:", ap.params["generation"], "\ttime 0", "\tID", genome.id, "\tgene", gene.id, marka, "\texpr: %.3f"%genome.gene_expr[ gene.id ][timeslice], "\t\t[" + gene.name + "]"
                        print ""
    
                
                """ Update all TF expression levels."""         
                for tf in ap.params["trlist"]:
                    tf_expr_level[ genome.genes[tf].id ] = genome.gene_expr[ genome.genes[tf].id ][timeslice]
                            
                """ For each gene, update its expression level based on the delta G of binding at its regulatory sequence."""               
                threads = []
                
                for gene in genome.genes:            
                    #
                    #
                    #
                    #self.temp_genome = genome
                    #t = threading.Thread(target = self.get_fitness_thread, args=(gene, ap, tf_expr_level, ko, timeslice))
                    #threads.append( t )
                    #t.start()

                    
                    if gene.id in ko:
                        pe = 0.0
                    else:
                        """ Calculate the delta G of binding on the cis-region for every gene.
                            See the notes in the code for the function "get_expr_modifier".
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
                
                for t in threads:
                    t.join()
          
                #genome = self.temp_genome
                for gene in genome.genes:
                    """ Print a report to stdout."""
                    if ap.params["verbosity"] > 5:                    
                        expr_delta = 0.0
                        print genome.gene_expr[ gene.id ]
                        expr_delta = genome.gene_expr[ gene.id ][timeslice+1] - genome.gene_expr[ gene.id ][timeslice]
                        marka = "      "
                        if gene.has_dbd and gene.is_repressor == False:
                            marka = "[act.]"
                        elif gene.has_dbd and gene.is_repressor:
                            marka = "[rep.]"
                        print "r:", rid, "\tgen:", ap.params["generation"], "\ttime", timeslice+1, "\tID", genome.id, "\tgene", gene.id, marka, "\texpr: %.3f"%genome.gene_expr[ gene.id ][timeslice+1], "\td: %.3f"%expr_delta, "\t[" + gene.name + "]" 
    
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
        
        #
        # September 2013: begin ctypes revision here
        #
        
        # 1. Instantiate ptables
        # ptables is a multidimensional matrix that holds the occupancy distribution
        timeon = datetime.utcnow()
        ptables = ProbTable( ap.params["numtr"], ap.params["maxgd"], gene.urs.__len__() )
        ap.params["sumtime_ptables"] += (datetime.utcnow() - timeon).total_seconds()
        
        # 2. Initialize ptables
        # Here we fill ptables with values, using a dynamic algorithm similar to the knapsack problem
        timeon = datetime.utcnow()
        ptables = self.calc_prob_tables(genome, gene, tf_expr_levels, ptables, ap)          
        ap.params["sumtime_calcprobtables"] += (datetime.utcnow() - timeon).total_seconds()  
        
        # 3. Use ptables
        # Here we sample from the ptables distribution:
        timeon = datetime.utcnow()
        pe = self.prob_expr(genome, ptables, gene, tf_expr_levels, ap)        
        ap.params["sumtime_probexpr"] += (datetime.utcnow() - timeon).total_seconds()
                
        # Note: pe is the probability of increasing the expression of the gene. 
        # It ranges from 0.0 to 1.0.  A value of 0.5 means the gene expression will remain unchanged.
        pe = pe - 0.5 # This will make pe range from -0.5 to +0.5.  
        if ap.params["verbosity"] >= 99:
            print "landscape.py 209 pe before/after", (pe + 0.5), pe
        return pe
    
    def calc_prob_tables(self, genome, gene, rel_tf_expr, ret, ap):
        """returns a ProbTable object, named ret."""
        # M is num TR
        M = ap.params["numtr"]
        D = ap.params["maxgd"]
        L = gene.urs.__len__()
        dim1 = D*M*(M+1)
        dim2 = D*(M+1)
        dim3 = D
        
        # D is max gamma distance
        # L is URS length       
        for x in xrange(0, L): # foreach site in gene's upstream region
            sum_cpr = 0.0
            for i in ap.params["trlist"]:   # foreach transcription factor
                #print genome.genes[i].gamma
                genei = genome.genes[i]
                seq = gene.urs[ x: ]
                this_aff = genei.pwm.affinity(gene.urs[x:x+genei.pwm.P.__len__()])
                #print "TF", i, "binds", gene.urs, "at site", x, "with %.3f"%this_aff, "this_aff."
                sum_cpt = 0.0
                for j in ap.params["trlist+"]: # +1 to also consider the empty case
                    for d in ap.params["rangegd"]:
                        #print x, i, j, d                        
                        if (L-x < self.r[i]): 
                            """ CASE 1: TF i's PWM is too wide to start binding at site x"""
                            ret.cpa[x*dim1 + i*dim2 + j*dim3 + d] = 0 # P @ x = 0
                            sum_cpt += ret.cpa[x*dim1 + i*dim2 + j*dim3 + d]
                            sum_cpr += ret.cpa[x*dim1 + i*dim2 + j*dim3 + d]
                            if ap.params["verbosity"] > 100:
                                print "case 1:", "TF", i, " and TF", j, "distance", d, "site", (x+1), "p= ", ret.cpa[x*dim1 + i*dim2 + j*dim3 + d]
                            continue
                        if (j < ap.params["numtr"]):  # TF j is real, not the empty slot.
                            if (d < MIN_TF_SEPARATION): 
                                # If d == 0, then this case should never be reached.
                                """ CASE 4: the distance between TFs i and j is too small."""
                                # we forbid TFs to bind this close together
                                ret.cpa[x*dim1 + i*dim2 + j*dim3 + d] = 0 # then P @ x = 0
                                if ap.params["verbosity"] > 100:
                                    print "case 4:", "TF", i, " and TF", j, "distance", d, "site", (x+1), "p= ", ret.cpa[x*dim1 + i*dim2 + j*dim3 + d]
                                continue
                            elif (L - self.r[i] - d - self.r[j] < 0):
                                """CASE 3: TF i and j cannot both fit on the sequence."""
                                ret.cpa[x*dim1 + i*dim2 + j*dim3 + d] = 0 # then P @ x = 0
                                if ap.params["verbosity"] > 100:
                                    print "case 3:", "TF", i, " and TF", j, "distance", d, "site", (x+1), "p= ", ret.cpa[x*dim1 + i*dim2 + j*dim3 + d]
                                continue
                            else:                                
                                """Case 5: TF i and TF j can fit with distance d between them"""
                                #print rel_tf_expr[i] * this_aff * ( (genome.genes[i].gamma[j, d] + genome.genes[j].gamma[i,d]) * 0.5 )
                                ret.cpa[x*dim1 + i*dim2 + j*dim3 + d] = rel_tf_expr[i] * this_aff * ( (genome.genes[i].gamma[j, d] + genome.genes[j].gamma[i,d]) * 0.5 )
                                sum_cpt += ret.cpa[x*dim1 + i*dim2 + j*dim3 + d]
                                sum_cpr += ret.cpa[x*dim1 + i*dim2 + j*dim3 + d]
                                if ap.params["verbosity"] > 100:
                                    print "case 5:", "TF", i, " and TF", j, "distance", d, "site", (x+1), "p= ", ret.cpa[x*dim1 + i*dim2 + j*dim3 + d]
                                continue
                        elif (j == ap.params["numtr"]): # j is the empty slot
                            """CASE 6: no binding here:"""
                            ret.cpa[x*dim1 + i*dim2 + j*dim3 + d] = rel_tf_expr[i] * this_aff;
                            sum_cpt += ret.cpa[x*dim1 + i*dim2 + j*dim3 + d]
                            sum_cpr += ret.cpa[x*dim1 + i*dim2 + j*dim3 + d]
                            if ap.params["verbosity"] > 100:
                                print "case 6:", "TF", i, " and TF", j, "distance", d, "site", (x+1), "p= ", ret.cpa[x*dim1 + i*dim2 + j*dim3 + d]
                            continue
                ret.cpt[i,x] = sum_cpt
            ret.cpr[x] = sum_cpr
        #if ap.params["verbosity"] > 99:
        #    print ret
        return ret
        
    def sample_cdf(self, site, ptables, ap):
        """Picks a random configuration starting at site, drawing IID samples from ptables.
        ptables.cpr must contain the cumulative marginal distributions.
        This method returns [i,j,d], where TF "i" is binding to site, followed by distance
        "d", and then TF "j" binds."""
        #return [0,0,0]
        
        if ptables.cpr[site] == 0.0:
            return None
        totp = ptables.cpr[site]
        randp = random.random() * totp
        
        x1 = site*ptables.dim1
        sump = 0.0                
        for i in ap.params["trlist"]:
            x2 = i*ptables.dim2
            for j in ap.params["trlist+"]:
                x3 = j*ptables.dim3
                for d in  ap.params["rangegd"]:
                    sump += ptables.cpa[ x1 + x2 + x3 + d]
                    if sump > randp:
                        return [i, j, d, randp]
        return [M-1, M, D, randp]
        
            
    #
    # pe_helper is the work required to process one sample from the CDF.
    # See prob_expr.
    #
    def pe_helper(self, genome, ptables, gene, ap):
        """"1. Build a configuration (i.e., the variable named this_config), by sampling 
        cells from ptables.cpa.  Each configuration is a linear collection of transcription
        factors arranged at various sites along the regulatory sequence. """          
        urslen = gene.urs.__len__()
        this_config = {} # key = site, value = the TF bound starting at this site.
        site = 0 # a counter we'll increment as we add stuff to this configuration
        while (site < urslen):
            if ptables.cpr[site] == 0.0:
                """Skip sites where nothing is binding...."""
                site += 1
            else:                
                totp = ptables.cpr[site]
                randp = random.random() * totp
                x1 = site*ptables.dim1
                sump = 0.0       
                reti = 0
                retj = 0
                retd = 0

                sump = 0.0
                for i in ap.params["trlist"]:
                    x2 = i*ptables.dim2
                    for j in ap.params["trlist+"]:
                        x3 = j*ptables.dim3
                        for d in ap.params["rangegd"]:
                            #sump += 10000
                            sump += ptables.cpa[ x1 + x2 + x3 + d]
                            if sump > randp:
                                reti = i
                                retj = j
                                retd = d
                                #print "401:", reti, retj, retd
                                break
                        if sump > randp:
                            break
                    if sump > randp:
                        break

                i = reti
                j = retj
                d = retd

                """Regardless of which sample CDF method we use, we next do this:"""
                this_config[site] = i
                
                """Save the configuration..."""
                #configurations[site].append( [i,j,d] )

                """Advance the site counter..."""
                #if i < ap.params["numtr"]:
                site += genome.genes[i].pwm.P.__len__()
                if site < urslen:
                    site += d
                if j < ap.params["numtr"]:
                    jpwmlen = genome.genes[j].pwm.P.__len__()
                    if (site + jpwmlen) < urslen:
                        this_config[site] = j
                        site += jpwmlen
            
        """"2. Calculate the binding energy of the configuration."""
        sum_lambda_act = 0.0
        sum_lambda_rep = 0.0
        for site in this_config:
            tf = this_config[site]
            if tf < ap.params["numtr"]: # there will be no binding energy for the empty configuration:
                """Get the strength of TF binding at this site..."""
                thispwm = genome.genes[tf].pwm
                this_aff = thispwm.affinity(gene.urs[site:site+thispwm.P.__len__()])
                
                if genome.genes[tf].is_repressor:
                    sum_lambda_rep += this_aff #* tf_expr_levels[tf] #Aug22,2013, Victor: add this multiplier clause tf_expr_levels[tf]
                else:
                    sum_lambda_act += this_aff #* tf_expr_levels[tf] #Aug22,2013, Victor: add this multiplier clause tf_expr_levels[tf]
        
        
        """ 3. Calculate the probability of this configuration."""         
        #print sum_lambda_act-sum_lambda_rep
        try:
            this_pe = (1/( 1+math.exp(-1*ap.params["pe_scalar"]*(sum_lambda_act-sum_lambda_rep) ) ))
        except OverflowError:
            this_pe = MAX_PE
        self.lock.acquire()
        self.pe_sum += this_pe
        self.lock.release()


    def prob_expr(self, genome, ptables, gene, tf_expr_levels, ap, lam=None):        
        """Returns a floating-point value corresponding to the expression level of gene,
        given the ProbTable ptables"""        
        #print "prob_expr, gene", gene.id, gene.urs
        
        pe_sum = 0
        #min_r = min( self.r )
        urslen = gene.urs.__len__()
        
        # configurations is depricated
        #configurations = {}
        #for site in xrange(0, urslen):
        #    configurations[site] = []
        """configurations: key = site, value = array of arrays, [i,j,d] samples"""
        
        # September 2013:
        maxi = ap.params["trlist"][ ap.params["trlist"].__len__()-1 ]
        maxj = maxi + 1
        maxd = ap.params["rangegd"][ ap.params["rangegd"].__len__()-1 ]
        
        for sample in xrange(0, ap.params["iid_samples"]):
            
            #self.pe_helper(genome, ptables, gene, ap)
            
            #t = threading.Thread(target = self.pe_helper, args=(genome, ptables, gene, ap))
            #t.start()
            
            """"1. Build a configuration (i.e., the variable named this_config), by sampling 
            cells from ptables.cpa.  Each configuration is a linear collection of transcription
            factors arranged at various sites along the regulatory sequence. """          
            this_config = {} # key = site, value = the TF bound starting at this site.
            site = 0 # a counter we'll increment as we add stuff to this configuration
            while (site < urslen):
                if ptables.cpr[site] == 0.0:
                    """Skip sites where nothing is binding...."""
                    site += 1
                else:
                    """Sample a configuration..."""     
                    #[i, j, d, randp] = self.sample_cdf(site, ptables, ap)                                  
                    #print "367:", i, j, d
                    
                    #[i,j,d] = [0,0,0]
                    
                    """An alternative to sample_cdf..."""
                    #if ptables.cpr[site] == 0.0:
                    #    return None
                    
                    totp = ptables.cpr[site]
                    randp = random.random() * totp
                    x1 = site*ptables.dim1
                    sump = 0.0       
                    reti = 0
                    retj = 0
                    retd = 0

#
#                    #C-based version:
#
                    #i = _chi2.chi2(site, totp, randp, ptables.cpa)
                    #j = 0
                    #d = 0

#                    #The while-loop version appears to be slower than the for-loop version.                    
#                    q = 0
#                    while q < ptables.dim1:
#                        sump += ptables.cpa[q]
#                        if sump > randp:
#                            reti = int(q / ptables.dim2)
#                            retj = int(q / ptables.dim3)
#                            retd = q - reti*ptables.dim2 - retj*ptables.dim3 
#                            #print "388:", q, reti, retj, retd
#                            q = ptables.dim1
#                        q += 1
                    
                    #
                    # V2:
                    #
                    sump = 0.0
                    for i in ap.params["trlist"]:
                        x2 = i*ptables.dim2
                        for j in ap.params["trlist+"]:
                            x3 = j*ptables.dim3
                            for d in ap.params["rangegd"]:
                                #sump += 10000
                                sump += ptables.cpa[ x1 + x2 + x3 + d]
                                if sump > randp:
                                    reti = i
                                    retj = j
                                    retd = d
                                    #print "401:", reti, retj, retd
                                    break
                            if sump > randp:
                                break
                        if sump > randp:
                            break

                    i = reti
                    j = retj
                    d = retd

                    """Regardless of which sample CDF method we use, we next do this:"""
                    this_config[site] = i
                    
                    # depricated
                    #"""Save the configuration..."""
                    #configurations[site].append( [i,j,d] )

                    """Advance the site counter..."""
                    #if i < ap.params["numtr"]:
                    site += genome.genes[i].pwm.P.__len__()
                    if site < urslen:
                        site += d
                    if j < ap.params["numtr"]:
                        jpwmlen = genome.genes[j].pwm.P.__len__()
                        if (site + jpwmlen) < urslen:
                            this_config[site] = j
                            site += jpwmlen
                
            """"2. Calculate the binding energy of the configuration."""
            sum_lambda_act = 0.0
            sum_lambda_rep = 0.0
            for site in this_config:
                tf = this_config[site]
                if tf < ap.params["numtr"]: # there will be no binding energy for the empty configuration:
                    """Get the strength of TF binding at this site..."""
                    thispwm = genome.genes[tf].pwm
                    this_aff = thispwm.affinity(gene.urs[site:site+thispwm.P.__len__()])
                    
                    if genome.genes[tf].is_repressor:
                        sum_lambda_rep += this_aff #* tf_expr_levels[tf] #Aug22,2013, Victor: add this multiplier clause tf_expr_levels[tf]
                    else:
                        sum_lambda_act += this_aff #* tf_expr_levels[tf] #Aug22,2013, Victor: add this multiplier clause tf_expr_levels[tf]
            
            
            """ 3. Calculate the probability of this configuration."""         
            #print sum_lambda_act-sum_lambda_rep
            try:
                this_pe = (1/( 1+math.exp(-1*ap.params["pe_scalar"]*(sum_lambda_act-sum_lambda_rep) ) ))
            except OverflowError:
                this_pe = MAX_PE
            pe_sum += this_pe
        
        if ap.params["verbosity"] > 3:
            """Print what we've sampled..."""
            self.print_binding_distribution( ptables, genome, gene, ap )
        return (self.pe_sum / ap.params["iid_samples"])
   
    def print_binding_distribution(self, ptables, genome, gene, ap):
        """configs is a hashtable of configurations...
            configs[site] = array of triples [gene i, gene j, distance]"""
        foutpath = ap.params["workspace"] + "/" + ap.params["runid"] + "/" + CONFIG_HISTORY + "/config.gen" + ap.params["generation"].__str__() + ".gid" + genome.id.__str__() + ".txt"
        try:
            fout = open( foutpath , "a")
        except IOError:
            fout = open(ap.params["workspace"] + "/" + ap.params["runid"] + "/LOGS/error.txt", "a")
            #fout.write("\n" + time.localtime() + "\n")
            fout.write("I had a problem opening " + foutpath + "\n")
            fout.close()
            print "\n. I had a problem opening", foutpath
            return
            #exit()

        if self.t_counter == 0 and gene.id == 0:        
            fout.write("===================================================\n")
            fout.write("This file expresses the occupancy distribution for\nindividual " + genome.id.__str__() + " at generation " + self.t_counter.__str__() + ".\n")
            #fout.write("\nKey:\n")
            #fout.write("[+] : activator\n")
            #fout.write("[-] : repressor\n")
            #fout.write("p   : proportion of IID samples\n")
            fout.write("===================================================\n\n")
            #fout.write("e   : energy bound by TF in this configuration\n\n")
        
        #fout.write("#\n# For example, \"site 3 :  1 [-] 0.672 0.120\" indicates that site 3 is bound by repressor 1\n# in 67.2 percent of IID samples with activation energy = 0.120.\n# Activation energies < 0.5 indicate repression, = 0.5 indicate no activation,\n# > 0.5 indicates activation.#\n#\n")
        fout.write(". Time " + self.t_counter.__str__() + " Gene " + gene.id.__str__() + "\tURS: " + gene.urs + "\n")
        
        for site in xrange(0, ptables.cpr.__len__()):
            sump = ptables.cpr[site]
            line = None
            for tf in ap.params["trlist"]:
                #print "438", (site+1).__str__(), ptables.cpt[i, site]
                sumt = ptables.cpt[tf, site]
                if line == None:
                    line = "site " + (site+1).__str__() + " : "
                #pe = genome.genes[tf].pwm.affinity(site, gene.urs)
                if genome.genes[tf].is_repressor:
                    tf_type = "[-]"
                else:
                    tf_type = "[+]"
                if sump == 0:
                    p = 0
                else:
                    p = sumt/sump
                line += "\t" + tf.__str__() + " " + tf_type + " %.3f"%(p)# " p= %.3f"%(p)# + " e= %.3f"%pe
            if line != None:
                fout.write(line + "\n")
        fout.write("\n")
        fout.close()
        
        