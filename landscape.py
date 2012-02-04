from configuration import *
from genome import *

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
        pass
    
    