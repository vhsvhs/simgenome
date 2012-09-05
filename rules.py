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
    """Each Fitness_RUle object defines a particular expression gate that
    must be achieved for maximal fitness.  See the class Rule_Collection, in which
    multiple Fitness_Rule instances are combined to define an overall fitness landscape."""
    timepoint = None
    expression_level = None
    reporter_id = None
    rule_type = None # See RULE_TYPES

    def __init__(self, t, e, id, r):
        self.timepoint = t
        self.expression_level = e
        self.reporter_id = id
        self.rule_type = r
    
    def __str__(self):
        l = "  + At time "
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

class Rule_Collection:
    """Each instance of the class Rule_Collection defines a unique set of temporal expression patterns
    for a subset of reporter genes in each genome.  See the parent class, Landscape, for how to
    calculate the fitness of a given Rule_Collection.  This class can be initialized with an a priori
    desired expression pattern, or the pattern can be randomly initialized."""
    collection_id = None
    rules = []
    
    def __init__(self, collection_id):
        self.collection_id = collection_id
        self.rules = []
        
    def collapse(self):
        data = []
        for r in self.rules:
            data.append( r.collapse() )
        return data
    
    def uncollapse(self, data):
        for rule in data:
            this_rule = Fitness_Rule(rule[0], rule[1], rule[2], rule[3])
            self.rules.append( this_rule )
    
    def init_random(self):
        """This method generates a random set of Fitness_Rule objects."""
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
        ret += "--> Fitness Goals:\n"
        for r in self.rules:
            ret += r.__str__()
        return ret
