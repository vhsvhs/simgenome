from configuration import *
from genome import *

class Fitness_Goal:
    """ideal_expression[gene x][time slice i] = 
        = desired expression level of x at i.
        """
    ideal_expression = {}
    
    def __init__(self, basal_gene_id, levels):
        pass

class Landscape:
    
    """An array of Fitness_Goals"""
    fitness_goals = []
    
    def __init__(self):
        pass
    
    def create_random_fitness_goal(self, genome):
        """Creates a random fitness goal, using genome as the seed."""    
        
        """First, select a basal TR to activate the goal."""
        randid = random.randint(0, N_TR)
        
        """Next, select a random set of desired expression pulses."""
        levels = {}
        for i in range(0, GOAL_COMPLEXITY):
            
        
        goal = Fitness_Goal( randid )
        self.fitness_goals.append( goal )
    
    def get_fitness(self, genome):
        """Calculates the fitness of the given genome.  Returns a floating-point value."""
        pass
    
    