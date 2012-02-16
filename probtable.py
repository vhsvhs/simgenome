from configuration import *

class ProbTable:
    """The ProbTable class wraps four arrays, cpa, cpt, cpr, and cpm.
        cpt - full cumulative probability table, all cells.
        cpr - marginal cumulative probability row"""
    
    cpa = None
    cpt = None
    cpr = None
    cpm = None
    
    def __init__(self, M, D, L):
        self.cpa = zeros( (M,(M+1),D,L), dtype=float)
        self.cpt = zeros( (M, L), dtype=float)
        self.cpr = zeros( (L), dtype=float)

    
        
  