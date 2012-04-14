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
        
    def __str__(self):        
        ret = ""
        ret += "CPA\nsite\tTF\tTF\td\tVALUE\n"
        dims = self.cpa.shape
        for x in range(0, dims[3]):
            for i in range(0, dims[0]):
                for j in range(0, dims[1]):
                    for d in range(0, dims[2]):
                        ret += x.__str__() + "\t"
                        ret += i.__str__() + "\t"
                        ret += j.__str__() + "\t"
                        ret += d.__str__() + "\t"
                        ret += self.cpa[i,j,d,x].__str__() + "\n"
        return ret

        
  