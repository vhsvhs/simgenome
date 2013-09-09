from configuration import *

class ProbTable:
    """The ProbTable class wraps four arrays, cpa, cpt, cpr, and cpm.
        cpt - full cumulative probability table, all cells.
        cpr - marginal cumulative probability row"""
    
    cpa = None
    cpt = None
    cpr = None
    cpm = None
    
    dim1 = None
    dim2 = None
    dim3 = None
    
    def __init__(self, M, D, L):
        # M is num TR
        # D is max gamma distance
        # L is URS length
        #self.cpa = zeros( (M,(M+1),D,L), dtype=float)
        #self.cpt = zeros( (M, L), dtype=float)
        #self.cpr = zeros( (L), dtype=float)

        #self.cpa = zeros((L, M,(M+1),D), dtype=float)
        self.cpa = [0] * (L*M*(M+1)*D) # a long array
        self.dim1 = M*(M+1)*D
        self.dim2 = (M+1)*D
        self.dim3 = D
        
        self.cpt = zeros( (M, L), dtype=float)
        self.cpr = zeros( (L), dtype=float)    
            
    def __str__(self):        
        ret = ""
#        ret += "CPA\nsite\tTF\tTF\td\tVALUE\n"
#        dims = self.cpa.shape
#        for x in range(0, dims[0]):
#            for i in range(0, dims[1]):
#                for j in range(0, dims[2]):
#                    for d in range(0, dims[3]):
#                        ret += x.__str__() + "\t"
#                        ret += i.__str__() + "\t"
#                        ret += j.__str__() + "\t"
#                        ret += d.__str__() + "\t"
#                        ret += self.cpa[x,i,j,d].__str__() + "\n"

        #
        # to-do: print cpt and cpr
        #
        ret += "\n"
        ret += "CPT\nsite\tTR\tvalue\n"
        dims = self.cpt.shape
        for x in range(0, dims[1]):
            for i in range(0, dims[0]):
                ret += x.__str__() + "\t"
                ret += i.__str__() + "\t"
                ret += self.cpt[i,x].__str__() + "\n"
        
        ret += "\n"
        ret += "CPR\nsite\n"
        dims = self.cpr.shape
        for i in range(0, dims[0]):
            ret += i.__str__() + "\t"
            ret += self.cpr[i].__str__() + "\n"        
        
        return ret