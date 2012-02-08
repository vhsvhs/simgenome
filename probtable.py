from configuration import *

#class Parray:
#    """ProbTable Array"""
#    dims = []
#    mem = []
#    
#    def __init__(self, len):
#        print "len = ", len
#        self.mem = [0]*len
#
#    def pos(self, nops):
#        """
#        Returns a cell index in mem.
#        nops are directions through the dimensions.
#        Go nops[0] units along dims[0], then nops[1] units along dims[1], etc.
#        """
#        
#        cell = 0
#        dimprod = 1
#        for i in nops:
#           cell += dimprod * i
#           dimprod *= self.dims[i] 
#        return cell
#
#    def posi(self, i):
#        return self.mem[i]

class ProbTable:
    """The ProbTable class wraps four arrays, cpa, cpt, cpr, and cpm."""
    cpa = None
    cpt = None
    cpr = None
    cpm = None
    
    def __init__(self, M, D, L):
        self.cpa = zeros( (M,(M+1),D,L), dtype=float)
        #cpalen = M * (M+1) * D * L + 2
        #self.cpa = Parray(cpalen)
        #self.cpa.dims = [M, (M+1), D, L]
        
        #cptlen = M * L
        #self.cpt = Parray(cptlen)
        self.cpt = zeros((M, L), dtype=float)

        #cprlen = L
        #self.cpr = Parray(cprlen)
        self.cpr = zeros((L), dtype=float)

        #cpmlen = (D*(M+1)+1) * M * L
        #self.cpm = Parray(cpmlen)
        self.cpm = zeros( ((D*(M+1)+1),M,L), dtype=float)   
  