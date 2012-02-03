from configuration import *

class PWM:
    P = None # index = site, value = hashtable, where key = char, value = bits
    
    def __init__(self):
        self.P = []
    
    def randomize(self, len):
        """generate a random P matrix with length len"""
        self.P = []
        for i in range(0, len=INIT_PWM_LEN):
            self.P.append( {} )
            sump = 0.0
            for c in ALPHABET:
                thisp = random.uniform(0.0, 1.0 - sump)
                self.P[i][c] = thisp
                sump += thisp
    
    def make_flat(self, len=INIT_PWM_LEN):
        """flattens the P matrix to be non-specific, with length len"""
        flat_p = 1.0 / ALPHABET.__len__()
        self.P = []
        for i in range(0, len):
            self.P.append( {} )
            for c in ALPHABET:
                self.P[i][c] = flat_p 
            