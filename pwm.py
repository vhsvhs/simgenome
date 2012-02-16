from configuration import *

class PWM:
    P = None # index = site, value = hashtable, where key = char, value = bits
    
    def __init__(self):
        self.P = []
    
    def randomize(self):
        """generate a random P matrix with length len"""
        self.P = []
        for i in range(0, INIT_PWM_LEN):
            self.P.append( {} )
            sump = 0.0
            shuffled_alphabet = ALPHABET
            random.shuffle(shuffled_alphabet)
            for c in shuffled_alphabet:
                thisp = random.uniform(0.0, 1.0 - sump)
                self.P[i][c] = thisp
                sump += thisp
    
    def make_flat(self):
        """flattens the P matrix to be non-specific, with length len"""
        flat_p = 1.0 / ALPHABET.__len__()
        self.P = []
        for i in range(0, INIT_PWM_LEN):
            self.P.append( {} )
            for c in ALPHABET:
                self.P[i][c] = flat_p
    
    def prob_binding(self, pos, urs):
        """What is the probability that the PWM in gene tf will bind the upstream 
        regulatory region (urs) sequence with the right edge of the PWM at site pos?"""
        a = ALPHABET.__len__()
        if (pos - self.P.__len__() + 1 < 0) or ( pos >= urs.__len__() ):
            return 0.0
        part = urs[pos:pos+self.P.__len__()] # the partial URS
        res1 = 1
        # In Kevin Bullaughey's original Rescape code he assumed 
        # the affinity of a TF is the sum of the affinity to both
        # strand at this location.
        # But in my code, I'm only dealing with single stranded DNA.
        for k in range(0, self.P.__len__()):
            res1 *= self.P[k][ part[k] ]
        res1 *= a**( self.P.__len__() )
        return res1