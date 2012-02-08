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
        """What is the probability that this PWM will bind the upstream 
        regulatory region (urs) sequence with the right edge of the PWM at site pos?"""
        
        #
        # to-do
        #
        
        return 1.0