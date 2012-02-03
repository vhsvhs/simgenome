from configuration import *

class Gene:
    self.urs = ""        # upstream regulatory sequence
    self.has_dbd = False # False = gene is reporter, True = gene is regulatory gene
    self.pwm = None
    self.id = None
    
    def __init__(self, id, urs = None, has_dbd = False):
        """Initialize the gene.  urs can be a specified sequence, or random"""
        if urs == None:
            """Create a random URS"""
            for i in range(0, URS_LEN):
                self.urs.append( random.sample(ALPHABET, 1) )
        self.has_dbd = has_dbd
        if (self.has_dbd):
            self.pwm = PWM()
            self.pwm.make_flat() # initialize all PWMs to be non-specific
