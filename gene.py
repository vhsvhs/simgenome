from configuration import *
from pwm import *

class Gene:
    urs = ""        # upstream regulatory sequence
    has_dbd = False # False = gene is reporter, True = gene is regulatory gene
    pwm = None
    id = None
    
    def __init__(self, id, urs = None, has_dbd = False):
        self.id = id
        """Initialize the gene.  urs can be a specified sequence, or random"""
        if urs == None:
            """Create a random URS"""
            for i in range(0, URS_LEN):
                """Sample a random character from the alphabet, where all chars are equally probable."""
                self.urs += random.sample(ALPHABET, 1)[0]
        self.has_dbd = has_dbd
        if (self.has_dbd):
            self.pwm = PWM()
            #self.pwm.make_flat() # initialize all PWMs to be non-specific
            self.pwm.randomize()

