from configuration import *
from pwm import *

class Gene:
    urs = ""        # upstream regulatory sequence
    has_dbd = False # False = gene is reporter, True = gene is regulatory gene
    pwm = None
    id = None
    is_repressor = False # True means it's an activator
        
    def __init__(self, id, urs = None, has_dbd = False, repressor = False, pwm = None):
        self.id = id
        """Initialize the gene.  urs can be a specified sequence, or random"""
        if urs == None:
            """Create a random URS"""
            for i in range(0, URS_LEN):
                """Sample a random character from the alphabet, where all chars are equally probable."""
                self.urs += random.sample(ALPHABET, 1)[0]
        self.has_dbd = has_dbd
        if (self.has_dbd and pwm != None):
            self.pwm = pwm
        elif self.has_dbd:
            self.pwm = PWM()
            #self.pwm.make_flat() # initialize all PWMs to be non-specific
            self.pwm.randomize()
        self.is_repressor = repressor
        
    def collapse(self):
        return [self.id, self.urs, self.has_dbd, self.is_repressor, self.pwm]

