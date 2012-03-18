from configuration import *
from pwm import *

class Gene:
    urs = ""        # upstream regulatory sequence
    has_dbd = False # False = gene is reporter, True = gene is regulatory gene
    pwm = None
    id = None
    is_repressor = False # True means it's an activator
        
    def __init__(self, id, urs = None, has_dbd = False, repressor = False, pwm = None):
        """1. id"""
        self.id = id
        
        """2. urs"""
        if urs == None:
            """Create a random URS"""
            for i in range(0, URS_LEN):
                """Sample a random character from the alphabet, where all chars are equally probable."""
                self.urs += random.sample(ALPHABET, 1)[0]
        else:
            self.urs = urs
        
        """3. has_dbd"""
        self.has_dbd = has_dbd
                
        """4. is_repressor"""
        self.is_repressor = repressor

        """5. pwm"""
        if (self.has_dbd and pwm != None):
            self.pwm = pwm
        elif self.has_dbd:
            self.pwm = PWM()
            #self.pwm.make_flat() # initialize all PWMs to be non-specific
            self.pwm.randomize()

    def collapse(self):
        return [self.id, self.urs, self.has_dbd, self.is_repressor, self.pwm]

