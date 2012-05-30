from configuration import *
from pwm import *

class Gene:
    urs = ""        # "urs" stands for upstream regulatory sequence
    has_dbd = False # False = gene is reporter, True = gene is regulatory gene
    pwm = None      # "pwm" stands for position specific weight matrix.
    id = None
    is_repressor = False
        
    def __init__(self, id, urs_len, urs = None, has_dbd = False, repressor = False, pwm = None):
        """1. id"""
        self.id = id
        
        """2. urs"""
        if urs == None:
            """Create a random URS"""
            for i in range(0, urs_len):
                """Sample a random character from the alphabet, where all chars are equally probable."""
                self.urs += random.sample(ALPHABET, 1)[0]
        else:
            self.urs = urs
        
        """3. has_dbd"""
        self.has_dbd = has_dbd
                
        """4. is_repressor"""
        self.is_repressor = repressor

        """5. pwm"""
        if (self.has_dbd == True and pwm != None):
            self.pwm = pwm
        elif self.has_dbd == True:
            self.pwm = PWM()
            self.pwm.randomize()

    def collapse(self):
        return [self.id, self.urs.__len__(), self.urs, self.has_dbd, self.is_repressor, self.pwm]

