from configuration import *
from pwm import *

class Gene:
    urs = ""        # "urs" stands for upstream regulatory sequence
    has_dbd = False # False = gene is reporter, True = gene is regulatory gene
    pwm = None      # "pwm" stands for position specific weight matrix.
    id = None
    is_repressor = False
    gamma = None # i.e., cofactor affinity.  It's a 2-d Numpy array. gamma[tf][distance].  Values: 1.0 = no effect on binding, <1.0 = decreases binding, >1.0 = strengthens binding.
        
    def __init__(self, id, urs_len, urs = None, has_dbd = False, repressor = False, pwm = None, gamma = None, ap = None):
        """gamma and ap are optional, but you need at least one of them."""
        
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
        
        """6. Co-factor affinity"""
        if self.has_dbd and gamma == None:
            self.set_gamma(ap)
        elif self.has_dbd and gamma != None:
            self.gamma = gamma

    def collapse(self):
        #print "gene.py47 collapsing into", [self.id, self.urs.__len__(), self.urs, self.has_dbd, self.is_repressor, self.pwm, self.gamma]
        return [self.id, self.urs.__len__(), self.urs, self.has_dbd, self.is_repressor, self.pwm, self.gamma]

    def coopfunc(self, g, d):
        """Calculates the degree of binding cooperativity between two TFs binding distance d apart.
        g is positive for synergistic interactions, and negative for antagonistic interactions.
        g ranges from -1 to +inifinity. g of zero means no effect on binding."""
        return 1 + g * math.exp( (-1)*(d**2)/V_RATE_OF_COOP_DECAY );

    def set_gamma(self, ap):
        if ap == None:
            if self.gamma == None:
                print "\n. Error gene.py line 56.  You need to specify ap= or gamma= when you call the function set_gamma."
            return
        
        # tfcoop is an intermediary matrix that holds random draws from the gamma distribution.
        # tfcoop becomes transformed into self.gamma
        tfcoop = None
        if ap.params["coopinit"] == "random":
            tfcoop = numpy.random.gamma(2.0,10.0, (ap.params["numtr"]) ) - 4.0
        else:
            tfcoop = zeros( (ap.params["numtr"]), dtype=float)

        self.gamma = zeros( (ap.params["numtr"], ap.params["maxgd"]), dtype=float)
        for i in ap.params["rangetrs"]:
            for d in ap.params["rangegd"]:
                self.gamma[i,d] = self.coopfunc( tfcoop[i], d)
