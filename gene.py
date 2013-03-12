from configuration import *
from pwm import *

class Gene:
    urs = ""        # "urs" stands for upstream regulatory sequence
    has_dbd = False # False = gene is reporter, True = gene is regulatory gene
    pwm = None      # "pwm" stands for position specific weight matrix.
    id = None
    name = ""
    is_repressor = False
    tfcoop = None # array of relative co-factor affinities, without consideration of distance. [-1,0) = anti-co-factor activity. 0 = no activity, (0,+infinity] = positive co-factor activity
    gamma = None # i.e., cofactor affinity at various distances.  This is calculated from self.tfcoop.  It's a 2-d Numpy array. gamma[tf][distance].  Values: 1.0 = no effect on binding, <1.0 = decreases binding, >1.0 = strengthens binding.
        
    def __init__(self, id, urs_len, urs = None, has_dbd = False, repressor = False, pwm = None, tfcoop = None, gamma = None, ap = None, name = None):
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
            #print "gene 27", pwm
        elif self.has_dbd == True:
            self.pwm = PWM()
            self.pwm.randomize()
        
        """6. Co-factor affinity"""
        if self.has_dbd and tfcoop == None:
            self.set_tfcoop(ap)
        elif self.has_dbd and tfcoop != None:
            self.tfcoop = tfcoop
            self.gamma = gamma
            
        if name != None:
            self.name = name

    def collapse(self):
        #print "gene.py47 collapsing into", [self.id, self.urs.__len__(), self.urs, self.has_dbd, self.is_repressor, self.pwm, self.gamma]
        return [self.id, self.urs.__len__(), self.urs, self.has_dbd, self.is_repressor, self.pwm, self.tfcoop, self.gamma, self.name]

    def coopfunc(self, g, d):
        """Calculates the degree of binding cooperativity between two TFs binding distance d apart.
        g is positive for synergistic interactions, and negative for antagonistic interactions.
        g ranges from -1 to +inifinity. g of zero means no effect on binding."""
        return 1 + g * math.exp( (-1)*(d**2)/V_RATE_OF_COOP_DECAY );

    def set_tfcoop(self, ap):
        if ap == None:
            if self.tfcoop == None:
                print "\n. Error gene.py line 56.  You need to specify ap= or gamma= when you call the function set_gamma."
                exit(1)
        
        # tfcoop is an intermediary matrix that holds random draws from the gamma distribution.
        # tfcoop becomes transformed into self.gamma
        if ap.params["coopinit"] == "random":
            self.tfcoop = numpy.random.gamma(2.0,5.0, (ap.params["numtr"]) ) - 4.0
        else:
            self.tfcoop = zeros( (ap.params["numtr"]), dtype=float)
        self.set_gamma(ap)


    def set_gamma(self, ap):
        self.gamma = zeros( (ap.params["numtr"], ap.params["maxgd"]), dtype=float)
        for i in ap.params["trlist"]:
            for d in ap.params["rangegd"]:
                self.gamma[i,d] = self.coopfunc( self.tfcoop[i], d)
        #print "debug gene.py 75", " rank=", comm.Get_rank()," gamma=", self.gamma
        
        
    def mutate_urs(self):
        """Pick a random site"""
        rand_site = random.randint(0, self.urs.__len__()-1)
        """Mutate!"""
        curr_state = self.urs[rand_site]
        ALPHABET.remove(curr_state)
        new_state = random.choice( ALPHABET )
        ALPHABET.append(curr_state)
        new_urs = ""
        for j in range(0, self.urs.__len__()):
            if j == rand_site:
                new_urs += new_state.__str__()
            else:
                new_urs += self.urs[j]
        self.urs = new_urs
    
    def mutate_urs_len(self):
        """Pick a random site at which we will insert or delete"""
        rand_site = random.randint(0, self.urs.__len__()-1)        
        mulen = random.randint(1, URS_LEN_INDEL_MAX)
        action = random.randint(0, 1)
        if self.urs.__len__() <= MIN_URS_LEN:
            action = 1
        if action == 0:
            # delete
            newurs = ""
            for j in range(0, rand_site): 
                newurs += self.urs[j]
            if rand_site+mulen < self.urs.__len__():
                for j in range(rand_site+mulen, self.urs.__len__()):
                    newurs += self.urs[j]
            self.urs = newurs
        elif action == 1:
            #insert
            newurs = ""
            for j in range(0, rand_site): 
                newurs += self.urs[j]
            for j in range(0,mulen):
                newurs += ALPHABET[random.randint(0, ALPHABET.__len__()-1 )]
            for j in range(rand_site, self.urs.__len__()):
                newurs += self.urs[j]
            self.urs = newurs            
    
    def mutate_gamma(self, ap):
        """Modify this TF's preference for one rand_tr_id gene."""
        rand_tr_id = random.randint(0, ap.params["numtr"]-1)
        delta = random.randint(-1,10)
        print "gene.py 100:", self.tfcoop
        #print "\t gene.py100 ", self.tfcoop[rand_tr_id], " + ", delta
        self.tfcoop[rand_tr_id] += delta
        if self.tfcoop[rand_tr_id] < -1:
            self.tfcoop[rand_tr_id] = -1.0
        self.set_gamma(ap)