from configuration import *

class PWM:
    P = None # index = site, value = hashtable, where key = char, value = bits
    rangesites = []
    
    def __init__(self, copyfrom=None):
        self.P = []
        self.rangesites = []
        if copyfrom != None:
            for i in range(0, copyfrom.P.__len__()):
                self.P.append({})
                for c in ALPHABET:
                    self.P[i][c] = copyfrom.P[i][c]
                self.rangesites.append(i)
        
    def __str__(self):
        str = ""
        for i in range(0, self.P.__len__()):
            str += "(site " + i.__str__() + ")"
            for c in ALPHABET:
                str += " " + c + ":%.3f"%self.P[i][c]                
            str += "\n"
        return str
    
    def randomize(self):
        self.P = []
        self.rangesites = []
        for i in range(0, INIT_PWM_LEN):
            self.P.append( {} )
            sump = 0.0
            shuffled_alphabet = ALPHABET
            random.shuffle(shuffled_alphabet)
            for c in shuffled_alphabet:
                thisp = random.uniform(0.0, 1.0 - sump)
                self.P[i][c] = thisp
                sump += thisp
            self.rangesites.append(i)
    
    def mutate(self, ap):
        rand_site = random.randint(0, self.P.__len__()-1)
        rand_state = random.randint(0,3)
        
        d = random.random() * ap.params["pwmmu"]
        new_p = (self.P[rand_site][ ALPHABET[rand_state] ] + d)%1.0
        print "pwm.46 mutating PWM site", rand_site, ALPHABET[rand_state], "%.3f"%self.P[rand_site][ ALPHABET[rand_state]], "%.3f"%new_p

        self.P[rand_site][ ALPHABET[rand_state] ] = new_p

        # . . . and then normalize the PWM. . . 
        sum_states = 0.0
        for c in ALPHABET:
            #print self.P[rand_site][c]
            sum_states += self.P[rand_site][c]
        for c in ALPHABET:
            self.P[rand_site][c] = self.P[rand_site][c] / sum_states
        #print self.P
    
    def read_from_file(self, path, id):
        self.P = []
        self.rangesites = []
        fin = open(path, "r")
        i = -1
        found_our_id = False
        for l in fin.readlines():
            if l.startswith("#"):
                continue
            if l.startswith("pwm"):
                tokens = l.split()
                if int(tokens[1]) == id:
                    found_our_id = True
                if int(tokens[1]) > id:
                    found_our_id = False
                continue
            if found_our_id == True:
                if l.__len__() > 2 and False == l.__contains__("#"):
                    i += 1
                    tokens = l.split()
                    cc = 0
                    self.P.append( {} )
                    self.rangesites.append(i)
                    #print tokens
                    for c in ALPHABET:
                        self.P[i][c] = float(tokens[cc])
                        cc += 1
        fin.close()
    
    def make_flat(self, ap):
        """flattens the P matrix to be non-specific, with length len"""
        flat_p = 1.0 / ALPHABET.__len__()
        self.P = []
        for i in range(0, ap.params["init_pwm_len"]):
            self.P.append( {} )
            for c in ALPHABET:
                self.P[i][c] = flat_p
    
    def specificity(self, pos, urs):
        """What is the probability that the PWM in gene tf will bind the upstream 
        regulatory region (urs) sequence with the right edge of the PWM at site pos?"""        
        if (urs.__len__() - pos < self.P.__len__() ):
            return 0.0
        if ( pos + 1 > urs.__len__() ):
            return 0.0
        a = ALPHABET.__len__()
        part = urs[pos:pos+self.P.__len__()] # the partial URS
        #print "pos=", pos, "part = ", part
        res1 = 0.0        
        
        """In Kevin Bullaughey's original Rescape code he assumed 
            the affinity of a TF is the sum of the affinity to both
            strand at this location.
            But in my code, I'm only dealing with single stranded DNA."""
        for k in self.rangesites:
            #if self.P[k][ part[k] ] == 0.0:
            #    continue
            #else:
            res1 += self.P[k][ part[k] ]
            #print "res1=", res1
        
        """In contrast, this is what Kevin's code does:"""
        #res1 *= a**( self.P.__len__() )
        
        return res1
    
    
    