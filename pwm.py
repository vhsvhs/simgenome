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
    
#    def collapse(self):
#        ret = {}
#        for site in range(0, self.P.__len__()):
#            ret[site] = {}
#            for state in ALPHABET:
#                ret[site] = self.P[site][state]
            
    
    
    def generate_random_site(self):
        p = {}
        sump = 0.0
        shuffled_alphabet = ALPHABET
        random.shuffle(shuffled_alphabet)
        for c in shuffled_alphabet:
            thisp = random.uniform(0.0, 1.0 - sump)
            p[c] = thisp
            sump += thisp
        return p
    
    def randomize(self):
        self.P = []
        self.rangesites = []
        for i in range(0, INIT_PWM_LEN):
            self.P[i] = self.generate_random_site()
            self.rangesites.append(i)
    
    def mutate(self, ap):
        rand_site = random.randint(0, self.P.__len__()-1)
        rand_state = random.randint(0,3)
        
        d = random.random() * ap.params["pwmdeltamax"]
        new_p = (self.P[rand_site][ ALPHABET[rand_state] ] + d)%1.0
        print "\t> mutating PWM site", rand_site, ". ML state =", ALPHABET[rand_state], ", old p = %.3f"%self.P[rand_site][ ALPHABET[rand_state]], ", new p = %.3f"%new_p

        self.P[rand_site][ ALPHABET[rand_state] ] = new_p

        """ Normalize the PWM"""
        sum_states = 0.0
        for c in ALPHABET:
            sum_states += self.P[rand_site][c]
        for c in ALPHABET:
            self.P[rand_site][c] = self.P[rand_site][c] / sum_states

    def mutate_len(self, ap):
        d = random.random() # mutate or not?
        if d < ap.params["pwmlenmu"]:
            size = int(random.random() * ap.params["pwmmulenmax"])
            if size < 1:
                size = 1
            e = random.random() # insert or deletion?
            if e > 0.5: #insertion
                if ap.params["verbosity"] >= 80:
                    print "insertion"
                newP = []
                rand_site = random.randint(0, self.P.__len__()-1)
                #print "82:", self.P, rand_site
                for i in range(0, rand_site):
                    #print i
                    newP.append( copy.deepcopy( self.P[i] ) )
                for i in range(0, size):
                    #print i
                    newP.append( self.generate_random_site() )
                for i in range(rand_site+size, self.P.__len__()):
                    newP.append( copy.deepcopy(self.P[i-size]) )

                #print dir(copy) 
                #self.P = copy.deepcopy( newP )
                self.P = newP

                self.rangesites = []
                for i in range(0, newP.__len__()):
                    self.rangesites.append(i)
                #print "pwm 81:", self.rangesites
            elif self.P.__len__() > 2: # deletion
                if ap.params["verbosity"] >= 80:
                    print "deletion"
                for count in range(0, size):
                    rand_site = random.randint(0, self.P.__len__()-1)
                    for i in range(rand_site+1, self.P.__len__()):
                        for c in ALPHABET:
                            self.P[i-1][c] = self.P[i][c]
                    self.P.pop(self.P.__len__()-1)
                    self.rangesites = self.rangesites[0:self.rangesites.__len__()-1]
                    #print "pwm 88:", self.rangesites
            return True
        return False
        
    def read_from_file(self, path, id):
        """Reads only those lines that correspond to the PWM for TF #id."""
        self.P = []
        self.rangesites = []
        fin = open(path, "r")
        i = -1
        found_our_id = False
        reg_mode = 0
        for l in fin.readlines():
            if l.__len__() < 2:
                continue
            if l.startswith("#"):
                continue
            if l.__contains__("pwm"):
                tokens = l.split()
                if tokens[1] == id:
                    found_our_id = True
                    reg_mode = int(tokens[2]) # 0 = activator, 1 = repressor
                else:
                    if found_our_id == True:
                        return (found_our_id, reg_mode)
                    found_our_id = False
                continue
            if found_our_id == True:
                if l.__len__() > 2 and False == l.__contains__("#"):
                    i += 1
                    tokens = l.split()
                    cc = 0
                    self.P.append( {} )
                    self.rangesites.append(i)
                    for c in ALPHABET:
                        self.P[i][c] = float(tokens[cc])
                        cc += 1
        fin.close()
        return (found_our_id, reg_mode)
    
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
    
    
    