"""
Notes on the class PWM.

PWM's are a data structure to express and operate on position-specific affinity matricies (PSAMs).
PSAMs can be calculated from data gathered with technology such as MITOMI.
See Foat et al. 2005 Bioinformatics for more information about PSAMs.
See Fordyce et al. 2010 Nature Biotech. for more information about MITOMI.

PSAMs express the relative affinity difference between nucleotides at DNA-binding sites.
The values for PSAMs are unbounded. Zero for nucleotide X means that the existence of X
at the binding site in question will have no effect on the affinity of the protein for the site.
A positive value indicates that affinity will increase if state X is present, while a negative
value indicates that affinity will decrease if state X is present.
In cases where affinity will increase or decrease, the magnitude of change is equal to:

   = e^(value)
   
The relative affinity of a sequence is therefore the product of magnitudes for each site.   
"""
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
        """Generates a hashtable with one key for each letter in the alphabet, and with random values for those
        keys.  Essentially, this generates one random site in a PSAM."""
        p = {}
        sump = 0.0
        shuffled_alphabet = ALPHABET
        random.shuffle(shuffled_alphabet)
        for c in shuffled_alphabet:
            p[c] = random.uniform(-1.0, 1.0)
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
        sign_rand = random.uniform(0.0,1.0)
        if sign_rand > 0.5:
            new_p *= -1
        print "\t\t> Site", rand_site, ": ML state =", ALPHABET[rand_state], ", old affinity = %.3f"%self.P[rand_site][ ALPHABET[rand_state]], ", new affinity = %.3f"%new_p

        self.P[rand_site][ ALPHABET[rand_state] ] = new_p

        #
        # Aug 2013: depricated:
        # normalization is not required for PSAMs
        #
        #""" Normalize the PWM"""
        #sum_states = 0.0
        #for c in ALPHABET:
        #    sum_states += self.P[rand_site][c]
        #for c in ALPHABET:
        #   self.P[rand_site][c] = self.P[rand_site][c] / sum_states

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
                    #print "111:", self.P[i-size]
                    #print "112:", copy.deepcopy(self.P[i-size])
                    #print "113:", newP
                    x = {}
                    for key in self.P[i]:
                        x[key] = self.P[i][key]
                    newP.append(x)
                for i in range(0, size):
                    #print i
                    newP.append( self.generate_random_site() )
                for i in range(rand_site+size, self.P.__len__()):
                    #print "116:", self.P[i-size]
                    #print "117:", copy.deepcopy(self.P[i-size])
                    #print "118:", newP
                    x = {}
                    for key in self.P[i-size]:
                        x[key] = self.P[i-size][key]
                    newP.append(x)
                    #newP.append( copy.deepcopy(self.P[i-size]) )

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
    
#    #
#    # Depricated.  Now use cli.py 
#    #
#    def read_from_file(self, path, id):
#        """Reads only those lines that correspond to the PWM for TF #id."""
#        self.P = []
#        self.rangesites = []
#        fin = open(path, "r")
#        i = -1
#        found_our_id = False
#        reg_mode = 0
#        for l in fin.readlines():
#            if l.__len__() < 2:
#                continue
#            if l.startswith("#"):
#                continue
#            if l.__contains__("pwm"):
#                tokens = l.split()
#                if tokens[1] == id:
#                    found_our_id = True
#                    reg_mode = tokens[2] # 0 = activator, 1 = repressor
#                    if reg_mode == "+" or reg_mode == "1":
#                        reg_mode = 0
#                    elif reg_mode == "-" or reg_mode == "0":
#                        reg_mode = 1
#                else:
#                    if found_our_id == True:
#                        return (found_our_id, reg_mode)
#                    found_our_id = False
#                continue
#            if found_our_id == True:
#                if l.__len__() > 2 and False == l.__contains__("#"):
#                    i += 1
#                    tokens = l.split()
#                    cc = 0
#                    self.P.append( {} )
#                    self.rangesites.append(i)
#                    for c in ALPHABET:
#                        self.P[i][c] = float(tokens[cc])
#                        cc += 1
#        fin.close()
#        return (found_our_id, reg_mode)
    
    def make_flat(self, ap):
        """flattens the P matrix to be non-specific, with length len"""
        flat_p = 1.0 / ALPHABET.__len__()
        self.P = []
        for i in range(0, ap.params["init_pwm_len"]):
            self.P.append( {} )
            for c in ALPHABET:
                self.P[i][c] = flat_p
    
    def affinity(self, pos, urs):
        """Returns the total affinity of this transcription factor for the
        upstream regulatory sequence (URS) bound at a given position (pos).
        Total affinity is the product of affinities at each site.
        """
        
        # If the PWM can't fit on the sequence:
        if (urs.__len__() - pos < self.P.__len__() ):
            return 0.0
        # If the position is beyond the edge of the sequence:
        if ( pos + 1 > urs.__len__() ):
            return 0.0
        
        a = ALPHABET.__len__()
        part = urs[pos:pos+self.P.__len__()] # part is the substring of the sequence

        retaff = 1.0 # the affinity, to be returned
        for k in self.rangesites:
            retaff *= e**( self.P[k][ part[k] ] )
        return retaff
    
    #
    # depricated.  Now use the function "affinity"
    #
#    def specificity(self, pos, urs):    
#        """What is the probability that the PWM in gene tf will bind the upstream 
#        regulatory region (urs) sequence with the right edge of the PWM at site pos?"""        
#        if (urs.__len__() - pos < self.P.__len__() ):
#            return 0.0
#        if ( pos + 1 > urs.__len__() ):
#            return 0.0
#        a = ALPHABET.__len__()
#        part = urs[pos:pos+self.P.__len__()] # the partial URS
#        #print "pos=", pos, "part = ", part
#        res1 = 0.0        
#        
#        """In Kevin Bullaughey's original Rescape code he assumed 
#            the affinity of a TF is the sum of the affinity to both
#            strand at this location.
#            But in my code, I'm only dealing with single stranded DNA."""
#        for k in self.rangesites:
#            #if self.P[k][ part[k] ] == 0.0:
#            #    continue
#            #else:
#            res1 += self.P[k][ part[k] ]
#            #print "res1=", res1
#        
#        """In contrast, this is what Kevin's code does:"""
#        #res1 *= a**( self.P.__len__() )
#        
#        return res1
    
    
    