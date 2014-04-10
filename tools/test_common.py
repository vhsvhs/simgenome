#
# Common functions used by test modules.
#
import math, os, random
from diverging_colors import *
from argparser import *
ap = ArgParser(sys.argv)

def get_rand_nt(len):
    nts = ["A", "C", "G", "T"]
    ret = ""
    for i in range(0, len):
        ret += random.sample(nts, 1)[0]
    return ret

def get_repeat(r, N):
    """Returns N copies of r, concat-ed together.""" 
    ret = ""
    for i in range(0, N):
        ret += r
    return ret

color = {}
lwd = {}

def lwd_for_gene(x):
    if x not in lwd:
        if x == 0:
            this_lwd = "1"
        else:
            this_lwd = "1"
        lwd[x] = this_lwd
    return lwd[x] 

# set is an array of floats
def mean(set):
    if set.__len__() == 0:
        return None
    sum = 0.0
    for x in set:
        sum += x
    return sum / float( set.__len__() )

# standard deviation
def sd(set):
    avg = mean(set)
    if avg == None:
        return None
    sumofsquares = 0.0
    for x in set:
        sumofsquares += (x - avg)**2
    return math.sqrt( float( sumofsquares ) / float( set.__len__() ) )

# calculates variance
def var(set):
    avg = mean(set)
    if avg == None:
        return None
    sumofsquares = 0.0
    for x in set:
        sumofsquares += (x - avg)**2
    return math.sqrt( float( sumofsquares ) / float( set.__len__() - 1 ) ) 

def stderr(set):
    return (sd(set) / math.sqrt( set.__len__() ) )

def get_maxfit_id(outdir, gen):
    """Returns an integer corresponding to the max. fit individual at generation 'gen'
    in the simulation located at 'outdir'"""
    if False == os.path.exists(outdir):
        print "\n. I can't find the simulation located at", outdir
    fin = open(outdir + "/FITNESS/fitness.gen" + gen.__str__() + ".txt", "r")
    lines = fin.readlines()
    fin.close()
    maxid = None
    minerr = None
    for l in lines:
        if l.__len__() > 2:
            tokens = l.split()
            id = int(tokens[0])
            err = float(tokens[2])
            if maxid == None:
                maxid = id
            if minerr == None:
                minerr = err
            if minerr > err:
                minerr = err
                maxid = id
    return maxid

#def get_reporters(outdir):
#    fin = open(outdir + "/DIR/pop.gen0.save.txt", "r")
    