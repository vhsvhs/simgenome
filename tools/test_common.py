#
# Common functions used by test modules.
#
import random
from diverging_colors import *


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