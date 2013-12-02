#
# Common functions used by test modules.
#
import random

def get_rand_nt(len):
    nts = ["A", "C", "G", "T"]
    ret = ""
    for i in range(0, len):
        ret += random.sample(nts, 1)[0]
    return ret