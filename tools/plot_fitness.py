#
# Input: the captured STDOUT stream from a sim-reg run.
# Output: an R plot of generations
#
import os, re, sys

fin = open(sys.argv[1], "r")

# data is  key = gene id, value = hash: key = time, value = expression level

tarr = []
color = {}

data = {}

def color_for_run(x):
    if x not in color:
        this_color = (color.__len__()+1).__str__()
        color[x] = this_color
    return color[x]