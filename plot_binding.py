#
# USAGE:
# python plot_binding.py --configpath <PATH> --timeslice X
#
# . . . where <PATH> is a filepath to a config.* file
#

import re, sys, os
import matplotlib.pyplot as plt
import math
import networkx as nx
from argparser import *

######################################################################
ap = ArgParser(sys.argv)

path = ap.getArg("--configpath")
timeslice = int( ap.getArg("--timeslice") )
gene = int( ap.getArg("--gene") )

#
# Find the lines that are relevant to this timeslice and gene
#
foundit = False
fin = open(path, "r")
relevant_lines = []
for l in fin.readlines():
    if l.startswith(". TIME"):
        tokens = l.split()
        this_time = int(tokens[2])
        this_gene = int(tokens[4])
        if this_time == timeslice and this_gene == gene:
            foundit = True
            continue
        else:
            foundit = False
    if foundit == True and l.__len__() > 2:
        relevant_lines.append(l)
fin.close()

#
# Fill site_tf_expre
#
site_tf_expr = {}
for l in relevant_lines:
    print l
    tokens = l.split()
    this_site = int(tokens[1])
    tokens = tokens[2:]
    i = 0
    while i < tokens.__len__():
        this_tf = int(tokens[i])
        this_p = float(tokens[3])
        if this_site not in site_tf_expr:
            site_tf_expr[this_site] = {}
        site_tf_expr[this_site][this_tf] = this_p
        i += 6

for site in site_tf_expr:
    print site, site_tf_expr[site]
     