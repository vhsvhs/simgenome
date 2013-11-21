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
tf_site_expr = {}
for l in relevant_lines:
    print l
    tokens = l.split()
    this_site = int(tokens[1])
    tokens = tokens[2:]
    i = 0
    while i < tokens.__len__():
        this_tf = int(tokens[i])
        this_p = float(tokens[i+3])
        if this_tf not in tf_site_expr:
            tf_site_expr[this_tf] = {}
        tf_site_expr[this_tf][this_site] = this_p
        
        if this_site not in site_tf_expr:
            site_tf_expr[this_site] = {}
        site_tf_expr[this_site][this_tf] = this_p
        
        i += 6

#
# for testing:
#
#for site in site_tf_expr:
#    print site, site_tf_expr[site]

#for tf in tf_site_expr:
#    print tf, tf_site_expr[tf]

maxy = 1.0
miny = 0.0
maxx = 0 # N sites
minx = 0
lout = "" # line out

# scan for the max timepoint
for tf in tf_site_expr:
    for site in tf_site_expr[tf]:
        if site > maxx:
            maxx = site

for tf in tf_site_expr:
    l = "tf" + tf.__str__() + " <-c("    
    lastval = None
    for site in range(1, maxx+1):
        if site in tf_site_expr[tf]:
            lastval = tf_site_expr[tf][site]
            l += "%.3f"%tf_site_expr[tf][site]
        else:
            lastval = lastval / 2 # decay function
            l += "%.3f"%lastval
        l += ","
    l = re.sub(",$", "", l)
    l += ");\n"
    lout += l

l = "x <-c("
for site in range(1, maxx+1):
    l += site.__str__() + ","
l = re.sub(",$", "", l)
l += ");\n"
lout += l

lout += "pdf('" + path + ".pdf', width=7, height=4);\n"
title = "Binding Probabilities"
xlab = "Sites"
ylab = "P"
lout += "plot(c(1," + maxx.__str__() + "), c(0," + maxy.__str__() + "), type='n', main='" + title + "', xlab='" + xlab + "', ylab='" + ylab + "');\n"

for tf in tf_site_expr:
    if tf == 5:
        lout += "points(x, tf" + tf.__str__() + ",type='l', col=" + tf.__str__() + ", lwd=1);\n"

fout_cran = open(path + ".cran", "w")
fout_cran.write(lout)
fout_cran.write("dev.off();\n")
fout_cran.close()

os.system("r --no-save < " + path + ".cran")

