#
# Measure H1 versus H2. See my lab notes for more info.
#
# Compare the K logs for two individuals from different generations.
#

import os, re, sys

from argparser import *
from test_common import *
from plot_includeme import *
ap = ArgParser(sys.argv)

indir = ap.getArg("--dir") #outputdir is the folder into which a SimGenome run placed output.
gen1 = int( ap.getArg("--gen1") )
gen2 = int( ap.getArg("--gen2") )
id1 = int( ap.getArg("--id1") )
id2 = int( ap.getArg("--id2") )
rids =  ap.getList("--rids") # desired regulatory problems
for i in range(0, rids.__len__()):
    rids[i] = int(rids[i])
times =  ap.getList("--times") # desired time slices to compare
for i in range(0, times.__len__()):
    times[i] = int(times[i])
reporters = ap.getList("--reporters")
for i in range(0, reporters.__len__()):
    reporters[i] = int(reporters[i])
regulators = ap.getList("--regulators")
for i in range(0, regulators.__len__()):
    regulators[i] = int(regulators[i])

def get_kpath(gen, id, gene, rid):
    return indir + "/OCCUPANCY/k.gen" + gen.__str__() + ".id" + id.__str__() + ".gene" + gene.__str__() + ".rid" + rid.__str__() + ".txt"


def read_klog(kpath):
    fin = open(kpath, "r")
    lines = fin.readlines()
    fin.close()
    time_tf_val = {}
    for l in lines:
        if l.__len__() > 2:
            tokens = l.split()
            time = int(re.sub("\:", "", tokens[0]))
            time_tf_val[time] = {}
            for ii in range(1, tokens.__len__()):
                val = float(tokens[ii])
                time_tf_val[time][ii-1] = val
    return time_tf_val
            

def compare_klogs(kpath1, kpath2):
    """Returns the distance between two K log matrices.
    A negative value indicates that the gene associated with kpath2
    has weaker regulatory interactions than the same gene in kpath1.
    A positive value indicates that the gene associated with kpath2
    gained stronger regulatory interactions compared to that gene
    in kpath1."""
    m1 = read_klog(kpath1)
    m2 = read_klog(kpath2)
    d = 0.0 #  distance between values
    for t in times:
        for tf in m1[t]:
            d += m2[t][tf] - m1[t][tf]
    return d

#
# Verify the files exist
#
kpaths = []
for rid in rids:
    for rep in reporters:
        kpaths.append( get_kpath( gen1, id1, rep, rid ) )
        kpaths.append( get_kpath( gen2, id2, rep, rid ) )
    for reg in regulators:
        kpaths.append( get_kpath( gen1, id1, reg, rid ) )
        kpaths.append( get_kpath( gen2, id2, reg, rid ) )
for kpath in kpaths:
    if False == os.path.exists(kpath):
        print "\n. I can't find the K log at", kpath
        exit()

sumd_reps = 0.0
for rid in rids:
    for rep in reporters:
        kpath1 = get_kpath( gen1, id1, rep, rid ) 
        kpath2 = get_kpath( gen2, id2, rep, rid )
        d = compare_klogs(kpath1, kpath2)
        sumd_reps += d

sumd_regs = 0.0
for rid in rids:
    for reg in regulators:
        kpath1 = get_kpath( gen1, id1, reg, rid ) 
        kpath2 = get_kpath( gen2, id2, reg, rid )
        d = compare_klogs(kpath1, kpath2)
        sumd_regs += d

mean_repd = sumd_reps / reporters.__len__()
mean_regd = sumd_regs / regulators.__len__()
print "mean rep. d = ", mean_repd
print "mean reg. d = ", mean_regd
    