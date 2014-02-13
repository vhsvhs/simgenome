"""
Plot a stereotype of an individual's circuit.

USAGE:
python plot_circuit.py --outdir X --gen X --id X

"""

import os, re, sys

from argparser import *
from test_common import *
from plot_includeme import *
ap = ArgParser(sys.argv)

outputdir = ap.getArg("--outdir") #outputdir is the folder into which a SimGenome run placed output.
generation = int( ap.getArg("--gen") )
id = int( ap.getArg("--id") )
rid =  int( ap.getArg("--rid") )
timeslice =  ap.getOptionalArg("--time")
if timeslice == False:
    timeslice = -1
else:
    timeslice = int(timeslice)



found_tfs = {}

def get_intx_submatrix(opath):
    tf_wt = {} # key = TF, value = interaction weight
    fin = open(opath, "r")
    last_time_seen = None
    last_rid_seen = None
    for l in fin.xreadlines():
        if l.startswith(". time"):
            tokens = l.split()
            last_time_seen = int(tokens[2])
            last_rid_seen = int(tokens[4])
        if l.startswith("site") and (timeslice == -1 or timeslice == last_time_seen) and (rid == last_rid_seen):
            #print last_time_seen
            tokens = l.split()
            toki = 3
            while(toki < tokens.__len__()):
                this_tf = int( tokens[toki] )
                this_mode = tokens[toki+1]
                this_p = float( tokens[toki+2] )
                if this_tf not in tf_wt:
                    tf_wt[this_tf] = 0.0
                if this_tf not in found_tfs.keys():
                    found_tfs[this_tf] = this_mode
                tf_wt[this_tf] += this_p
                toki += 3
    fin.close()
    return tf_wt

def get_intx_matrix():
    gene_tf_wt = {} # key = target gene, value = hash; key = TF, value = interaction weight
    
    fkey = "occ.gen" + generation.__str__() + ".id" + id.__str__() + ".gene"
    files = os.listdir(outputdir + "/OCCUPANCY")
    if files.__len__() == 0:
        print "\n. I found no occupancy logs for those settings."
    
    for f in files:
        if f.startswith(fkey):
            this_gene = re.sub(fkey, "", f)
            this_gene = re.sub(".txt", "", this_gene)
            this_gene = int( this_gene.split(".")[0] )
            inpath = outputdir + "/OCCUPANCY/" + f
            gene_tf_wt[this_gene] = get_intx_submatrix( inpath )
    return gene_tf_wt

m = get_intx_matrix()

if found_tfs.__len__() == 0:
    print "\n. I found no occupancy logs for those settings."
    exit()

tflist = found_tfs.keys()
tflist.sort()
line = "Gene\t"
for tf in tflist:
    line += tf.__str__() + found_tfs[tf] + "\t"
print line
for gene in m:
    line = gene.__str__() + "\t"
    for tf in tflist:
        line += "%.2f"%m[gene][tf] + "\t"
    print line


