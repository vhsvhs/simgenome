"""
Use this script as a basic test of simreg functionality.  It is mainly intended
to test biochemical/transcriptional aspects of the unified model, rather than
mutational/evolutionary properties.

USAGE:
%> python tools/test_basic.py --run_id X

OUTPUT:
A directory at X will be created, and all output from the analysis will be written into X.
This means that you should be careful about the directory in which you launch this analysis.

NOTE: This script requires that you install 'simreg' as a system-wide link.
"""


import sys, os
from test_common import *
from argparser import *
ap = ArgParser(sys.argv)
RUNID = ap.getArg("--run_id")

if sys.argv.__len__() < 2:
    print "USAGE:"
    print "  $> python test_basic.py OUTDIR"
    print ""
    exit()

def get_command(id = "", random = False, popsave = False, popsize = 2, randseed = 12345, pe_scalar=0.001):
    OUTDIR = os.path.abspath( RUNID + id )
    command = "simreg "
    command += " --outdir " + OUTDIR
    command += " --maxgen 1 "
    command += " --niid 15000 "
    command += " --rulepath " + OUTDIR + "/" + RUNID + ".rules "
    command += " --popsize " + popsize.__str__()
    command += " --verbosity 5 "
    command += " --maxgd 1 "
    command += " --psam_mu 0.05 " #1.0 " #0.05 "
    command += " --urs_mu 0.005 " #0.1 " #0.005"
    command += " --growth_scalar 1.0 "
    command += " --decay_scalar 1.0 "
    command += " --f_scalar -2.0 "
    command += " --pe_scalar " + pe_scalar.__str__()
    command += " --elite_prop 0.29 "
    #command += " --tran_cdf "
    command += " --randseed " + randseed.__str__() + " "
    if popsave == False and random == False:
        command += "--psampath " + OUTDIR + "/" + RUNID + ".psam "
        command += "--urspath " + OUTDIR + "/" + RUNID + ".urs "
    if random:
        command += " --randompop"
    if popsave:
        command += " --poppath pop.gen0.save.backup"
    return command

def print_psams(outdir = ""):
    fout = open(outdir + "/" + RUNID + ".psam", "w")
    fout.write("# inducer signal binds G\n")
    fout.write("Gene 0 0\n")
    # dec4:
    #fout.write("0.0 0.0 80 0\n")    
    # dec2:
    fout.write("0.0 0.0 7.0 0\n")
    fout.write("0.0 0.0 8.0 0\n")
    fout.write("0.0 0.0 7.0 0\n")
    fout.write("0.0 0.0 8.0 0\n")
    fout.write("0.0 0.0 7.0 0\n")
    fout.write("0.0 0.0 8.0 0\n")

    fout.write("# non-spec activation\n")
    fout.write("Gene 1 1\n")
    fout.write("1 1 1 1\n")
    fout.write("1 1 1 1\n")
    fout.write("1 1 1 1\n")
    fout.write("1 1 1 1\n")

    fout.write("# repressor binds TTTT\n")
    fout.write("Gene 2 0\n")
    fout.write("0.0 0.0 0.0 18\n")
    fout.write("0.0 0.0 0.0 18\n")
    fout.write("0.0 0.0 0.0 10\n")
    fout.write("0.0 0.0 0.0 10\n")
    fout.write("0.0 0.0 0.0 10\n")

    
def print_urss( outdir=""):
    fout = open(outdir + "/" + RUNID + ".urs", "w")
    fout.write(">0\n")
    fout.write( "AAAAA\n")
    fout.write(">1\n")
    fout.write( "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n")
    fout.write(">2\n")
    fout.write( "AAAAAAAAAGGGGGGAAAAAAAAAAAAAAAAAA\n")
    fout.write(">3\n")
    fout.write( "AAAAAAAAATTAAAAAAAAAAAAAAAAAAAAAA\n")
    fout.write(">4\n")
    fout.write( "AAAAAAAAATTTAAAAAAAAAAAAAAAAAAAAA\n")
    fout.write(">5\n")
    fout.write( "AAAAAAAAATTTTAAAAAAAAAAAAAAAAAAAA\n")
    fout.write(">6\n")
    fout.write( "AAAAAAAAATTTTTAAAAAAAAAAAAAAAAAAA\n")
    fout.write(">7\n")
    fout.write( "AAAAAAAAATTTTTTAAAAAAAAAAAAAAAAAA\n")
    fout.write(">8\n")
    fout.write( "AAAAAAAAATTTTTTTAAAAAAAAAAAAAAAAA\n")
    fout.write(">9\n")
    fout.write( "AAAAAAAAATTTTTTTTAAAAAAAAAAAAAAAA\n")
    fout.write(">10\n")
    fout.write( "AAAAAAAAATTTTTTTTTAAAAAAAAAAAAAAA\n")
    fout.close()

def print_rules( outdir = ""):
    fout = open( outdir + "/" + RUNID + ".rules", "w" )
    fout.write("# inducer pulse\n")
    fout.write("INPUT 0     0     0     10      0.00001\n")
    fout.write("INPUT 0     0     11     12     0.0001\n")
    fout.write("INPUT 0     0     12     13     0.001\n")
    fout.write("INPUT 0     0     13     14     0.01\n")
    fout.write("INPUT 0     0     14     25     0.1\n")
    fout.write("INPUT 0     0     25     27     0.01\n")
    fout.write("INPUT 0     0     27     29     0.001\n")
    fout.write("INPUT 0     0     29     31     0.0001\n")
    fout.write("INPUT 0     0     31     50     0.00001\n")
    
    fout.write("# low-level activation always on\n")
    fout.write("INPUT 0     1     0     50     0.00047\n")


    #fout.write("RULE 0     5     4     1.0     0    1.0\n")
    #fout.write("RULE 0     5     7     0.0001  1    1.0\n")
    fout.close()

# analyze for Kd
def analyze_kd(outdir = ""):
    """THis function gets data for a kinetics experiment on Dec. 5 2013.
    It contains some ad-hoc things, and should be considered depricated."""
    print "\n================================================="
    print ". Analyzing for Kd..."
    fin = open(outdir + "/LOGS/expression.txt", "r")

    # dec2:
    timepoint = 32
    print ". Expression at Time", timepoint,"..."
    gene_expr = {}
    for l in fin.xreadlines():
        if l.__len__() > 10:
            tokens = l.split()
            time = int(tokens[5])
            if time == timepoint:
                gene = int(tokens[9])
                if tokens[11].startswith("expr"):
                    gene_expr[ gene ] = float(tokens[12])
                else:
                    gene_expr[ gene ] = float(tokens[11])
    for gene in gene_expr:
        print gene_expr[gene]
    fin.close()
    
    gene_pgene2 = {}
    for gene in gene_expr:
        #if gene == 0:
        #    continue
        gene_pgene2[gene] = 0.0
        fin = open(outdir + "/OCCUPANCY/occ.gen0.id0.gene" + gene.__str__() + ".txt", "r")
        foundit = False
        for l in fin.xreadlines():
            #print gene, l
            if l.__len__() > 10:
                if l.__contains__("time " + timepoint.__str__()):
                    foundit = True
                elif l.__contains__("time") and False == l.__contains__("time " + timepoint.__str__()):
                    foundit = False
                elif foundit:
                    tokens = l.split()
                    #print tokens
                    pgene2 = float(tokens[11])        
                    if gene_pgene2[gene] < pgene2:
                        gene_pgene2[gene] = pgene2
        fin.close()
    print ". Occupancy P(gene2)..."
    for gene in gene_expr:
        print gene_pgene2[gene]

#################################
#
# main
#
commands = []

pe_scalars = [0.01]
count = 0
for pe in pe_scalars:
    count += 1
    thisid = "a" + count.__str__()
    OUTDIR = os.path.abspath( RUNID + thisid )
    if False == os.path.exists(OUTDIR):
        os.system( "mkdir " + OUTDIR )
    print_rules(outdir = OUTDIR)
    print_psams(outdir = OUTDIR)
    print_urss(outdir = OUTDIR)

    # This version uses the pre-built PSAMs and URSs
    commands.append( get_command(id = thisid, popsize=1, randseed=123456789, pe_scalar=pe) )

# This vesion builds random populations
#commands.append( get_command(id = "rand1", random = True, popsize=48, randseed=123456789) )

fout = open(OUTDIR + "/" + RUNID + ".commands.txt", "w")
for c in commands:
    fout.write(c + "\n")
fout.close()
os.system("source " + OUTDIR + "/" + RUNID + ".commands.txt" )

# post analysis:
analyze_kd(outdir = OUTDIR)

