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
    command += " --niid 10000 "
    command += " --rulepath " + OUTDIR + "/" + RUNID + ".rules "
    command += " --popsize " + popsize.__str__()
    command += " --verbosity 5 "
    command += " --maxgd 1 "
    command += " --psam_mu 0.05 "
    command += " --urs_mu 0.005"
    command += " --growth_scalar 1.0 "
    command += " --decay_scalar 1.0 "
    command += " --f_scalar -2.0 "
    command += " --pe_scalar " + pe_scalar.__str__()
    command += " --elite_prop 0.29 "
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
    fout.write("# activator binds GGGG\n")
    fout.write("Gene 0 1\n")
    fout.write("0.0 0.1 30.0 1\n")
    fout.write("0.0 0.1 30.0 1\n")
    fout.write("0.0 0.1 30.0 1\n")
    fout.write("0.0 0.1 30.0 1\n")

    fout.write("# activator binds CCCC\n")
    fout.write("Gene 1 1\n")
    fout.write("0.0 30 0.1 0.1\n")
    fout.write("0.0 30 0.1 0.1\n")
    fout.write("0.0 30 0.1 0.1\n")
    fout.write("0.0 30 0.1 0.1\n")
    
    fout.write("# repressor binds TTTT\n")
    fout.write("Gene 2 0\n")
    fout.write("0.0 0.1 0.1 80.0\n")
    fout.write("0.0 0.1 0.1 20.0\n")
    fout.write("0.0 0.1 0.1 20.0\n")
    fout.write("0.0 0.1 0.1 20.0\n")
    fout.write("0.0 0.1 0.1 20.0\n")
    
    fout.write("# non-spec repression\n")
    fout.write("Gene 3 0\n")
    fout.write("1 1 1 1\n")
    fout.write("1 1 1 1\n")
    fout.write("1 1 1 1\n")
    
    fout.write("# non-spec activation\n")
    fout.write("Gene 4 1\n")
    fout.write("1 1 1 1\n")
    fout.write("1 1 1 1\n")
    fout.write("1 1 1 1\n")
    fout.close()
    
def print_urss( outdir=""):
    fout = open(outdir + "/" + RUNID + ".urs", "w")
    fout.write(">0\n")
    fout.write( get_repeat("A", 20) + "\n")
    fout.write(">1\n")
    fout.write( get_repeat("A", 20) + "\n")
    fout.write(">2\n")
    fout.write( get_repeat("A", 20) + "GGGG" + "\n")
    fout.write(">3\n")
    fout.write( get_repeat("A", 20) + "\n")
    fout.write(">4\n")
    fout.write( get_repeat("A", 20) + "\n")
    fout.write(">5\n")
    fout.write( "AAAAAAGGGGAAAAAATTTTTTAAAA\n")
   
    fout.close()

def print_rules( outdir = ""):
    fout = open( outdir + "/" + RUNID + ".rules", "w" )
    fout.write("RULE 0     5     4     1.0     0    1.0\n")
    fout.write("RULE 0     5     7     0.0001  1    1.0\n")
    fout.write("INPUT 0     0     0     7     0.1\n")
    fout.write("INPUT 0     1     0     7     0.0001\n")
    #fout.write("INPUT 0     2     0     7     0.0001\n")
    fout.write("INPUT 0     3     0     7     0.0001\n")
    fout.write("INPUT 0     4     0     7     0.0001\n")
    #fout.write("INPUT 0     1     0     7     0.1\n")
    fout.close()


#################################
#
# main
#
commands = []

pe_scalars = [0.005]
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

