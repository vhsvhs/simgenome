##############################################################
#
#
#
#
#
# USAGE:
# %> python tools/test_mu.py P
# . . . where P is the path to a desired output directory
#
# 
# Screens values of mutation rates versus N generations to solution,
# for the simple problem of assembling a feed-forward loop.
#
# NOTE: This script requires that you install simreg as an exe available from anywhere on the path.
#
import sys, os
from test_common import *

if sys.argv.__len__() < 2:
    print "USAGE:"
    print "  $> python test_mu.py OUTDIR"
    print ""
    exit()

OUTDIR = sys.argv[1]
#OUTDIR = os.path.abspath(OUTDIR)
PSAM = [0.1]
URS = [0.005,0.01,0.05,0.1]
nreps = 1

nprocs = PSAM.__len__() * URS.__len__() * nreps # how many MPI processes?

def make_runid(psammu, ursmu, rep):
    return "p" + psammu.__str__() + "-u" + ursmu.__str__() + "-r" + rep.__str__()

def make_directory(psammu, ursmu, rep):
    runid = make_runid(psammu, ursmu, rep)
    print OUTDIR + "/out." + runid
    if os.path.exists(OUTDIR + "/out." + runid):
        os.system("rm -rf " + OUTDIR + "/out." + runid)
    os.system("mkdir " + OUTDIR + "/out." + runid)        

def get_command(psammu, ursmu, rep):
    runid = make_runid(psammu, ursmu, rep)
    command = "simreg "
    command += " --outdir " + OUTDIR + "/out." + runid + " " 
    command += "--psampath " + OUTDIR + "/test_mu.psam "
    command += "--urspath " + OUTDIR + "/test_mu.urs "
    command += "--rulepath " + OUTDIR + "/test_mu.rules "
    command += "--maxgen 500 "
    command += "--niid 10000 "
    command += "--popsize 10 "
    command += "--verbosity 21 "
    command += "--maxgd 1 "
    command += " --urs_mu " + ursmu.__str__()
    command += " --psam_mu " + psammu.__str__()
    command += " --psamlenmu 0.1"
    command += " --pe_scalar 0.001"
    command += " --growth_scalar 0.3"
    command += " --decay_scalar 0.3"
    return command
def print_psams():
    fout = open(OUTDIR + "/test_mu.psam", "w")
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
    
def print_urss():
    fout = open(OUTDIR + "/test_mu.urs", "w")
    fout.write(">0\n")
    fout.write("AAAAA\n")
    fout.write(">1\n")
    fout.write( get_rand_nt(1000) + "\n")
    #fout.write("AAAAAGGGGGAAAAAAAAACCCCAAAAAAAAACCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n")
    fout.write(">2\n")
    fout.write( get_rand_nt(1000) + "\n" )
    #fout.write("AAAAAAAAAACCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n")
    fout.write(">3\n")
    fout.write("AAAA\n")
    fout.write(">4\n")
    fout.write("AAAA\n")
    fout.write(">5\n")
    fout.write( get_rand_nt(1000) + "\n" )
    #fout.write("AAAACCCCCCAAAAAAAAAAATTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n")
   
    fout.close()

def print_rules():
    fout = open( OUTDIR + "/test_mu.rules", "w" )
    fout.write("RULE 0 5 3 0.1 0 8.0\n")
    fout.write("RULE 0 5 8 0.0001 1 1.0\n")
    fout.write("INPUT 0 0 0 1 0.6\n")
    fout.write("INPUT 0 0 2 10 0.00001\n")
    fout.write("INPUT 0 3 0 10 0.5\n")
    fout.write("INPUT 0 4 0 10 0.4\n")
    fout.close()
    
    
    
def write_get_data_tool():
    fout = open(OUTDIR + "/get_data.py", "w")
    fout.write("import os, sys\n")
    fout.write("for f in os.listdir(\"./\"):\n")
    fout.write("    if os.path.isdir(f):\n")
    fout.write("        genlog = f + \"/LOGS/generations.txt\"\n")
    fout.write("        if False == os.path.exists(genlog):\n")
    fout.write("            print \"I can't find the file\", genlog\n")
    fout.write("            continue\n")
    fout.write("        fin = open(genlog, \"r\")\n")
    fout.write("        lines = fin.readlines()\n")
    fout.write("        lastline = lines[ lines.__len__()-1 ]\n")
    fout.write("        tokens = lastline.split()\n")
    fout.write("        gen = tokens[1]\n")
    fout.write("        max = tokens[3]\n")
    fout.write("        min = tokens[5]\n")
    fout.write("        mean = tokens[7]\n")
    fout.write("        print f, gen, max, min, mean\n")
    fout.write("        fin.close()\n")
    fout.close()

#################################
#
# main
#

if False == os.path.exists( OUTDIR ):
    os.system("mkdir " + OUTDIR)
print_rules()
print_psams()
print_urss()
commands = []
for ii in PSAM:
    for jj in URS:
        print ii, jj
        for rep in range(0, nreps):
            make_directory(ii, jj, rep)
            c = get_command(ii, jj, rep)
            commands.append( c )


fout = open(OUTDIR + "/test_mu.commands.txt", "w")
for c in commands:
    fout.write(c + "\n")
fout.close()

print "I wrote a BATCH script to 'test_mu.commands.txt'"

print "Now run the following command, or some variant that will work on your machine:"
print "mpirun -np " + nprocs.__str__() + " --machinefile hosts.txt /common/bin/mpi_dispatch test_mu.commands.txt"

write_get_data_tool()
print "When that completes, run the script " + OUTDIR + "/get_data.py"