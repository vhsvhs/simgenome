"""
This script creates a text file with a list of batched simreg jobs.
You can then use mpi_dispatch to distribute the jobs to slave nodes."
Customize the variables at the top of the loop with the values
for your analysis.

USAGE:
python test_job_composer.py --run_id X

"""

import os, re, sys
from argparser import *
ap = ArgParser(sys.argv)
cpath = "run." + ap.getArg("--run_id") + ".sh"

#
# Customize these values. . .
#
fscalar = [-4]
pscalar = [0.1,0.01]
psammu = [0.01,0.05]
ursmu = [0.01,0.005]
growth = [1.0,0.2]
decay = [1.0,0.2]
MAXGEN = 50
POPSIZE = 30
NREPS = 3

def print_rules():
    fout = open( ap.getArg("--run_id") + ".rules", "w" )
    fout.write("RULE     0     5     4     1.0         0    1.0\n")
    fout.write("RULE     0     5     7     0.0001      1    1.0\n")
    fout.write("INPUT    0     0     0     7         0.5\n")
    fout.write("INPUT    0     1     0     7         0.5\n")
    fout.close()

print_rules()
fout = open(cpath, "w")
count = 0
for f in fscalar:
    for p in pscalar:
        for ps in psammu:
            for us in ursmu:
                for g in growth:
                    for d in decay:
                        for rep in range(0, NREPS):
                            count += 1
                            c = "simreg "
                            c += " --outdir ~/Applications/simgenome-c/experimental_runs/out." + ap.getArg("--run_id") + "." + count.__str__()
                            c += " --maxgen " + MAXGEN.__str__()
                            c += " --niid 6000 "
                            c += " --rulepath " + ap.getArg("--run_id") + ".rules"
                            c += " --popsize " + POPSIZE.__str__()
                            c += " --verbosity 5 "
                            c += " --maxgd 1 "
                            c += " --psam_mu " + ps.__str__()
                            c += " --urs_mu  " + us.__str__()
                            c += " --psamlenmu 0.0 "
                            c += " --f_scalar " + f.__str__()
                            c += " --pe_scalar " + p.__str__()
                            c += " --elite_prop 0.15 "
                            c += " --randompop"
                            c += " --growth_scalar " + g.__str__()
                            c += " --decay_scalar " + d.__str__()
                            fout.write(c + "\n")

fout.close()