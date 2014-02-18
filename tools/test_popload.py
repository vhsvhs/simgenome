"""
USAGE:
%> python tools/test_popload.2014.py --run_id X

. . . where X is the unique name of this run of the experiment.

"""

import sys, os
from test_common import *
from pwm_library import *
from argparser import *
ap = ArgParser(sys.argv)
RUNID = ap.getArg("--run_id")

if sys.argv.__len__() < 2:
    print "USAGE:"
    print "  $> python test_basic.py OUTDIR"
    print ""
    exit()

def get_command(id = "", random = False, poppath = None, popsize = 2, randseed = 12345, pe_scalar=0.001):
    OUTDIR = os.path.abspath( RUNID + id )
    command = "simreg "
    command += " --outdir " + OUTDIR
    #if poppath != None:
    #    command += " --run_clean "
    command += " --maxgen 11 "
    command += " --niid 50000 "
    command += " --no_sex"
    command += " --rulepath " + OUTDIR + "/" + RUNID + ".rules "
    command += " --popsize " + popsize.__str__()
    command += " --verbosity 5 "
    command += " --maxgd 1 "
    #command += " --nomu"
    command += " --psam_mu 0.1 " #1.0 " #0.05 "
    command += " --urs_mu 0.1 " #0.1 " #0.005"
    command += " --growth_scalar 1.0 "
    command += " --decay_scalar 1.0 "
    command += " --f_scalar -0.0005 "
    command += " --pe_scalar " + pe_scalar.__str__()
    #command += " --elite_prop 0.1 "
    command += " --randseed " + randseed.__str__() + " "
    if poppath == None and random == False:
        command += "--psampath " + OUTDIR + "/" + RUNID + ".psam "
        command += "--urspath " + OUTDIR + "/" + RUNID + ".urs "
    if random:
        command += " --randompop"
    if poppath != None:
        command += " --poppath " + poppath
    return command

def print_psams(outdir = ""):
    fout = open(outdir + "/" + RUNID + ".psam", "w")
    fout.write("# activation signal binds Gs\n")
    fout.write("Gene 0 1\n")
    fout.write("0.0 0.0 4.0 0.0\n")
    fout.write("0.0 0.0 4.0 0.0\n")
    fout.write("0.0 0.0 4.0 0.0\n")
    fout.write("0.0 0.0 4.0 0.0\n")
    fout.write("0.0 0.0 4.0 0.0\n")
    fout.write("0.0 0.0 4.0 0.0\n")

    fout.write("# non-specific repression\n")
    fout.write("Gene 1 0\n")
    fout.write("1 1 1 1\n")

    
def print_urss( outdir=""):
    fout = open(outdir + "/" + RUNID + ".urs", "w")
    
    # Transcription Factors:
    fout.write(">0\n")
    fout.write( "AAAAA\n")
    fout.write(">1\n")
    fout.write( "AAAAA\n")
    # A
    fout.write(">2\n")
    fout.write( get_repeat("A", 1) + "GGGGGG" + get_repeat("A", 1) + "\n")
    # B
    fout.write(">3\n")
    fout.write( get_repeat("A", 1) + "GGGGG" + get_repeat("A", 1) +  "\n")

    fout.close()

def print_rules( outdir = ""):
    fout = open( outdir + "/" + RUNID + ".rules", "w" )
    fout.write("# problem 0.\n")
    fout.write("INPUT 0     0     0     8      0.1\n")
    fout.write("INPUT 0     1     0     8      0.00001\n")
    fout.write("RULE  0     2     1     0.1      0    1.0\n")
    #fout.write("RULE  0     3     1     0.1      0    1.0\n")
    fout.close()



#################################
#
# MAIN part 1: run a normal SimGenome simulation. . .
#
commands = []

pe = 0.005
thisid = "erg"
OUTDIR = os.path.abspath( RUNID + thisid )
if os.path.exists(OUTDIR):
    os.system("rm -rf " + OUTDIR)
os.system( "mkdir " + OUTDIR )
print_rules(outdir = OUTDIR)
print_psams(outdir = OUTDIR)
print_urss(outdir = OUTDIR)
commands.append( get_command(id = thisid, popsize=4, randseed=123456789, pe_scalar=pe) )

ppath = OUTDIR + "/POPS/pop.gen10.save.txt"

thisid = "copy"
OUTDIR = os.path.abspath( RUNID + thisid )
if os.path.exists(OUTDIR):
    os.system("rm -rf " + OUTDIR)
os.system( "mkdir " + OUTDIR )
#ppath = None
print_rules(outdir = OUTDIR)
#print_psams(outdir = OUTDIR)
#print_urss(outdir = OUTDIR)
commands.append( get_command(id = thisid, popsize=4, randseed=123456789, pe_scalar=pe, poppath = ppath ) )

thisid = "copy2"
OUTDIR = os.path.abspath( RUNID + thisid )
if os.path.exists(OUTDIR):
    os.system("rm -rf " + OUTDIR)
os.system( "mkdir " + OUTDIR )
#ppath = None
print_rules(outdir = OUTDIR)
#print_psams(outdir = OUTDIR)
#print_urss(outdir = OUTDIR)
commands.append( get_command(id = thisid, popsize=4, randseed=123456789, pe_scalar=pe, poppath = ppath ) )


print commands

# This vesion builds random populations
#commands.append( get_command(id = "rand1", random = True, popsize=48, randseed=123456789) )

fout = open(OUTDIR + "/" + RUNID + ".commands.txt", "w")
for c in commands:
    fout.write(c + "\n")
fout.close()
os.system("source " + OUTDIR + "/" + RUNID + ".commands.txt" )

print ". Writing commands to " + OUTDIR + "/" + RUNID + ".commands.txt"

#os.system("python tools/plot_expression.py " + RUNID + "a1 --skip_binding True")

