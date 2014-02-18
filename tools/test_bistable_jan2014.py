"""
This script performs the experiment. . 

TBD.

USAGE:
%> python tools/test_jan15.2014.py --run_id X

. . . where X is the unique name of this run of the experiment.

OUTPUT:
A directory at X will be created, and all output from the analysis will be written into X.

NOTE: This script requires that you install 'simreg' and make it available in your PATH
environment variable.
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

def get_command(id = "", random = False, popsave = False, popsize = 2, randseed = 12345, pe_scalar=0.001):
    OUTDIR = os.path.abspath( RUNID + id )
    command = "simreg "
    command += " --outdir " + OUTDIR
    command += " --maxgen 75 "
    command += " --niid 8500 "
    command += " --rulepath " + OUTDIR + "/" + RUNID + ".rules "
    command += " --popsize " + popsize.__str__()
    command += " --verbosity 5 "
    command += " --maxgd 1 "
    command += " --psam_mu 0.05 " #1.0 " #0.05 "
    command += " --urs_mu 0.01 " #0.1 " #0.005"
    command += " --growth_scalar 1.0 "
    command += " --decay_scalar 1.0 "
    command += " --f_scalar -0.0005 "
    command += " --pe_scalar " + pe_scalar.__str__()
    command += " --elite_prop 0.1 "
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
    fout.write("# activation signal binds Gs\n")
    fout.write("Gene 0 1\n")
    fout.write("0.0 0.0 4.0 0.0\n")
    fout.write("0.0 0.0 4.0 0.0\n")
    fout.write("0.0 0.0 4.0 0.0\n")
    fout.write("0.0 0.0 4.0 0.0\n")
    fout.write("0.0 0.0 4.0 0.0\n")
    fout.write("0.0 0.0 4.0 0.0\n")

    fout.write("# activation signal binds Ts\n")
    fout.write("Gene 1 1\n")
    fout.write("0.0 0.0 0.0 4.0\n")
    fout.write("0.0 0.0 0.0 4.0\n")
    fout.write("0.0 0.0 0.0 4.0\n")
    fout.write("0.0 0.0 0.0 4.0\n")
    fout.write("0.0 0.0 0.0 4.0\n")
    fout.write("0.0 0.0 0.0 4.0\n")

    fout.write("# non-specific repression\n")
    fout.write("Gene 2 0\n")
    fout.write("1 1 1 1\n")

    fout.write("# Gene A is a represser, binds CCTCCG\n")
    fout.write("Gene 3 0\n")
    fout.write("0.0 4.0 0.0 0.0\n")
    fout.write("0.0 4.0 0.0 0.0\n")
    fout.write("0.0 0.0 0.0 4.0\n")
    fout.write("0.0 4.0 0.0 0.0\n")
    fout.write("0.0 4.0 0.0 0.0\n")
    fout.write("0.0 0.0 4.0 0.0\n")

    fout.write("# Gene B is an activator, TATAAA\n")
    fout.write("Gene 4 1\n")
    fout.write("0.0 0.0 0.0 3\n")
    fout.write("3.0 0.0 0.0 0\n")
    fout.write("0.0 0.0 0.0 3\n")
    fout.write("3.0 0.0 0.0 0\n")
    fout.write("3.0 0.0 0.0 0\n")
    fout.write("3.0 0.0 0.0 0\n")

    fout.write("# Gene C is a represser, bindg CGTTAC\n")
    fout.write("Gene 5 0\n")
    fout.write("0.0 4.0 0.0 0.0\n")
    fout.write("0.0 0.0 4.0 0.0\n")
    fout.write("0.0 0.0 0.0 4.0\n")
    fout.write("0.0 0.0 0.0 4.0\n")
    fout.write("4.0 0.0 0.0 0.0\n")
    fout.write("0.0 4.0 0.0 0.0\n")

    fout.write("# Gene D is an activator, binds CTGGAC\n")
    fout.write("Gene 6 1\n")
    fout.write("0.0 3.0 0.0 0\n")
    fout.write("0.0 0.0 0.0 3\n")
    fout.write("0.0 0.0 3.0 0\n")
    fout.write("0.0 0.0 3.0 0\n")
    fout.write("3.0 0.0 0.0 0\n")
    fout.write("0.0 3.0 0.0 0\n")
    
    fout.write("Gene 7 0\n")
    fout.write( get_random_pwm(6) )

    fout.write("Gene 8 0\n")
    fout.write( get_random_pwm(6) )
    
    fout.write("Gene 9 1\n")
    fout.write( get_random_pwm(6) )
    
    fout.write("Gene 10 1\n")
    fout.write( get_random_pwm(6) )

    fout.write("Gene 11 1\n")
    fout.write( get_random_pwm(6) )
    
    fout.write("Gene 12 1\n")
    fout.write( get_random_pwm(6) )
    
def print_urss( outdir=""):
    fout = open(outdir + "/" + RUNID + ".urs", "w")
    
    # Transcription Factors:
    fout.write(">0\n")
    fout.write( "AAAAA\n")
    fout.write(">1\n")
    fout.write( "AAAAA\n")
    fout.write(">2\n") # background repression
    fout.write( "AAAAA\n")
    
    # A
    fout.write(">3\n")
    fout.write( get_rand_nt(100) + "GGGGGGGG" + get_rand_nt(192) + "\n")
    # B
    fout.write(">4\n")
    fout.write( get_rand_nt(100) + "GGGGG" + get_rand_nt(100) + "CCTCCG" + get_rand_nt(93) +  "\n")
    # C
    fout.write(">5\n")
    fout.write( get_rand_nt(100) + "TTTTTTTT" + get_rand_nt(192) + "\n")
    # D
    fout.write(">6\n")
    fout.write( get_rand_nt(100) + "TTTTT" + get_rand_nt(100) + "CGTTAC" + get_rand_nt(93) + "\n")
    # E
    fout.write(">7\n")
    fout.write( get_rand_nt(300) + "\n")
    # F
    fout.write(">8\n")
    fout.write( get_rand_nt(300) + "\n")
    # G
    fout.write(">9\n")
    fout.write( get_rand_nt(300) + "\n")
    # H
    fout.write(">10\n")
    fout.write( get_rand_nt(300) + "\n")
    # G
    fout.write(">11\n")
    fout.write( get_rand_nt(300) + "\n")
    # H
    fout.write(">12\n")
    fout.write( get_rand_nt(300) + "\n")
    
    # A Reporters:    
    fout.write(">13\n")
    fout.write( get_rand_nt(100) + "TATAAA" + get_rand_nt(4) + "TATAAA" + get_rand_nt(190) +  "\n")
    fout.write(">14\n")
    fout.write( get_rand_nt(100) + "TATAAA" + get_rand_nt(4) + "TATAAA" + get_rand_nt(190) + "\n")
    fout.write(">15\n")
    fout.write( get_rand_nt(100) + "TATAAA" + get_rand_nt(4) + "TATAAA" + get_rand_nt(190) + "\n")
    fout.write(">16\n")
    fout.write( get_rand_nt(100) + "TATAAA" + get_rand_nt(4) + "TATAAA" + get_rand_nt(190) + "\n")
    
    # B Reporters:        
    fout.write(">17\n")
    fout.write( get_rand_nt(100) + "CTGGAC" + get_rand_nt(4) + "CTGGAC" + get_rand_nt(190) +  "\n")
    fout.write(">18\n")
    fout.write( get_rand_nt(100) + "CTGGAC" + get_rand_nt(4) + "CTGGAC" + get_rand_nt(190) + "\n")
    fout.write(">19\n")
    fout.write( get_rand_nt(100) + "CTGGAC" + get_rand_nt(4) + "CTGGAC" + get_rand_nt(190) + "\n")
    fout.write(">20\n")
    fout.write( get_rand_nt(100) + "CTGGAC" + get_rand_nt(4) + "CTGGAC" + get_rand_nt(190) + "\n")
    fout.close()

def print_rules( outdir = ""):
    fout = open( outdir + "/" + RUNID + ".rules", "w" )
    fout.write("# problem 0.\n")
    fout.write("INPUT 0     0     0     8      0.1\n")
    fout.write("INPUT 0     1     0     8      0.00001\n")
    fout.write("INPUT 0     2     0     8     0.0005\n")

    fout.write("RULE  0     13     3     0.01      0    1.0\n")
    fout.write("RULE  0     14     3     0.01      0    1.0\n")
    fout.write("RULE  0     15     3     0.01      0    1.0\n")
    fout.write("RULE  0     16     3     0.01      0    1.0\n")
    
    fout.write("RULE  0     13     4     0.1      0    1.0\n")
    fout.write("RULE  0     14     4     0.1      0    1.0\n")
    fout.write("RULE  0     15     4     0.1      0    1.0\n")
    fout.write("RULE  0     16     4     0.1      0    1.0\n")

    fout.write("RULE  0     13     5     0.1      0    1.0\n")
    fout.write("RULE  0     14     5     0.1      0    1.0\n")
    fout.write("RULE  0     15     5     0.1      0    1.0\n")
    fout.write("RULE  0     16     5     0.1      0    1.0\n")

    fout.write("RULE  0     13     7     0.01      1    1.0\n")
    fout.write("RULE  0     14     7     0.01      1    1.0\n")
    fout.write("RULE  0     15     7     0.01      1    1.0\n")
    fout.write("RULE  0     16     7     0.01      1    1.0\n")
    
    fout.write("RULE  0     13     8     0.0001      1    1.0\n")
    fout.write("RULE  0     14     8     0.0001      1    1.0\n")
    fout.write("RULE  0     15     8     0.0001      1    1.0\n")
    fout.write("RULE  0     16     8     0.0001      1    1.0\n")

    fout.write("RULE  0     17     1     0.0001      1    0.25\n")
    fout.write("RULE  0     18     1     0.0001      1    0.25\n")
    fout.write("RULE  0     19     1     0.0001      1    0.25\n")
    fout.write("RULE  0     20     1     0.0001      1    0.25\n")
    
    fout.write("RULE  0     17     3     0.0001      1    0.25\n")
    fout.write("RULE  0     18     3     0.0001      1    0.25\n")
    fout.write("RULE  0     19     3     0.0001      1    0.25\n")
    fout.write("RULE  0     20     3     0.0001      1    0.25\n")
    
    fout.write("RULE  0     17     5     0.0001      1    0.25\n")
    fout.write("RULE  0     18     5     0.0001      1    0.25\n")
    fout.write("RULE  0     19     5     0.0001      1    0.25\n")
    fout.write("RULE  0     20     5     0.0001      1    0.25\n")

    fout.write("RULE  0     17     7     0.0001      1    0.25\n")
    fout.write("RULE  0     18     7     0.0001      1    0.25\n")
    fout.write("RULE  0     19     7     0.0001      1    0.25\n")
    fout.write("RULE  0     20     7     0.0001      1    0.25\n")


    fout.write("# problem 1.\n")
    fout.write("INPUT 1     0     0     8      0.00001\n")
    fout.write("INPUT 1     1     0     8      0.1\n")
    fout.write("INPUT 1     2     0     8     0.0005\n")

    fout.write("RULE  1     17     3     0.01      0    1.0\n")
    fout.write("RULE  1     18     3     0.01      0    1.0\n")
    fout.write("RULE  1     19     3     0.01      0    1.0\n")
    fout.write("RULE  1     20     3     0.01      0    1.0\n")
    
    fout.write("RULE  1     17     4     0.1      0    1.0\n")
    fout.write("RULE  1     18     4     0.1      0    1.0\n")
    fout.write("RULE  1     19     4     0.1      0    1.0\n")
    fout.write("RULE  1     20     4     0.1      0    1.0\n")

    fout.write("RULE  1     17     5     0.1      0    1.0\n")
    fout.write("RULE  1     18     5     0.1      0    1.0\n")
    fout.write("RULE  1     19     5     0.1      0    1.0\n")
    fout.write("RULE  1     20     5     0.1      0    1.0\n")

    fout.write("RULE  1     17     7     0.01      1    1.0\n")
    fout.write("RULE  1     18     7     0.01      1    1.0\n")
    fout.write("RULE  1     19     7     0.01      1    1.0\n")
    fout.write("RULE  1     20     7     0.01      1    1.0\n")
    
    fout.write("RULE  1     17     8     0.0001      1    1.0\n")
    fout.write("RULE  1     18     8     0.0001      1    1.0\n")
    fout.write("RULE  1     19     8     0.0001      1    1.0\n")
    fout.write("RULE  1     20     8     0.0001      1    1.0\n")
    
    fout.write("RULE  1     13     1     0.0001      1    0.25\n")
    fout.write("RULE  1     14     1     0.0001      1    0.25\n")
    fout.write("RULE  1     15     1     0.0001      1    0.25\n")
    fout.write("RULE  1     16     1     0.0001      1    0.25\n")
    
    fout.write("RULE  1     13     3     0.0001      1    0.25\n")
    fout.write("RULE  1     14     3     0.0001      1    0.25\n")
    fout.write("RULE  1     15     3     0.0001      1    0.25\n")
    fout.write("RULE  1     16     3     0.0001      1    0.25\n")
    
    fout.write("RULE  1     13     5     0.0001      1    0.25\n")
    fout.write("RULE  1     14     5     0.0001      1    0.25\n")
    fout.write("RULE  1     15     5     0.0001      1    0.25\n")
    fout.write("RULE  1     16     5     0.0001      1    0.25\n")

    fout.write("RULE  1     13     7     0.0001      1    0.25\n")
    fout.write("RULE  1     14     7     0.0001      1    0.25\n")
    fout.write("RULE  1     15     7     0.0001      1    0.25\n")
    fout.write("RULE  1     16     7     0.0001      1    0.25\n")
    

    fout.close()



#################################
#
# MAIN part 1: run a normal SimGenome simulation. . .
#
commands = []

pe_scalars = [0.005]
count = 0
for pe in pe_scalars:
    count += 1
    thisid = ".bi" + count.__str__()
    OUTDIR = os.path.abspath( RUNID + thisid )
    if False == os.path.exists(OUTDIR):
        os.system( "mkdir " + OUTDIR )
    print_rules(outdir = OUTDIR)
    print_psams(outdir = OUTDIR)
    print_urss(outdir = OUTDIR)
    commands.append( get_command(id = thisid, popsize=40, randseed=123456789, pe_scalar=pe) )

# This vesion builds random populations
#commands.append( get_command(id = "rand1", random = True, popsize=48, randseed=123456789) )

fout = open(OUTDIR + "/" + RUNID + ".commands.txt", "w")
for c in commands:
    fout.write(c + "\n")
fout.close()
os.system("source " + OUTDIR + "/" + RUNID + ".commands.txt" )

#os.system("python tools/plot_expression.py " + RUNID + "a1 --skip_binding True")

