#
# Screens values of mutation rates versus N generations to solution,
# for the simple problem of assembling a feed-forward loop.
#
# NOTE: This script requires that you install simreg as a system-wide link.
#
import sys, os

OUTDIR = os.path.abspath( sys.argv[1] )

def get_command():
    command = "simreg "
    command += " --outdir " + OUTDIR + "/out.test_basic "
    command += "--psampath " + OUTDIR + "/test_basic.psam "
    command += "--urspath " + OUTDIR + "/test_basic.urs "
    command += "--rulepath " + OUTDIR + "/test_basic.rules "
    command += "--maxgen 1 "
    command += "--niid 10000 "
    command += "--popsize 1 "
    command += "--verbosity 10 "
    command += "--maxgd 1 "
    command += " "
    return command

def print_psams():
    fout = open(OUTDIR + "/test_basic.psam", "w")
    # Gene 0 likes CCC
    fout.write("Gene 0 1\n")
    fout.write("0.0    30.0    0.0    0.0\n")
    fout.write("0.0    30.0    0.0    0.0\n")
    fout.write("0.1    30.0    0.0    0.0\n")
    # Gene 1 likes GGG
    fout.write("Gene 1 0\n")
    fout.write("0.0    0.0    30.0    0.0\n")
    fout.write("0.0    0.0    30.0    0.0\n")
    fout.write("0.1    0.0    30.0    0.0\n")
    fout.close()
    
def print_urss():
    fout = open(OUTDIR + "/test_basic.urs", "w")
    fout.write(">0\n")
    fout.write("AAAAAAAAAAAAAAAAAAAAAAAAAAA\n")
    fout.write(">1\n")
    fout.write("AAAAAAAAAAAAAAAAAAAAAAAAAAA\n")
    fout.write(">2\n")
    fout.write("AAAAAAAAAAAAAAAAAAAAAAAAAAA\n")
    fout.close()

def print_rules():
    fout = open( OUTDIR + "/test_basic.rules", "w" )
    fout.write("RULE 0 2 4 0.5 0 1\n")
    fout.write("RULE 0 2 8 0.01 1 1\n")

    # timepatternID, basal_gene_id, time_start, time_stop, expr_level
    fout.write("INPUT 0 0 0 8 1\n")
    fout.write("INPUT 0 1 0 8 1\n")
    fout.close()


#################################
#
# main
#

if os.path.exists(OUTDIR + "/out.test_basic"):
    os.system("rm " + OUTDIR + "/out.test_basic")
if False == os.path.exists(OUTDIR):
    os.system( "mkdir " + OUTDIR )
os.system("mkdir " + OUTDIR + "/out.test_basic")
print_rules()
print_psams()
print_urss()
commands = []
c = get_command()
commands.append( c )

fout = open(OUTDIR + "/test_basic.commands.txt", "w")
for c in commands:
    fout.write(c + "\n")
fout.close()
