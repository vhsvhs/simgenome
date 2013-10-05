#
# USAGE:
# %> python tools/test_basic.py P
# . . . where P is the path to a desired output directory
#
#
#
# Screens values of mutation rates versus N generations to solution,
# for the simple problem of assembling a feed-forward loop.
#
# NOTE: This script requires that you install simreg as a system-wide link.
#
import sys, os

if sys.argv.__len__() < 2:
    print "USAGE:"
    print "  $> python test_basic.py OUTDIR"
    print ""
    exit()

OUTDIR = os.path.abspath( sys.argv[1] )

def get_command(random = False, popsave = False):
    command = "simreg "
    command += " --outdir " + OUTDIR + "/out.test_basic "
    command += "--maxgen 100 "
    command += "--niid 10000 "
    command += "--rulepath " + OUTDIR + "/test_basic.rules "
    command += "--popsize 20 "
    command += "--verbosity 5 "
    command += "--maxgd 1 "
    command += " --psam_mu 0.05 "
    command += " --urs_mu 0.005"
    command += " --f_scalar -2.0 "
    command += " --pe_scalar 0.001"
    command += " --elite_prop 0.29 "
    if popsave == False and random == False:
        command += "--psampath " + OUTDIR + "/test_basic.psam "
        command += "--urspath " + OUTDIR + "/test_basic.urs "
    if random:
        command += " --randompop"
    if popsave:
        command += " --poppath pop.gen0.save.backup"
    return command

def print_psams():
    fout = open(OUTDIR + "/test_basic.psam", "w")
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
    fout = open(OUTDIR + "/test_basic.urs", "w")
    fout.write(">0\n")
    fout.write("AAAAA\n")
    fout.write(">1\n")
    fout.write("AAAAAGGGGGAAAAAAAAACCCCAAAAAAAAACCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n")
    fout.write(">2\n")
    fout.write("AAAAAAAAAACCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n")
    fout.write(">3\n")
    fout.write("AAAA\n")
    fout.write(">4\n")
    fout.write("AAAA\n")
    fout.write(">5\n")
    fout.write("AAAACCCCCCAAAAAAAAAAATTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n")
   
    fout.close()

def print_rules():
    fout = open( OUTDIR + "/test_basic.rules", "w" )
    fout.write("RULE 0     5     4     1.0     0    1.0\n")
    fout.write("RULE 0     5     7     0.0001  1    1.0\n")
    fout.write("INPUT 0     0     0     7     0.5\n")
    fout.write("INPUT 0     1     0     7     0.5\n")
    #fout.write("INPUT 0 3 0 10 0.5\n")
    #fout.write("INPUT 0 4 0 10 0.4\n")
    fout.close()


#################################
#
# main
#

if False == os.path.exists(OUTDIR):
    os.system( "mkdir " + OUTDIR )
if os.path.exists(OUTDIR + "/out.test_basic"):
    os.system("rm " + OUTDIR + "/out.test_basic")
os.system("mkdir " + OUTDIR + "/out.test_basic")
print_rules()
print_psams()
print_urss()
commands = []
#commands.append( get_command() )
#commands.append( "cp " + OUTDIR + "/out.test_basic/POPS/pop.gen0.save.txt" + " pop.gen0.save.backup")
#commands.append( get_command(popsave = True) )
#commands.append( get_command(popsave = True) )
#commands.append( get_command(popsave = True) )
commands.append( get_command(random = True) )

fout = open(OUTDIR + "/test_basic.commands.txt", "w")
for c in commands:
    fout.write(c + "\n")
fout.close()
os.system("source " + OUTDIR + "/test_basic.commands.txt" )

