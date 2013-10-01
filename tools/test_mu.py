#
# Screens values of mutation rates versus N generations to solution,
# for the simple problem of assembling a feed-forward loop.
#
# NOTE: This script requires that you install simreg as a system-wide link.
#
import sys, os

OUTDIR = sys.argv[1]
#OUTDIR = os.path.abspath(OUTDIR)
PSAM = [0.0]
URS = [0.001,0.005,0.01,0.05,0.1]
#URS = [0.05]
nreps = 4

nprocs = PSAM.__len__() * URS.__len__() * nreps # how many MPI processes?

def make_runid(psammu, ursmu, rep):
    return "p" + psammu.__str__() + "-u" + ursmu.__str__() + "-r" + rep.__str__()

def make_directory(psammu, ursmu, rep):
    runid = make_runid(psammu, ursmu, rep)
    print OUTDIR + "/out." + runid
    if False == os.path.exists(OUTDIR + "/out." + runid):
        os.system("mkdir " + OUTDIR + "/out." + runid)

def get_command(psammu, ursmu, rep):
    runid = make_runid(psammu, ursmu, rep)
    command = "/common/bin/simreg "
    command += " --outdir " + OUTDIR + "/out." + runid + " " 
    command += "--psampath " + OUTDIR + "/test_mu.psam "
    command += "--urspath " + OUTDIR + "/test_mu.urs "
    command += "--rulepath " + OUTDIR + "/test_mu.rules "
    command += "--maxgen 500 "
    command += "--niid 10000 "
    command += "--popsize 10 "
    command += "--verbosity 10 "
    command += "--run_clean "
    command += "--maxgd 1 "
    command += " --urs_mu " + ursmu.__str__()
    command += " --psam_mu " + psammu.__str__()
    command += " --pe_scalar 0.001"
    command += " --growth_scalar 0.3"
    command += " --decay_scalar 0.3"
    return command

def print_psams():
    fout = open(OUTDIR + "/test_mu.psam", "w")
    fout.write("# activator, GGG\n")
    fout.write("Gene 0 1\n")
    fout.write("0.0 0.1 30.0 0.1\n")
    fout.write("0.0 0.1 30.0 0.1\n")
    fout.write("0.0 0.1 30.0 0.1\n")
    fout.write("0.0 0.1 30.0 0.1\n")
    fout.write("# repressor, TTTT\n")
    fout.write("Gene 1 0\n")
    fout.write("0.0 0.1 0.1 30.0\n")
    fout.write("0.0 0.1 0.1 30.0\n")
    fout.write("0.0 0.1 0.1 30.0\n")
    fout.write("0.0 0.1 0.1 30.0\n")
    fout.write("0.0 0.1 0.1 30.0\n")
    fout.write("# non-spec binding\n")
    fout.write("Gene 2 0\n")
    fout.write("1 1 1 1\n")
    fout.write("1 1 1 1\n")
    fout.write("# non-spec binding\n")
    fout.write("Gene 3 1\n")
    fout.write("1 1 1 1\n")
    fout.write("1 1 1 1\n")
    fout.close()
    
def print_urss():
    fout = open(OUTDIR + "/test_mu.urs", "w")
    fout.write(">0\n")
    fout.write("AAAAA\n")
    fout.write(">1\n")
    fout.write("AAAAAGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n")
    fout.write(">2\n")
    fout.write("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n")
    fout.write(">3\n")
    fout.write("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n")
    fout.write(">4\n")
    fout.write("AAAAGGGGAAAAAAAAAAATTTTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n")
    fout.close()

def print_rules():
    fout = open( OUTDIR + "/test_mu.rules", "w" )
    fout.write("RULE 0 4 4 0.1 0 1.0\n")
    fout.write("RULE 0 4 4 0.0001 1 1.0\n")
    fout.write("INPUT 0 0 0 10 0.6\n")
    fout.write("INPUT 0 2 0 10 0.5\n")
    fout.write("INPUT 0 3 0 10 0.4\n")
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


fout = open("test_mu.commands.txt", "w")
for c in commands:
    fout.write(c + "\n")
fout.close()

print "I wrote a BATCH script to 'test_mu.commands.txt'"

print "mpirun -np " + nprocs.__str__() + " /common/bin/mpi_dispatch test_mu.commands.txt"
