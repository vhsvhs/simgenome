#
# Screens values of mutation rates versus N generations to solution,
# for the simple problem of assembling a feed-forward loop.
#
# NOTE: This script requires that you install simreg as a system-wide link.
#
import sys, os

OUTDIR = sys.argv[1]


def make_runid(psammu, ursmu):
    return "p" + psammu.__str__() + "--u" + ursmu.__str__()

def get_command(psammu, ursmu):
    command = "simreg "
    command += " --outdir out." + make_runid(psammu, ursmu) + "/ " 
    command += "--psampath test_mu.psam "
    command += "--urspath test_mu.urs "
    command += "--rulepath test_my.rules "
    command += "--maxgen 5000 "
    command += "--niid 5000 "
    command += "--popsize 1 "
    command += "--verbosity 10 "
    command += "--run_clean "
    command += "--maxgd 3 "

def print_psams():
    fout = open(OUTDIR + "/test_mu.psam", "w")
    fout.write("Gene 0 1\n")
    fout.write("1.0 0.0 0.0 -1.0\n")
    fout.write("3.0 0.0 0.0 -1.0\n")
    fout.write("3.1 -0.1 -0.1 4.0\n")
    fout.write("Gene 1 0\n")
    fout.write("1.0 30.0 0.0 0.0\n")
    fout.write("1 30.0 0 0\n")
    fout.write("1 30 0.0 0.1\n")
    fout.close()
    
def print_urss():
    fout = open(OUTDIR + "/test_mu.urs", "w")
    ax1000 = ""
    for i in range(0, 1000):
        ax1000 += "A"
    fout.write(">0\n")
    fout.write(ax1000 + "\n")
    fout.write(">1\n")
    fout.write(ax1000 + "\n")
    fout.write(">2\n")
    fout.write(ax1000 + "\n")
    fout.close()

