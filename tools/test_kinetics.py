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
    fout.write("0.0 0.0 4.0 0\n")
    fout.write("0.0 0.0 4.0 0\n")
    fout.write("0.0 0.0 4.0 0\n")
    fout.write("0.0 0.0 4.0 0\n")
    fout.write("0.0 0.0 4.0 0\n")
    fout.write("0.0 0.0 4.0 0\n")

    fout.write("# non-spec activation\n")
    fout.write("Gene 1 1\n")
    fout.write("1 1 1 1\n")
    fout.write("1 1 1 1\n")
    fout.write("1 1 1 1\n")
    fout.write("1 1 1 1\n")

    fout.write("# repressor binds TTTT\n")
    fout.write("Gene 2 0\n")
    fout.write("0.0 0.0 0.0 3\n")
    fout.write("0.0 0.0 0.0 3\n")
    fout.write("0.0 0.0 0.0 3\n")
    fout.write("0.0 0.0 0.0 3\n")
    fout.write("0.0 0.0 0.0 3\n")
#     fout.write("0.0 0.0 0.0 3\n")
#     fout.write("0.0 0.0 0.0 3\n")
#     fout.write("0.0 0.0 0.0 3\n")
#     fout.write("0.0 0.0 0.0 3\n")

    
def print_urss( outdir=""):
    fout = open(outdir + "/" + RUNID + ".urs", "w")
    fout.write(">0\n")
    fout.write( "AAAAA\n")
    fout.write(">1\n")
    fout.write( "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n")
    fout.write(">2\n")
    fout.write( "AAAAAAAAAGGGGGGGGGAAAAAAAAAAAAAAA\n")
    fout.write(">3\n")
    fout.write( "AAAAAAAAATAAAAAAAAAAAAAAAAAAAAAAA\n")
    fout.write(">4\n")
    fout.write( "AAAAAAAAATTAAAAAAAAAAAAAAAAAAAAAA\n")
    fout.write(">5\n")
    fout.write( "AAAAAAAAATTTAAAAAAAAAAAAAAAAAAAAA\n")
    fout.write(">6\n")
    fout.write( "AAAAAAAAATTTTAAAAAAAAAAAAAAAAAAAA\n")
    fout.write(">7\n")
    fout.write( "AAAAAAAAATTTTTAAAAAAAAAAAAAAAAAAA\n")
    fout.write(">8\n")
    fout.write( "AAAAAAAAATTTTTTAAAAAAAAAAAAAAAAAA\n")
    fout.write(">9\n")
    fout.write( "AAAAAAAAATTTTTTTAAAAAAAAAAAAAAAAA\n")
    fout.write(">10\n")
    fout.write( "AAAAAAAAATTTTTTTTAAAAAAAAAAAAAAAA\n")
    fout.write(">11\n")
    fout.write( "AAAAAAAAATTTTTTTTTAAAAAAAAAAAAAAA\n")
    fout.write(">12\n")
    fout.write( "AAAAAAAAATTTTTTTTTTAAAAAAAAAAAAAA\n")
    fout.write(">13\n")
    fout.write( "AAAAAAAAATTTTTTTTTTTAAAAAAAAAAAAA\n")
    fout.write(">14\n")
    fout.write( "AAAAAAAAATTTTTTTTTTTTAAAAAAAAAAAA\n")
    fout.write(">15\n")
    fout.write( "AAAAAAAAATTTTTTTTTTTTTAAAAAAAAAAA\n")
    fout.write(">16\n")
    fout.write( "AAAAAAAAATTTTTTTTTTTTTTAAAAAAAAAA\n")
    fout.write(">17\n")
    fout.write( "AAAAAAAAATTTTTTTTTTTTTTTAAAAAAAAA\n")
    fout.write(">18\n")
    fout.write( "AAAAAAAAATTTTTTTTTTTTTTTTAAAAAAAA\n")
    fout.write(">19\n")
    fout.write( "AAAAAAAAATTTTTTTTTTTTTTTTTAAAAAAA\n")
    fout.close()

def print_rules( outdir = ""):
    fout = open( outdir + "/" + RUNID + ".rules", "w" )
    fout.write("# inducer pulse\n")
    fout.write("INPUT 0     0     0     9      0.00001\n")
    fout.write("INPUT 0     0     10     12     0.0001\n")
    fout.write("INPUT 0     0     12     14     0.001\n")
    fout.write("INPUT 0     0     14     16     0.01\n")
    fout.write("INPUT 0     0     16     25     0.1\n")
    fout.write("INPUT 0     0     25     27     0.01\n")
    fout.write("INPUT 0     0     27     29     0.001\n")
    fout.write("INPUT 0     0     29     31     0.0001\n")
    fout.write("INPUT 0     0     31     50     0.00001\n")
    
    fout.write("# low-level activation always on\n")
    fout.write("INPUT 0     1     0     50     0.0001\n")


    #fout.write("RULE 0     5     4     1.0     0    1.0\n")
    #fout.write("RULE 0     5     7     0.0001  1    1.0\n")
    fout.close()

# analyze for Kd
def analyze_kd(outdir = ""):
    """This function gets data for a kinetics experiment, December 2013.
    It contains some ad-hoc things, and should be considered depricated
    for future/general use."""
    print "\n================================================="
    print ". Analyzing for Kd..."



    t_g_o = {} # key = time, value = hashtable, where key = gene ID, value = occupancy by regulator gene 2 on this gene
    t_g_e = {} # key = time, value = hashtable, where key = gene ID, value = expression

    # dec2:
    timepoints = [31,32,34,35,36]
    for t in timepoints:
        t_g_o[t] = {}
        t_g_e[t] = {}
    
    found_genes = []
    print ". Occupancy of Regulator 2, and Expression of Gene X, at timepoints:", timepoints
    fin = open(outdir + "/LOGS/expression.txt", "r")
    for l in fin.xreadlines():
        if l.__len__() > 10:
            #print l
            tokens = l.split()
            time = int(tokens[5])
            #print time, timepoints
            if time in timepoints:
                gene = int(tokens[9])
                if gene not in found_genes:
                    found_genes.append( gene )
                expr = 0
                if tokens[11].startswith("expr"):
                    expr = float(tokens[12])
                else:
                    expr = float(tokens[11])
                t_g_e[time][gene] = expr
                #print time, gene, expr
    fin.close()
    
    for gene in found_genes:
        fin = open(outdir + "/OCCUPANCY/occ.gen0.id0.gene" + gene.__str__() + ".txt", "r")
        foundit = False
        for l in fin.xreadlines():
            if l.__len__() > 10:
                for t in timepoints:
                    if l.__contains__("time " + t.__str__()):
                        foundit = True
                    elif l.__contains__("time") and False == l.__contains__("time " + t.__str__()):
                        foundit = False
                    elif foundit:
                        tokens = l.split()
                        #print tokens
                        o = float(tokens[11])
                        if gene not in t_g_o[t]:
                            t_g_o[t][gene] = o
                        elif t_g_o[t][gene] < o:
                            t_g_o[t][gene] = o
                        """
                        if gene in t_g_o[t]:
                            if t_g_o[t][gene] < o:
                              t_g_o[t][gene] = o
                        else:    
                            t_g_o[t][gene] = o
                        """
        fin.close()
    
    print ". Occupancy P(gene2)..."
    out = "Gene\t"
    for t in timepoints:
        out += "O" + t.__str__() + "\tE" + t.__str__() + "\t|\t"
    print out
    for gene in found_genes:
        out = gene.__str__() + "\t"
        for t in timepoints:
            out += t_g_o[t][gene].__str__() + "\t"
            out += t_g_e[t][gene].__str__() + "\t"
            out += "|\t"
        print out
    
    
    maxx = None
    maxy = None
    fout = open(outdir + "/PLOTS/occ-exp.rscript", "w")
    for t in timepoints:
        xs = "x" + t.__str__() + " <-c("
        ys = "y" + t.__str__() + " <-c("
        for gene in found_genes:
            x = t_g_o[t][gene]
            y = t_g_e[t][gene]
            if x > maxx:
                maxx = x
            if y > maxy:
                maxy = y
            xs += x.__str__() + ","
            ys += y.__str__() + ","
        xs = re.sub(",$", "", xs)
        xs += ");\n"
        ys = re.sub(",$", "", ys)
        ys += ");\n"
        fout.write(xs + "\n")
        fout.write(ys + "\n")
    fout.write("pdf('" + outdir + "/PLOTS/occ-exp.pdf');\n")
    fout.write("plot(c(0.0," + maxx.__str__() + "), c(0," + maxy.__str__() + "), type='n', main='', xlab='Occupancy of Repressor', ylab='Expression');\n")
    for t in timepoints:
        fout.write("points(x" + t.__str__() + ",y" + t.__str__() + ",type='l');\n" )
    fout.write("dev.off();\n")
    os.system("r --no-save < " + outdir + "/PLOTS/occ-exp.rscript")

#################################
#
# MAIN part 1: run a normal SimGenome simulation. . .
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

#
# MAIN part 2: post-analysis of transcription DNA-binding kinetics
analyze_kd(outdir = OUTDIR)

os.system("python tools/plot_expression.py " + RUNID + "a1 --skip_binding True")

