import os
import random
import re
from random import choice
import scipy.stats as ss
import math

from pwm_library import *

def get_random_nt(n):
    """Returns a random nucleotide sequence n sites long."""
    line = ""
    for i in range(0, n):
        line += choice(["A", "C", "G", "T"])
    return line


def run_segal2008():    
    baseid = "segal2008"
    os.system("mkdir UTESTS")    
        
    write_files( "./UTESTS/pwm." + baseid + ".txt" , "./UTESTS/urs." + baseid + ".fasta" )
     
    fout = open("./UTESTS/rules." + baseid +".txt", "w")
    fout.write("# timepatternID, reporter_gene_id, timepoint, expression_level, rule_type, multiplier\n")
    fout.write("RULE 0 6 1 1.0 le\n")

    
    fout.write("# timepatternID, basal_gene_id, time_start, time_stop, expr_level\n")
    fout.write("INPUT 0 0 0 1 0.00\n") #hunchback
    fout.write("INPUT 0 1 0 1 0.00\n") #tailless
    fout.write("INPUT 0 2 0 1 0.00\n") #giant
    fout.write("INPUT 0 3 0 1 0.00\n") #kruppel
    fout.write("INPUT 0 4 0 1 0.350\n") #knirps
    fout.write("INPUT 0 5 0 1 0.875\n") #bicoid
    fout.write("INPUT 0 6 0 1 0.00\n") #caudal
    fout.write("INPUT 0 7 0 1 0.946\n") #torRE
    fout.close()
    
    
    #print "utests.py 456", n, rep
    cismu = 0.00
    dbdmu = 0.0
    runid = baseid
    scriptname = "UTESTS/" + baseid + ".run.sh"
    fout = open(scriptname, "w")
    line = "mpirun -np 2 --machinefile hosts.txt python runme.py "
    line += " --numtr 8 --numreporter 44 --popsize 1 "
    line += " --patternpath ./UTESTS/rules." + baseid + ".txt "
    line += " --urspath ./UTESTS/urs." + baseid + ".fasta "
    line += " --verbose 91 "
    line += " --runid " + runid 
    line += " --growth_rate 1.0 --decay_rate 1.0 "
    line += " --pe_scalar 0.001 " 
    line += " --pwmpath ./UTESTS/pwm." + baseid + ".txt "
    line += " --workspace ./UTESTS --tfcoop zeros "
    line += " --maxgd 1 "
    line += " --iid_samples 3000" 
    line += " --eliteproportion 0.0 "            
    line += " --elitemu 0.0 "
    line += " --mu 0 "
    line += " --cismu " + cismu.__str__() # VARIABLE
    line += " --dbdmu " + dbdmu.__str__() # VARIABLE
    line += " --pwmdeltamax 0.5 "
    line += " --pwmlenmu 0.0 --pwmmulenmax 2"
    line += " --urslenmu 0.0 "
    line += " --p2pmu 0.0"
    line += " --p2pmudelta 1.0"
    line += " --tfcoop zeros " # disable CO-op interactions
    line += " --sexual_ratio 0.8"
    line += " --maxgd 1"
    line += " --stop_early 0.8" # stop when the mean pop. f reaches 0.8
    line += " --maxgenerations 1 --maxtime 1 "

    fout.write(line + "\n")
    fout.close()
    os.system("source " + scriptname)




def write_files(pwmpath, urspath):
    fin = open("segal2008/tf_weight_matrices.gxw.txt", "r")
    lines = fin.readlines()
    fin.close()
    
    tf_pwmlines = {}
    lasttf = None
    for l in lines:
        if l.startswith("<WeightMatrix"):
            tokens = l.split()
            tf = re.sub("Name=\"", "", tokens[1])
            tf = re.sub("\"", "", tf)
            lasttf = tf
            tf_pwmlines[tf] = ""
        elif l.__contains__("<Position Weights"):
            tokens = l.split()
            wts = re.sub("Weights=\"", "", tokens[1])
            wts = re.sub("\"></Position>", "", wts)
            wts.strip()
            wts = re.sub(";", "  ", wts)
            tf_pwmlines[lasttf] += wts + "\n"
    
    fin = open("segal2008/module_sequence.tab.txt", "r")
    lines = fin.readlines()
    fin.close()
    
    gene_urs = {}
    lastgene = None
    for l in lines:
        if l.startswith(">"):
            l = l.strip()
            gene = re.sub(">", "", l)
            lastgene = gene
        elif lastgene != None:
            l = l.strip()
            gene_urs[lastgene] = l

    counttf = 0
    fout = open(pwmpath, "w")
    for tf in tf_pwmlines:
        #print tf, tf_pwmlines[tf]
        fout.write("# " + tf + "\n")
        fout.write("pwm " + counttf.__str__() + "  0\n")
        fout.write(tf_pwmlines[tf] + "\n")
        counttf += 1
    fout.close()
    
    countgene = 0
    fout = open(urspath, "w")
    for tf in tf_pwmlines:
        fout.write("# " + tf + "\n")
        fout.write(">" + countgene.__str__() + " " + tf + "\n")
        fout.write("AAAAAAAAAA\n")
        countgene += 1
    for gene in gene_urs:
        fout.write("# " + gene + "\n")
        fout.write(">" + countgene.__str__() + " " + gene + "\n")
        fout.write(gene_urs[gene] + "\n")
        countgene += 1
    fout.close()
    #exit()

run_segal2008()

    
    
    