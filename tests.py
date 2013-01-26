import os
import random
from random import choice
import scipy.stats as ss
import math

def gen_nonspec_nt(n):
    count = 0
    line = ""
    while(True):
        for a in ["A", "C", "G", "T"]:
            line += a
            count += 1
            if count == n:
                return line

def get_random_nt(n):
    line = ""
    for i in range(0, n):
        line += choice(["A", "C", "G", "T"])
    return line

def get_random_pwm(n):
    line = ""
    for i in range(0, n):
        vals = []
        sum = 0.0
        for i in range(0,4):
            x = random.random()
            vals.append( x )
            sum += x
        for val in vals:
            line += (val/sum).__str__() + " "
        line += "\n"
    return line

def get_nonspec_pwm(n):
    line = ""
    for i in range(0, n):
        line += "0.25 0.25 0.25 0.25\n"
    return line

def get_Aspec(n):
    line = ""
    for i in range(0, n):
        line += "1.00 0.0 0.0 0.0\n"
    return line

def get_Cspec(n):
    line = ""
    for i in range(0, n):
        line += "0.00 1.0 0.0 0.0\n"
    return line

def get_Gspec(n):
    line = ""
    for i in range(0, n):
        line += "0.00 0.0 1.0 0.0\n"
    return line

def get_Tspec(n):
    line = ""
    for i in range(0, n):
        line += "0.00 0.0 0.0 1.0\n"
    return line

def get_matalpha2():
    # ATTTACATG
    line = ""
    line += "0.97 0.01 0.01 0.01\n"
    line += "0.60 0.01 0.01 0.38\n"
    line += "0.10 0.01 0.09 0.80\n"
    line += "0.01 0.01 0.01 0.97\n"
    line += "0.97 0.01 0.01 0.01\n"
    line += "0.01 0.97 0.01 0.01\n"
    line += "0.97 0.01 0.01 0.01\n"
    line += "0.01 0.20 0.01 0.78\n"
    line += "0.08 0.01 0.90 0.01\n"    
    return line

def get_matalpha2scram():
    # ATAATTCGA
    line = ""
    line += "0.97 0.01 0.01 0.01\n"
    line += "0.10 0.01 0.09 0.80\n"
    line += "0.97 0.01 0.01 0.01\n"
    line += "0.97 0.01 0.01 0.01\n"
    line += "0.01 0.01 0.01 0.97\n"
    line += "0.01 0.20 0.01 0.78\n"
    line += "0.01 0.97 0.01 0.01\n"
    line += "0.08 0.01 0.90 0.01\n"    
    line += "0.60 0.01 0.01 0.38\n"
    return line

def get_ste12():
    line = ""
    line += "0.01 0.01 0.01 0.97\n"
    line += "0.01 0.01 0.97 0.01\n"
    line += "0.97 0.01 0.01 0.01\n"
    line += "0.97 0.01 0.01 0.01\n"
    line += "0.97 0.01 0.01 0.01\n"
    line += "0.01 0.97 0.01 0.01\n"
    line += "0.70 0.01 0.28 0.01\n"  
    return line

def get_efg1():
    line = ""
    line += "0.036599    0.270976    0.002458    0.689964\n"
    line += "0.817856    0.177225    0.002458    0.002458\n"
    line += "0.020974    0.005349    0.002458    0.971217\n"
    line += "0.020974    0.036599    0.939967    0.002458\n"
    line += "0.005349    0.942857    0.018083    0.033709\n"
    line += "0.974108    0.005349    0.002458    0.018083\n"
    line += "0.13035     0.427228    0.002458    0.439962\n"
    line += "0.349102    0.13035     0.471212    0.049334\n"
    return line

def get_efg1_seq():
    return ""

# From Zhu
def get_tec1():
    line = ""
    line += "0.160117    0.430781    0.198437    0.210664\n"
    line += "0.500104    0.218635    0.073431    0.207828\n"
    line += "0.753885    0.00989    0.225484    0.01074\n"
    line += "0.079572    0.917027    0.002496    0.000904\n"
    line += "0.963133    0.000654    0.034699    0.001513\n"
    line += "0.013733    0.000967    0.000818    0.98448\n"
    line += "0.049433    0.00043    0.001742    0.948393\n"
    line += "0.000344    0.980815    0.007453    0.011387\n"
    line += "0.001567    0.482371    0.008104    0.507956\n"
    line += "0.007707    0.399204    0.01281    0.580277\n"
    line += "0.14392    0.163985    0.563246    0.128847\n"
    return line

# From Zhu
def get_ndt80():
    line = ""
    line += "0.215888    0.403718    0.140864    0.239529\n"
    line += "0.201577    0.430686    0.17625    0.191484\n"
    line += "0.0994    0.114555    0.495278    0.290765\n"
    line += "0.188009    0.020196    0.675723    0.116069\n"
    line += "0.444225    0.483732    0.041969    0.030072\n"
    line += "0.012211    0.948956    0.004543    0.034287\n"
    line += "0.940416    0.007801    0.04533    0.006451\n"
    line += "0.018448    0.967691    0.00696    0.0069\n"
    line += "0.965515    0.004107    0.009598    0.020778\n"
    line += "0.97028    0.005013    0.01588    0.008826\n"
    line += "0.945264    0.003527    0.014206    0.037001\n"
    line += "0.713859    0.024864    0.236744    0.024531\n"
    line += "0.496837    0.178778    0.082445    0.241937\n"
    line += "0.123789    0.630913    0.09573    0.149566\n"
    line += "0.127851    0.32957    0.440099    0.102478\n"
    line += "0.314093    0.330225    0.277595    0.078086\n"
    line += "0.360795    0.177409    0.230168    0.231626\n"
    line += "0.477799    0.192623    0.174716    0.15486\n"
    return line



#
# A simple FFL
#
def run_ffl():
    #runid = "ffl.12-14-2012"
    runid = "ffl.12-14-2012ko"
    scriptname = "UTESTS/" + runid + ".run.sh"
    os.system("mkdir UTESTS")
    fout = open(scriptname, "w")

    # KO:
    fout.write("mpirun -np 25 --machinefile hosts.txt python runme.py --ko 2 --start_generation 100 --pop_path ./UTESTS/ffl.12-14-2012/POP_HISTORY/population.gen82.pickle --numtr 4 --numreporter 2 --popsize 48 --patternpath ./UTESTS/rules." + runid + ".txt --urspath ./UTESTS/urs." + runid + ".txt --verbose 8 --runid " + runid + " --growth_rate 0.5 --decay_rate 0.5 --pwmpath ./UTESTS/pwm." + runid + ".txt --maxgenerations 1000 --maxtime 1 --workspace ./UTESTS --tfcoop zeros --maxgd 1 --iid_samples 10000\n")

    # d - same as c, but pe_scaler changed
    #fout.write("mpirun -np 25 --machinefile hosts.txt python runme.py --numtr 4 --numreporter 2 --popsize 48 --patternpath ./UTESTS/rules." + runid + ".txt --urspath ./UTESTS/urs." + runid + ".txt --verbose 8 --runid " + runid + " --growth_rate 0.5 --decay_rate 0.5 --pwmpath ./UTESTS/pwm." + runid + ".txt --maxgenerations 1000 --maxtime 1 --workspace ./UTESTS --tfcoop zeros --maxgd 1 --iid_samples 10000 --mu 1 --eliteproportion 0.3 --elitemu 0.0 --p2pmu 0.0 --dbdmu 0.2 --pwmmu 1.0 --pwmlenmu 0.5 --pwmmulenmax 2 --urslenmu 0.1 --cismu 0.5 --sexual_ratio 0.7 --pe_scalar 0.5\n")

    fout.close()
    
    fout = open("./UTESTS/rules." + runid + ".txt", "w")
    fout.write("# timepatternID, basal_gene_id, timepoint, expression_level, reporter_gene_ID, rule_type, multiplier\n")
    fout.write("RULE 0 0 4 0.01 5 le\n")
    fout.write("RULE 0 0 8 0.01 5 le\n")
    fout.write("RULE 0 0 4 0.4 4 ge 4\n")
    fout.write("RULE 0 0 8 0.01 4 le\n")
    
    fout.write("# timepatternID, basal_gene_id, time_start, time_stop, expr_level\n")
    fout.write("INPUT 0 0 0 8 1.00\n")
    fout.close()
    
    fout = open("./UTESTS/urs." + runid + ".txt", "w")
    fout.write("# ID, has_dbd (0 or 1), is_repressor (0 or 1), URS\n")
    fout.write("0 1 0 " + get_random_nt(30) + "\n")
    fout.write("1 1 0 " + get_random_nt(30) + "\n")
    fout.write("2 1 1 " + get_random_nt(30) + "\n")
    fout.write("3 1 1 " + get_random_nt(30) + "\n")
    fout.write("4 0 0 " + get_random_nt(30) + "\n")
    fout.write("5 0 0 " + get_random_nt(30) + "\n")
    
    #fout.write("0 1 0 AAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n")
    #fout.write("1 1 1 GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n")
    #fout.write("2 0 0 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n")
    
    fout.close()
    
    fout = open("./UTESTS/pwm." + runid + ".txt", "w")
    fout.write("pwm 0\n")
    #fout.write( get_nonspec_pwm(5) )
    fout.write( get_random_pwm(6) )
    fout.write("pwm 1\n")
    #fout.write( get_nonspec_pwm(5) )   
    fout.write( get_random_pwm(6) )
    fout.write("pwm 2\n")
    #fout.write( get_nonspec_pwm(5) )   
    fout.write( get_random_pwm(6) )
    fout.write("pwm 3\n")
    #fout.write( get_nonspec_pwm(5) )   
    fout.write( get_random_pwm(6) )
    fout.close()
    
    os.system("source " + scriptname)
    exit()

#
# A simple FFL
#
def run_ko():
    os.system("mkdir UTESTS")
    scriptpath = "UTESTS/ko.12-14-2012.run.sh"
    fout = open(scriptpath, "w")
    
    #fout.write("mpirun -np 2 --machinefile hosts.txt python runme.py --popsize 4  --patternpath ./UTESTS/ffl_rules.txt --urspath ./UTESTS/ffl_urs.txt --verbose 8 --runid ko.12-14-2012 --growth_rate 0.5 --decay_rate 0.5 --pwmpath ./UTESTS/ffl_pwm.txt --maxgenerations 200 --maxtime 1 --workspace ./UTESTS --tfcoop zeros --maxgd 1 --iid_samples 500 --mu 0.1 --elitemu 0.0 --p2pmu 0.0 --dbdmu 0.5 --pwmmu 1.0 --urslenmu 0.0 --cismu 0.0 --sexual_ratio 0.5\n")
    fout.write("mpirun -np 2 --machinefile hosts.txt python runme.py --popsize 4 --ko 1 0 --patternpath ./UTESTS/ffl_rules.txt --urspath ./UTESTS/ffl_urs.txt --verbose 8 --runid ko.12-14-2012 --growth_rate 0.5 --decay_rate 0.5 --pwmpath ./UTESTS/ffl_pwm.txt --maxgenerations 200 --maxtime 1 --workspace ./UTESTS --tfcoop zeros --maxgd 1 --iid_samples 1000\n")

    fout.close()
    
    fout = open("./UTESTS/ffl_rules.txt", "w")
    fout.write("# timepatternID, basal_gene_id, timepoint, expression_level, reporter_gene_ID, rule_type\n")
    fout.write("RULE 0 0 2 0.3 2 ge\n")
    fout.write("RULE 0 0 4 0.01 2 le\n")
    
    fout.write("# timepatternID, basal_gene_id, time_start, time_stop, expr_level\n")
    fout.write("INPUT 0 0 0 1 1.00\n")
    fout.close()
    
    fout = open("./UTESTS/ffl_urs.txt", "w")
    fout.write("# ID, has_dbd (0 or 1), is_repressor (0 or 1), URS\n")    
    fout.write("0 1 0 AAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n")
    fout.write("1 1 1 GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n")
    fout.write("2 0 0 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n")
    fout.close()
    
    fout = open("./UTESTS/ffl_pwm.txt", "w")
    fout.write("pwm 0\n")
    fout.write( get_Cspec(5) )
    fout.write("pwm 1\n") 
    fout.write( get_Tspec(5) )
    fout.close()
    
    os.system("source " + scriptpath)
    exit()


def run_tune_mu():
    #baseid = "mu.12-20-2012noelt"
    baseid = "mu.01-01-2013noelt"

    os.system("mkdir UTESTS")
    
    fout = open("./UTESTS/rules." + baseid + ".txt", "w")
    fout.write("# timepatternID, basal_gene_id, timepoint, expression_level, reporter_gene_ID, rule_type, multiplier\n")
    fout.write("RULE 0 0 4 0.01 5 le\n")
    fout.write("RULE 0 0 8 0.01 5 le\n")
    fout.write("RULE 0 0 4 0.4 4 ge 4\n")
    fout.write("RULE 0 0 8 0.01 4 le\n")
    
    fout.write("# timepatternID, basal_gene_id, time_start, time_stop, expr_level\n")
    fout.write("INPUT 0 0 0 8 1.00\n")
    fout.close()
    
    fout = open("./UTESTS/urs." + baseid + ".txt", "w")
    fout.write("# ID, has_dbd (0 or 1), is_repressor (0 or 1), URS\n")
    fout.write("0 1 0 " + get_random_nt(100) + "\n")
    fout.write("1 1 0 " + get_random_nt(100) + "\n")
    fout.write("2 1 1 " + get_random_nt(100) + "\n")
    fout.write("3 1 1 " + get_random_nt(100) + "\n")
    fout.write("4 0 0 " + get_random_nt(100) + "\n")
    fout.write("5 0 0 " + get_random_nt(100) + "\n")
    fout.close()
    
    fout = open("./UTESTS/pwm." + baseid + ".txt", "w")
    fout.write("pwm 0\n")
    #fout.write( get_nonspec_pwm(5) )
    fout.write( get_random_pwm(6) )
    fout.write("pwm 1\n")
    #fout.write( get_nonspec_pwm(5) )   
    fout.write( get_random_pwm(6) )
    fout.write("pwm 2\n")
    #fout.write( get_nonspec_pwm(5) )   
    fout.write( get_random_pwm(6) )
    fout.write("pwm 3\n")
    #fout.write( get_nonspec_pwm(5) )   
    fout.write( get_random_pwm(6) )
    fout.close()
    
    scriptname = "UTESTS/" + baseid + ".run.sh"
    
    cismus = [0.0, 0.001, 0.01, 0.1]
    dbdmus = [0.0, 0.01, 0.1, 0.2, 0.5, 1.0]
    #cismus = [1.0]
    #dbdmus = [1.5,2.0]
    #dbdmus = [3.0, 4.0]
    for cismu in cismus:
        for dbdmu in dbdmus:
            if cismu == 0.0 and dbdmu == 0.0:
                continue
            runid = baseid + "cis" + cismu.__str__() + ".dbd" + dbdmu.__str__()
            fout = open(scriptname, "w")
            line = "mpirun -np 25 --machinefile hosts.txt python runme.py "
            line += " --numtr 4 --numreporter 2 --popsize 48 "
            line += " --patternpath ./UTESTS/rules." + baseid + ".txt "
            line += " --urspath ./UTESTS/urs." + baseid + ".txt "
            line += " --verbose 8 "
            line += " --runid " + runid 
            line += " --growth_rate 0.5 --decay_rate 0.5 "
            line += " --pwmpath ./UTESTS/pwm." + baseid + ".txt "
            line += " --workspace ./UTESTS --tfcoop zeros "
            line += " --maxgd 1 "
            line += " --iid_samples 5000 "
            #line += " --eliteproportion 0.2 "
            line += " --eliteproportion 0.0 "            
            line += " --elitemu 0.0 "
            line += " --mu 1 "
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
            line += " --maxgenerations 100 --maxtime 1 "
            fout.write(line + "\n")
            fout.close()
            os.system("source " + scriptname)
    exit()


#
# Null test:
#
# This test should not activate gene 2.  The error of the system can be estimated by running this
# example and observing the fluctions in pe values for gene 2.
#

def run_3genes():
    os.system("mkdir UTESTS")
    fout = open("UTESTS/3genes.run.sh", "w")
    fout.write("mpirun -np 25 --machinefile hosts.txt python runme.py --popsize 24 --stop_early 0.9 --patternpath ./UTESTS/3genes_rules.txt --urspath ./UTESTS/3genes_urs.txt --verbose 8 --runid 3genes.test --growth_rate 10.0 --decay_rate 10.0 --popsize 10 --pwmpath ./UTESTS/3genes_pwm.txt --maxgenerations 20 --maxtime 5 --workspace ./UTESTS --tfcoop zeros --maxgd 2 --iid_samples 10000  --mu 1 --eliteproportion 0.0 --elitemu 0.0 --p2pmu 0.0 --dbdmu 0.2 --pwmmu 1.0 --pwmlenmu 0.0 --pwmmulenmax 2 --urslenmu 0.0 --cismu 0.1 --sexual_ratio 0.7 --pe_scalar 0.5\n")
    fout.close()
    
    fout = open("./UTESTS/3genes_rules.txt", "w")
    fout.write("# timepatternID, basal_gene_id, timepoint, expression_level, reporter_gene_ID, rule_type, multiplier\n")
    fout.write("RULE 0 0 5 0.6 2 ge 1\n")
    fout.write("RULE 0 0 10 0.05 2 le 1\n")
    fout.write("# basal_gene_id, time_start, time_stop, expr_level\n")
    fout.write("INPUT 0 0 0 10 1.00\n")
    fout.close()
    
    fout = open("./UTESTS/3genes_urs.txt", "w")
    fout.write("# ID, has_dbd (0 or 1), is_repressor (0 or 1), URS\n")
    fout.write("0 1 0 TTTTTTTT\n")
    fout.write("1 1 1 AAAAAAAA\n")
    fout.write("2 0 0 AAAAAAAA\n")
    fout.close()
    
    fout = open("./UTESTS/3genes_pwm.txt", "w")
    fout.write("pwm 0\n") # wants A
    fout.write(get_random_pwm(3) + "\n")
    #fout.write("1.0 0.0 0.0 0.0\n")
    #fout.write("1.0 0.0 0.0 0.0\n")
    #fout.write("1.0 0.0 0.0 0.0\n")
    fout.write("pwm 1\n") # wants A
    fout.write(get_random_pwm(2) + "\n")
    #fout.write("1.0 0.0 0.0 0.0\n")
    #fout.write("1.0 0.0 0.0 0.0\n")
    fout.close()
    
    os.system("source UTESTS/3genes.run.sh")
    exit()

def run_tune_iid():
    """The idea here is to study how IID sample size affects error in estimating P(e).
    Given that the true value for P(e) is very difficult to estimate, our approach in SimGenome
    has been to estimate P(e) by taking many IID samples from the CDF of configurations.
    
    In this test, we plot IID sample size versus variance in the estimate of P(e) for all six
    genes in this setup, at generation 3.  This means that any variance or error in the estimate of P(e),
    due to poor IID samples, will propagate through generations 0 through 3.
    """
    baseid = "mu.01-15-2013iid"

    os.system("mkdir UTESTS")    
    fout = open("./UTESTS/rules." + baseid + ".txt", "w")
    fout.write("# timepatternID, basal_gene_id, timepoint, expression_level, reporter_gene_ID, rule_type, multiplier\n")
    fout.write("RULE 0 0 0 0.01 5 le\n")
    
    fout.write("# timepatternID, basal_gene_id, time_start, time_stop, expr_level\n")
    fout.write("INPUT 0 0 0 3 1.00\n")
    fout.close()
    
    fout = open("./UTESTS/urs." + baseid + ".txt", "w")
    fout.write("# ID, has_dbd (0 or 1), is_repressor (0 or 1), URS\n")
    fout.write("0 1 0 TCTGAGTTAACGAGTGGTTCATTAATAAGGTTATTTGACCGCTGGTGTGCATATGTGTGGTTTATAGCACACGGACCGTGTAGTTATACTGTATATACTA\n")
    fout.write("1 1 0 AGGCGCGGATGACCCCTGACCAATCGAGTACGATCGTCGAGTGCAAAAGAATCGAATATACCGAGTGGTGTGGGGCATGTTGCCTGTATATCCCATAGAA\n")
    fout.write("2 1 1 CGTCAAATATAATCGCGAGTCTGCGTCTTCACCGTCGGGGCCGGACCAACTATGAGGGTGCAATCAGTTTACCCAGTGGCGTTAATAATCCTGAAATGTG\n")
    fout.write("3 1 1 AGACAATAGACACCGCATATATGTCCCCCCCGATTTCGAGTCGTAGCACTCCATCTGGAACTTTTCGTACACACTAACTGATCGCTGGCATAAATTTCCT\n")
    fout.write("4 0 0 TGTCCTGCTTAAGTCCGCATTATGACTACGTTTCGTCCAACCCCCCCCAACAGATAAACTTTACATTATACACGGCGCTTGCTTACCAAAGAGCATCCTA\n")
    fout.write("5 0 0 TGACCTTACGGGAACATGAGAGCCAGGGGGAGTTATACAGGCCACCTCATAGTTTTGATGCATTCTCTGAACGGCACCCACCGATCCTGAGGTACTCGTA\n")
    fout.close()
    
    fout = open("./UTESTS/pwm." + baseid + ".txt", "w")
    fout.write("pwm 0\n")
    fout.write( "0.362475979737 0.20999765836 0.069210386986 0.358315974917 \n0.307385773957 0.262295147377 0.213728870626 0.21659020804 \n0.0623330092843 0.625133215743 0.169535710651 0.142998064322 \n0.290513390678 0.306784854161 0.135719312963 0.266982442198 \n0.104907631565 0.365796712482 0.301295803176 0.227999852778 \n0.135966865545 0.29601311261 0.385380226771 0.182639795073 \n" )
    fout.write("pwm 1\n")  
    fout.write( "0.28578029618 0.295508947432 0.0518462380452 0.366864518343 \n0.330437528702 0.0688417398484 0.376829880467 0.223890850983 \n0.354602173001 0.290280773331 0.193174824888 0.16194222878 \n0.0684562988648 0.430876488373 0.4725139575 0.0281532552624 \n0.27430529701 0.339132068053 0.185752012549 0.200810622387 \n0.230475154228 0.216958690077 0.297448309651 0.255117846044 \n" )
    fout.write("pwm 2\n")  
    fout.write( "0.194801808013 0.496655453248 0.00348184776187 0.305060890978 \n0.22196781785 0.0973674634052 0.277836329218 0.402828389527 \n0.0031653620689 0.52346521493 0.0814927748891 0.391876648112 \n0.218587466819 0.349034276796 0.180038181327 0.252340075059 \n0.19033068297 0.127085083428 0.522736735861 0.159847497741 \n0.333502058895 0.335992531304 0.317089768379 0.0134156414222 \n" )
    fout.write("pwm 3\n")   
    fout.write( "0.297540056275 0.358539293015 0.0710458253845 0.272874825325 \n0.266649436287 0.440621634587 0.0312018384081 0.261527090718 \n0.247725468452 0.222368134331 0.333033072885 0.196873324333 \n0.161450111065 0.0451144390746 0.441393429739 0.352042020122 \n0.10048915746 0.357337188367 0.222898669931 0.319274984242 \n0.236956817171 0.376004659325 0.0691984236151 0.317840099889 \n")
    fout.close()
    

    
    sample_sizes = [1, 10, 100,500,1000,2000,4000,6000,8000,10000]
    nreps = 10
    
    for n in sample_sizes:
        for rep in range(0, nreps):
            #print "utests.py 456", n, rep
            cismu = 0.01
            dbdmu = 0.5
            runid = baseid + "iid" + n.__str__() + "rep" + rep.__str__()
            scriptname = "UTESTS/" + baseid + ".run.sh"
            fout = open(scriptname, "w")
            line = "mpirun -np 2 --machinefile hosts.txt python runme.py "
            line += " --numtr 4 --numreporter 2 --popsize 1 "
            line += " --patternpath ./UTESTS/rules." + baseid + ".txt "
            line += " --urspath ./UTESTS/urs." + baseid + ".txt "
            line += " --verbose 8 "
            line += " --runid " + runid 
            line += " --growth_rate 0.5 --decay_rate 0.5 "
            line += " --pwmpath ./UTESTS/pwm." + baseid + ".txt "
            line += " --workspace ./UTESTS --tfcoop zeros "
            line += " --maxgd 1 "
            line += " --iid_samples " + n.__str__()
            #line += " --eliteproportion 0.2 "
            line += " --eliteproportion 0.0 "            
            line += " --elitemu 0.0 "
            line += " --mu 1 "
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
    
    fin = open("catch_iid.txt", "r")
    lines = fin.readlines()
    fin.close()
    
    n_gene_estimates = {}
    for n in sample_sizes:
        n_gene_estimates[n] = {}
        for i in range(0,6):
            n_gene_estimates[n][i] = []
    
    pt = 0
    for n in sample_sizes:
        for rep in range(0, nreps):
            #print "n",n,"rep", rep
            foundit = False
            while pt < lines.__len__() and foundit == False:
                if lines[pt].startswith("r:"):
                    tokens = lines[pt].split()
                    #print lines[pt]
                    #print tokens.__len__()
                    t = tokens[5]
                    if t == "3":
                        gene = int(tokens[9])
                        if tokens.__len__() == 17:
                            pe = float(tokens[12])
                        else:
                            pe = float(tokens[11])
                        #print "gene",gene,"pe", pe
                        n_gene_estimates[n][gene].append( pe )
                        if gene == 5:
                            foundit = True
                #if foundit == True:
                #    continue
                pt += 1
                
    #print n_gene_estimates

#    for i in range(0,6):
#        for n in sample_sizes:
#            if n_gene_estimates[n][i].__len__() > 0:
#                print "\t", i, n, n_gene_estimates[n][i]
#                s = ss.describe( n_gene_estimates[n][i] )
#                print "gene", i, "n",n, "size", s[0], "min %.3f"%s[1][0], "max %.3f"%s[1][1], "mean %.3f"%s[2], "var %.3e"%s[3], "stderr %.3e"%(s[3]/math.sqrt(s[0]))
#            else:
#                print i, n, "X"
    
    header = "\t"
    for i in range(0,6):
        header += i.__str__() + "\t"
    print header
    for n in sample_sizes:
        line = n.__str__() + "\t"
        for i in range(0,6):
            s = ss.describe( n_gene_estimates[n][i] )
            line += "%.3e"%(s[1][1] - s[1][0]) + "\t"
        print line

    header = "\t"
    for i in range(0,6):
        header += i.__str__() + "\t"
    print header
    for n in sample_sizes:
        line = n.__str__() + "\t"
        for i in range(0,6):
            s = ss.describe( n_gene_estimates[n][i] )
            line += "%.3e"%(s[2]) + "\t"
        print line

    header = "\t"
    for i in range(0,6):
        header += i.__str__() + "\t"
    print header
    for n in sample_sizes:
        line = n.__str__() + "\t"
        for i in range(0,6):
            s = ss.describe( n_gene_estimates[n][i] )
            line += "%.3e"%(s[1][0]) + "\t"
        print line

    header = "\t"
    for i in range(0,6):
        header += i.__str__() + "\t"
    print header
    for n in sample_sizes:
        line = n.__str__() + "\t"
        for i in range(0,6):
            s = ss.describe( n_gene_estimates[n][i] )
            line += "%.3e"%(s[1][1]) + "\t"
        print line
    
    exit()

def run_merge():
    """A production experiment.
    """
    nreps = 1
    n_vals = [16]
    m_vals = [16]
    
    nacts = 4
    nreps = 4
    
    for rep in range(0, nreps):
        for n in n_vals:
            m = n
            baseid = "merge.01-22-2013.n=%d"%n
            baseid = baseid + ".r" + rep.__str__()
        
            actids = []
            for i in range(2,2+nacts):
                actids.append(i)
            repids = []
            for i in range(2+nacts, 2+nacts+nreps):
                repids.append(i)
            nids = []
            for i in range(2+nacts+nreps, 2+nacts+nreps+n):
                nids.append(i)
            mids = []
            for i in range(2+nacts+nreps+n, 2+nacts+nreps+n+m):
                mids.append(i)
            
            os.system("mkdir UTESTS")    
            fout = open("./UTESTS/rules." + baseid + ".txt", "w")
            fout.write("# timepatternID, basal_gene_id, time_start, time_stop, expr_level\n")
            fout.write("INPUT 0 0 0 4 1.00\n")
            fout.write("INPUT 1 1 0 4 1.00\n")
            
            fout.write("# timepatternID, basal_gene_id, timepoint, expression_level, reporter_gene_ID, rule_type, multiplier\n")
            for nid in nids:
                fout.write("RULE 0 0 4 0.3 " + nid.__str__() + " ge\n")
                fout.write("RULE 0 0 8 0.01 " + nid.__str__() + " le\n")
            for mid in mids:
                fout.write("RULE 1 1 4 0.3 " + mid.__str__() + " ge\n")
                fout.write("RULE 1 1 8 0.01 " + mid.__str__() + " le\n")
            fout.close()
            
            fout = open("./UTESTS/urs." + baseid + ".txt", "w")
            fout.write("# ID, has_dbd (0 or 1), is_repressor (0 or 1), URS\n")
            fout.write("0 1 0 " + get_random_nt(200) + "\n")
            fout.write("1 1 0 " + get_random_nt(200) + "\n")
            for actid in actids:
                fout.write(actid.__str__() + " 1 0 " + get_random_nt(200) + "\n")            
            for repid in repids:
                fout.write(repid.__str__() + " 1 1 " + get_random_nt(200) + "\n")  
            for nid in nids:
                fout.write(nid.__str__() + " 0 0 " + get_random_nt(200) + "\n")
            for mid in mids:
                fout.write(mid.__str__() + " 0 0 " + get_random_nt(200) + "\n")
            fout.close()
            
            fout = open("./UTESTS/pwm." + baseid + ".txt", "w")
            fout.write("pwm 0\n")
            fout.write( get_efg1()  + "\n")
            fout.write("pwm 1\n")  
            fout.write( get_tec1()  + "\n")
            for actid in actids:
                fout.write("pwm " + actid.__str__() + "\n")
                fout.write( get_random_pwm(6) )
            for repid in repids:
                fout.write("pwm " + repid.__str__() + "\n")
                fout.write( get_random_pwm(6) )
            fout.close()
            
            cismu = 0.1
            dbdmu = 0.5
            runid = baseid
            scriptname = "UTESTS/" + baseid + ".run.sh"
            fout = open(scriptname, "w")
            line = "mpirun -np 25 --machinefile hosts.txt python runme.py "
            line += " --popsize 48 "
            line += " --patternpath ./UTESTS/rules." + baseid + ".txt "
            line += " --urspath ./UTESTS/urs." + baseid + ".txt "
            line += " --verbose 8 "
            line += " --runid " + runid 
            line += " --growth_rate 1.0 --decay_rate 1.0 "
            line += " --pwmpath ./UTESTS/pwm." + baseid + ".txt "
            line += " --workspace ./UTESTS --tfcoop zeros "
            line += " --maxgd 1 "
            line += " --iid_samples  8000"
            line += " --eliteproportion 0.2 "            
            line += " --elitemu 0.0 "
            line += " --mu 1 "
            line += " --cismu " + cismu.__str__() # VARIABLE
            line += " --dbdmu " + dbdmu.__str__() # VARIABLE
            line += " --pwmdeltamax 0.5 "
            line += " --pwmlenmu 0.0 --pwmmulenmax 2"
            line += " --urslenmu 0.0 "
            line += " --p2pmu 0.0"
            line += " --p2pmudelta 1.0"
            line += " --tfcoop zeros " # disable Co-op interactions
            line += " --sexual_ratio 0.8"
            line += " --maxgd 1"
            line += " --stop_early 0.9" # stop when the mean pop. f reaches 0.8
            line += " --maxgenerations 1000 --maxtime 1 "
            line += " --pe_scalar 0.01"
            fout.write(line + "\n")
            fout.close()
            os.system("source " + scriptname)
            #print scriptname

def run_kd_test():
    baseid = "mu.01-17-2013kd"

    os.system("mkdir UTESTS")    
    fout = open("./UTESTS/rules." + baseid + ".txt", "w")
    fout.write("# timepatternID, basal_gene_id, timepoint, expression_level, reporter_gene_ID, rule_type, multiplier\n")
    fout.write("RULE 0 0 0 0.01 5 le\n")
    
    fout.write("# timepatternID, basal_gene_id, time_start, time_stop, expr_level\n")
    fout.write("INPUT 0 0 0 10 1.00\n")
    fout.write("INPUT 0 1 0 10 1.00\n")
    fout.write("INPUT 0 2 0 10 1.00\n")
    fout.close()
    
    fout = open("./UTESTS/urs." + baseid + ".txt", "w")
    fout.write("# ID, has_dbd (0 or 1), is_repressor (0 or 1), URS\n")
    fout.write("0 1 0 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n")
    fout.write("1 1 0 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n")
    fout.write("2 1 0 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n")
    fout.write("3 0 0 AAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n")
    fout.write("4 0 0 TTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n")
    fout.write("5 0 0 GGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n")
    fout.close()
    
    fout = open("./UTESTS/pwm." + baseid + ".txt", "w")
    fout.write("pwm 0\n")
    fout.write(get_Aspec(10))
    fout.write("pwm 1\n")  
    fout.write(get_Tspec(5))
    fout.write("pwm 2\n")  
    fout.write(get_Gspec(2))
    fout.close()
    
    #print "utests.py 456", n, rep
    cismu = 0.00
    dbdmu = 0.0
    runid = baseid + "iid"
    scriptname = "UTESTS/" + baseid + ".run.sh"
    fout = open(scriptname, "w")
    line = "mpirun -np 2 --machinefile hosts.txt python runme.py "
    line += " --numtr 3 --numreporter 2 --popsize 1 "
    line += " --patternpath ./UTESTS/rules." + baseid + ".txt "
    line += " --urspath ./UTESTS/urs." + baseid + ".txt "
    line += " --verbose 8 "
    line += " --runid " + runid 
    line += " --growth_rate 1.0 --decay_rate 1.0 "
    line += " --pwmpath ./UTESTS/pwm." + baseid + ".txt "
    line += " --workspace ./UTESTS --tfcoop zeros "
    line += " --maxgd 1 "
    line += " --iid_samples 10000"
    #line += " --eliteproportion 0.2 "
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
    line += " --maxgenerations 1 --maxtime 10 "
    line += " --pe_scalar 0.1 " 
    fout.write(line + "\n")
    fout.close()
    os.system("source " + scriptname)
    


#run_ko()
#run_ffl()
#run_3genes()
#run_tune_mu()
#run_tune_iid()
#run_kd_test()
run_merge()
