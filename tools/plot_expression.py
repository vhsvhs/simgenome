"""
This script plots the expression (time vs. concentration) for all genes.
By default, it will generate one of these plots for all individuals in all
generations, but this can be restricted by using --gen and --id

USAGE:
python plot_expression.py DIR
...where DIR is a simreg output directory with folder LOGS, FITNESS, etc.

Output: an R plot of time versus expression level, written into the PLOTS folder of DIR.

"""
import os, re, sys

from argparser import *
from test_common import *
from plot_includeme import *
ap = ArgParser(sys.argv)

outputdir = sys.argv[1] #outputdir is the folder into which a SimGenome run placed output.
if False == os.path.exists(outputdir + "/PLOTS"):
    os.system("mkdir " + outputdir + "/PLOTS")

generation = ap.getOptionalArg("--gen")
if generation != False:
    generation = int(generation)

indi = ap.getOptionalArg("--id")
if indi != False:
    indi = int(indi)
else:
    indi = -1


def get_expression_data(dir):
    tarr = []
    # data is  key = gene id, value = hash: key = time, value = expression level
    data = {}
    gene_mode = {}
    """Parses a text file containing the STDOUT from SimGenome."""
    expr_fin = open(dir + "/LOGS/expression.txt", "r")
    for l in expr_fin.xreadlines():
        if l.startswith("r:"):
            tokens = l.split()
            if tokens.__len__() < 11:
                continue
            rid = int( tokens[1] )
            genr = int( tokens[3] )
            if generation != False and generation != genr:
                continue
            t = int( tokens[5] )
            id = int( tokens[7] )
            
            #print id, indi
            if indi != -1 and id != indi: # skip this individual, possibly
                continue
            
            gid = int( tokens[9] )
            mode = tokens[10]
            if mode.startswith("expr:"):
                mode = ""
            if gid not in gene_mode:
                gene_mode[gid] = mode

            expr = None

            if t not in tarr:
                tarr.append(t)
            
            for ii in range(10, tokens.__len__()):
                if tokens[ii].startswith("expr:"):
                    expr = float( tokens[ii+1] )
            if generation == False or generation == genr:         
                if genr not in data:
                    data[genr] = {}
                if id not in data[genr]:
                    data[genr][id] = {}
                
                if rid not in data[genr][id]:
                    data[genr][id][rid] = {}
                if gid not in data[genr][id][rid]:
                    data[genr][id][rid][gid] = {}
                data[genr][id][rid][gid][t] = expr
    expr_fin.close()
    return (data,tarr,gene_mode)

def getr_expression( data, gen, id, tarr, gene_mode):
    tarr.sort()
    mint = None
    maxt = None
    ts = ""
    ts += "ts <-c("
    for t in tarr:
        ts += t.__str__() + ","
        if mint == None:
            mint = t
        elif mint > t:
            mint = t
        if maxt == None:
            maxt = t
        elif maxt < t:
            maxt = t
            
    ts = re.sub(",$", "", ts)
    ts += ");\n"
    
    miny = 0.01
    maxy = 10.0
    
    var_str = {}
    var_gene = {}
    for r in data[gen][id]:
        for g in data[gen][id][r]:
            var = "Gene" + g.__str__() + gene_mode[g]
            str = var + " <-c("
    
            for t in tarr:
                val = data[gen][id][r][g][t]
                if miny > val:
                    miny = val
                if maxy < val:
                    maxy = val
                str += val.__str__() + ","
            str = re.sub(",$", "", str)
            str += ");\n"
            var_str[var] = str
            var_gene[var] = g
    
    lstr = ""
    lstr += ts
    lstr +="par(xpd=TRUE);\n"
    #fout.write("par(fig=c(0,1.0,0,0.8), new=TRUE)\n")
    
    lstr +="x <-c(" + mint.__str__() + "," + maxt.__str__() + ");\n" 
    lstr +="y <-c(" + miny.__str__() + "," + maxy.__str__() + ");\n" 
    lstr +="plot(x,y,type='n', las=1, log='y', xlab='time', ylab='');\n"
    for var in var_str:
        lstr += var_str[var]
        lstr +="points(ts," + var + ",type='l',lwd='" + lwd_for_gene(var_gene[var]) + "', col=\"" + gene_color[var_gene[var]] + "\");\n"
        
    #fout.write("par(fig=c(0,1.0,0.55,1), new=TRUE)\n")
    
    #lstr = "legend(\"topleft\", c("
    lstr += "legend(0.0,1500.0, cex=0.7, ncol=5, c("
    varkeys = var_str.keys()
    varkeys.sort()
    for var in varkeys:
        lstr += "'" + var + "',"
    lstr = re.sub(",$", "", lstr)
    lstr += "), "
    lstr += "col=c("
    for var in varkeys:
        lstr += "\"" + gene_color[var_gene[var]] + "\","
    lstr = re.sub(",$", "", lstr)
    lstr += "), "
    lstr += "pch=c("
    for var in varkeys:
        lstr += "15,"
    lstr = re.sub(",$", "", lstr)
    lstr += ")"
    lstr += ");\n"
    
    return lstr


#########################################################################
#
# main
#
(expr_data,tarr,gene_mode) = get_expression_data(outputdir)
init_colors(gene_mode.keys().__len__())

for genr in expr_data:
    for id in expr_data[genr]:
        #plot_data(data, genr, id)
        foutpath = outputdir + "/PLOTS/" + "gen" + genr.__str__() + ".id" + id.__str__()
        fout = open(foutpath + ".rscript", "w")
        fout.write("pdf('" + foutpath + ".pdf', width=7, height=3.5);\n")
        
        # Time vs. Occupancy
        expr_string = getr_expression(expr_data, genr, id, tarr, gene_mode)
        fout.write(expr_string)
        
        if False == ap.getOptionalArg("--skip_binding"):
            # Sites vs. P(binding)
            for gene in gene_mode:
                for t in tarr:
                    timeslice = t
                    rid = 0
                    occupancy_path = outputdir + "/OCCUPANCY/occ.gen" + genr.__str__() + ".id" + id.__str__() + ".gene" + gene.__str__() + ".txt"
                    bind_string = getr_binding(occupancy_path, timeslice, rid)
                    fout.write(bind_string)
        
        fout.write("dev.off();\n")
        fout.write("par(xpd=TRUE);\n")
        fout.close()
        os.system("r --no-save < " + foutpath + ".rscript")