"""
This script plots the expression (time vs. concentration) for all genes.
By default, it will generate one of these plots for all individuals in all
generations, but this can be restricted by using --gen and --id

USAGE:
python plot_expression.py --outdir X
...where DIR is a simreg output directory with folder LOGS, FITNESS, etc.

Output: an R plot of time versus expression level, written into the PLOTS folder of DIR.

"""
import os, re, sys

from argparser import *
from test_common import *
from plot_includeme import *
ap = ArgParser(sys.argv)

outputdir = ap.getArg("--outdir") #outputdir is the folder into which a SimGenome run placed output.
if False == os.path.exists(outputdir + "/PLOTS"):
    os.system("mkdir " + outputdir + "/PLOTS")

#####################################################
#
# Get command-line arguments:
#
generation = ap.getArg("--gen")
if generation != False:
    generation = int(generation)
else:
    print "You need to specify a generation with --gen"
    exit()

indi = ap.getArg("--id")
if indi != False:
    indi = int(indi)
else:
    print "You need to specify an individual with --id"
    exit()

gene = ap.getArg("--gene")
if gene != False:
    gene = int(gene)
else:
    gene = -1

rid = ap.getOptionalArg("--rid")
if rid != False:
    rid = int(rid)
else:
    rid = 0

# mode specifies if the output goes to an R plot, or to the command-line
# mode can equal 'cran' or 'cli'
mode = ap.getOptionalArg("--mode")
if mode == False:
    mode = "cran"
elif mode != "cran" and mode != "cli":
    mode = "cran"
    


def get_expression_data(dir):
    tarr = [] # array of time points in data
    data = {} # key = generation, value = hash; key = gene id, value = hash: key = time, value = expression level
    gene_mode = {} # key = gene id, value = "a" or "r" or None
    
    """Parses a text file containing the STDOUT from SimGenome."""
    expr_fin = open(dir + "/LOGS/expression.txt", "r")
    for l in expr_fin.xreadlines():
        if l.startswith("r:"):
            tokens = l.split()
            if tokens.__len__() < 11:
                continue
            
            this_rid = int( tokens[1] )
            if this_rid != rid:
                    continue
            
            genr = int( tokens[3] )
            if generation != genr:
                continue
                        
            id = int( tokens[7] )
            if id != indi:
                continue
            
            gid = int( tokens[9] )
            if gid != gene and gene != -1:
                continue
            
            mode = tokens[10]
            if mode.startswith("expr:"):
                mode = ""
            if gid not in gene_mode:
                gene_mode[gid] = mode

            t = int( tokens[5] )
            if t not in tarr:
                tarr.append(t)

            expr = None
            
            for ii in range(10, tokens.__len__()):
                if tokens[ii].startswith("expr:"):
                    expr = float( tokens[ii+1] )
            if generation == genr:         
                if genr not in data:
                    data[genr] = {}
                if id not in data[genr]:
                    data[genr][id] = {}
                
                if this_rid not in data[genr][id]:
                    data[genr][id][this_rid] = {}
                if gid not in data[genr][id][this_rid]:
                    data[genr][id][this_rid][gid] = {}
                data[genr][id][this_rid][gid][t] = expr
                #print "87:", genr, id, this_rid, gid, t
    expr_fin.close()
    return (data,tarr,gene_mode)

def get_rstring( data, gen, id, tarr, gene_mode):
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
    
    var_str = {} # key = variable name for gene, R string
    var_lwd = {} # key = variable name for gene, value = line weight
    var_col = {} # key = variable name for gene, value = color
    for r in data[gen][id]:
        for g in data[gen][id][r]:
            var = "Gene" + g.__str__() + gene_mode[g]
            str = var + " <-c("
    
            for t in tarr:
                #print "119:", gen, id, r, g, t
                val = data[gen][id][r][g][t]

                if miny > val:
                    miny = val
                if maxy < val:
                    maxy = val
                str += val.__str__() + ","
            str = re.sub(",$", "", str)
            str += ");\n"
            var_str[var] = str
            if gene_mode[g] == "a" or gene_mode[g] == "r":
                var_lwd[var] = "1"
            else:
                var_lwd[var] = "2.1"
            var_col[var] = (g+1).__str__()
    
    lstr = ""
    lstr += ts
    lstr +="par(xpd=TRUE);\n"
    #fout.write("par(fig=c(0,1.0,0,0.8), new=TRUE)\n")
    
    lstr +="x <-c(" + mint.__str__() + "," + maxt.__str__() + ");\n" 
    lstr +="y <-c(" + miny.__str__() + "," + maxy.__str__() + ");\n" 
    lstr +="plot(x,y,type='n', las=1, log='y', xlab='time', ylab='');\n"
    for var in var_str:
        lstr += var_str[var]
        lstr +="points(ts," + var + ",type='l',lwd='" + var_lwd[ var ] + "', col=\"" + var_col[var] + "\");\n"
        
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
        lstr += "\"" + var_col[var] + "\","
    lstr = re.sub(",$", "", lstr)
    lstr += "), "
    lstr += "pch=c("
    for var in varkeys:
        lstr += "15,"
    lstr = re.sub(",$", "", lstr)
    lstr += ")"
    lstr += ");\n"
    
    return lstr


def plot_cran():
    (expr_data,tarr,gene_mode) = get_expression_data(outputdir)
    init_colors(gene_mode.keys().__len__())
    for genr in expr_data:
        for id in expr_data[genr]:
            #plot_data(data, genr, id)
            foutpath = outputdir + "/PLOTS/" + "gen" + genr.__str__() + ".id" + id.__str__()
            fout = open(foutpath + ".rscript", "w")
            fout.write("pdf('" + foutpath + ".pdf', width=7, height=3.5);\n")
            # Time vs. Occupancy
            expr_string = get_rstring(expr_data, genr, id, tarr, gene_mode)
            fout.write(expr_string)
            
            if False == ap.getOptionalArg("--skip_binding"):
                # Sites vs. P(binding)
                for g in gene_mode:
                    for t in tarr:
                        timeslice = t
                        this_rid = 0
                        occupancy_path = outputdir + "/OCCUPANCY/occ.gen" + genr.__str__() + ".id" + id.__str__() + ".gene" + g.__str__() + ".txt"
                        if os.path.exists(occupancy_path):
                            bind_string = getr_binding(occupancy_path, timeslice, this_rid)
                            fout.write(bind_string)    
            fout.write("dev.off();\n")
            fout.write("par(xpd=TRUE);\n")
            fout.close()
            os.system("r --no-save < " + foutpath + ".rscript")


def plot_cli():
    (expr_data,tarr,gene_mode) = get_expression_data(outputdir)
    
    #
    #expr_data[genr][id][this_rid][gid][t]
    #
    ybins = [1.0, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001,0.0000001]    
    
    for genr in expr_data:
        for id in expr_data[genr]:   
            for this_rid in expr_data[genr][id]:
                for gid in expr_data[genr][id][rid]:
                    yvals = []
                    for t in tarr:
                        yvals.append( float(expr_data[genr][id][this_rid][gid][t]) )
                    print "\nGenr:", genr, "ID:", id, "Problem:", this_rid, "Gene:", gid
                    #print yvals
                    for rowi in range(0,ybins.__len__()):
                        line = ybins[rowi].__str__() + "\t"
                        for timei in range(0, tarr.__len__()):
                            val = yvals[timei]
                            #print val
                            if rowi == 0:
                                if val >= ybins[0]:
                                    line += "x"
                                else:
                                    line += "."
                                #print val, rowi
                            elif val >= ybins[rowi] and val < ybins[rowi-1]:
                                line += "x"
                                #print val, rowi
                            else:
                                #print val
                                line += "."
                        print line 
                    print ""
                    
                        
                            
                        #print genr, id, this_rid, gid, t, y


#########################################################################
#
# main
#
if mode == "cran":
    plot_cran()
elif mode == "cli":
    plot_cli()
        
        
