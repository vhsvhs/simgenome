"""
This script plots the expression (time vs. concentration) for all genes.
By default, it will generate one of these plots for all individuals in all
ap.params["generation"]s, but this can be restricted by using --gen and --id

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

ap.params["outputdir"] = ap.getArg("--outdir") #ap.params["outputdir"] is the folder into which a SimGenome run placed output.
if False == os.path.exists(ap.params["outputdir"] + "/PLOTS"):
    os.system("mkdir " + ap.params["outputdir"] + "/PLOTS")

#####################################################
#
# Get command-line arguments:
#
ap.params["generation"] = ap.getOptionalArg("--gen")
if ap.params["generation"] != False:
    ap.params["generation"] = int(ap.params["generation"])

ap.params["id"] = ap.getOptionalArg("--id")
if ap.params["id"] != False:
    ap.params["id"] = int(ap.params["id"])
#else:
#    print "You need to specify an individual with --id"
#    exit()

ap.params["genelist"] = ap.getList("--gene", type=int)

ap.params["rid"] = ap.getOptionalArg("--rid")
if ap.params["rid"] != False:
    ap.params["rid"] = int(ap.params["rid"])
else:
    print "\n. You didn't specify a rule ID, using --rid, so by default I'm using rule ID 0."
    ap.params["rid"] = 0


rulepath = ap.getOptionalArg("--rulepath") # path to rule file
gene_time_expr_type = {}
CALC_ERROR = False
if rulepath != False:
    CALC_ERROR = True
    fin = open(rulepath, "r")
    for l in fin.xreadlines():
        if l.startswith("RULE"):
            tokens = l.split()
            this_rid = int(tokens[1])
            if this_rid != ap.params["rid"]:
                continue
            
            this_gene = int(tokens[2])
            if this_gene not in gene_time_expr_type:
                gene_time_expr_type[this_gene] = {}
            this_time = int(tokens[3])
            if this_time not in gene_time_expr_type[this_gene]:
                gene_time_expr_type[this_gene][this_time] = []
            this_expr = float(tokens[4])
            this_type = int(tokens[5])
            this_weight = float(tokens[6])
            gene_time_expr_type[this_gene][this_time] = (this_expr, this_type, this_weight)


# ap.params["mode"] specifies if the output goes to an R plot, or to the command-line
# ap.params["mode"] can equal 'cran' or 'cli'
ap.params["mode"] = ap.getOptionalArg("--mode")
if ap.params["mode"] == False:
    ap.params["mode"] = "cran"
elif ap.params["mode"] != "cran" and ap.params["mode"] != "cli":
    ap.params["mode"] = "cran"
    

def resolve_cli_options():
    """If the user didn't specify all the necessary commands, then auto-determine
    the best options for those missing commands."""
    
    if ap.params["generation"] == False:
        fin = open(ap.params["outputdir"] + "/LOGS/generations.txt", "r")
        lines = fin.readlines()
        lastline = lines[ lines.__len__()-1 ]
        ap.params["generation"] = int( lastline.split()[1] )
        fin.close()
        print "\n. You didn't specify --gen for generation. I'm defaulting to the last generation (" + ap.params["generation"].__str__() + ")."
    if ap.params["genelist"] == False:
        fin = open(ap.params["outputdir"] + "/POPS/pop.gen" + ap.params["generation"].__str__() + ".save.txt", "r")
        for l in fin.xreadlines():
            if l.startswith("Genome:"):
                ngenes = int( l.split()[3] )
                ap.params["genelist"] = []
                for ii in range(0, ngenes):
                    ap.params["genelist"].append( ii )
                fin.close()
                print "\n. You didn't specify --gene for a list of genes. I'm defaulting to these genes:", ap.params["genelist"]
                
                break
    if ap.params["rid"] == False:
       ap.params["rid"] = 0 
    
    if ap.params["id"] == False:
        fin = open(ap.params["outputdir"] + "/FITNESS/fitness.gen" + ap.params["generation"].__str__() + ".txt", "r")
        minid = None
        minerr = None
        for l in fin.xreadlines():
            if l.__len__() > 5:
                tokens = l.split()
                this_id = int( tokens[0] )
                this_err = float( tokens[2] )
                if minerr == None:
                    minid = this_id
                    minerr = this_err
                if minerr > this_err:
                    minid == this_id
                    minerr = this_err
        fin.close()
        ap.params["id"] = minid
        print "\n. You didn't specify --id for genome ID. I'm defaulting to the max. fit individual" + ap.params["id"].__str__() + "."


def get_expression_data(dir):
    tarr = [] # array of time points in data
    data = {} # key = ap.params["generation"], value = hash; key = gene id, value = hash: key = time, value = expression level
    gene_mode = {} # key = gene id, value = "a" or "r" or None
    
    if CALC_ERROR == True:
        gene_error = {} # key = gene, value = error summed over all timepoints, using the rules defined in the rule path file.
    
    """Parses a text file containing the STDOUT from SimGenome."""
    expr_fin = open(dir + "/LOGS/expression.txt", "r")
    for l in expr_fin.xreadlines():
        if l.startswith("r:"):
            tokens = l.split()
            if tokens.__len__() < 11:
                continue
            if tokens.__len__() > 15:
                continue
            
            this_rid = int( tokens[1] )
            if this_rid != ap.params["rid"]:
                    continue
            
            genr = int( tokens[3] )
            if ap.params["generation"] != genr:
                continue
                        
            id = int( tokens[7] )
            if id != ap.params["id"]:
                continue
            
            gid = int( tokens[9] )
            if gid not in ap.params["genelist"] and ap.params["genelist"][0] != -1:
                continue
            
            ap.params["mode"] = tokens[10]
            if ap.params["mode"].startswith("expr:"):
                ap.params["mode"] = ""
            if gid not in gene_mode:
                gene_mode[gid] = ap.params["mode"]

            t = int( tokens[5] )
            if t not in tarr:
                tarr.append(t)

            expr = None
            
            for ii in range(10, tokens.__len__()):
                if tokens[ii].startswith("expr:"):
                    expr = float( tokens[ii+1] )
            if ap.params["generation"] == genr:         
                if genr not in data:
                    data[genr] = {}
                if id not in data[genr]:
                    data[genr][id] = {}
                
                if this_rid not in data[genr][id]:
                    data[genr][id][this_rid] = {}
                if gid not in data[genr][id][this_rid]:
                    data[genr][id][this_rid][gid] = {}
                data[genr][id][this_rid][gid][t] = expr
                            
                if CALC_ERROR == True:
                    if gid not in gene_error:
                        gene_error[gid] = 0.0
                    if gid in gene_time_expr_type:
                        if t in gene_time_expr_type[gid]:
                            if gene_time_expr_type[gid][t][1] == 0:
                                gene_error[gid] +=  gene_time_expr_type[gid][t][0] / expr - 1.0
                            elif gene_time_expr_type[gid][t][1] == 1:
                                gene_error[gid] +=  expr / gene_time_expr_type[gid][t][0] - 1.0
                    #
                    # to-do: normalize gene error by weight
                    #
                
    if CALC_ERROR == True:
        print "Gene\tSum of Error"
        for gene in gene_error:
            print gene, "\t%.3f"%gene_error[gene]
    
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
    (expr_data,tarr,gene_mode) = get_expression_data(ap.params["outputdir"])
    init_colors(gene_mode.keys().__len__())
    for genr in expr_data:
        for id in expr_data[genr]:
            #plot_data(data, genr, id)
            foutpath = ap.params["outputdir"] + "/PLOTS/" + "gen" + genr.__str__() + ".id" + id.__str__() + ".rid" + ap.params["rid"].__str__()
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
                        #this_ap.params["rid"] = 0
                        occupancy_path = ap.params["outputdir"] + "/OCCUPANCY/occ.gen" + genr.__str__() + ".id" + id.__str__() + ".gene" + g.__str__() + ".rid" + ap.params["rid"].__str__() + ".txt"
                        if os.path.exists(occupancy_path):
                            bind_string = getr_binding(occupancy_path, timeslice, ap.params["rid"])
                            fout.write(bind_string)    
            fout.write("dev.off();\n")
            fout.write("par(xpd=TRUE);\n")
            fout.close()
            os.system("r --no-save < " + foutpath + ".rscript")


def plot_cli():
    (expr_data,tarr,gene_mode) = get_expression_data(ap.params["outputdir"])
    
    #
    #expr_data is : [genr][id][this_ap.params["rid"]][gid][t]
    #
    ybins = [1.0, 0.1, 0.01, 0.001,0.0001,0.00001,0.000001, 0.0000001]    
    
    for genr in expr_data:
        for id in expr_data[genr]:   
            for this_rid in expr_data[genr][id]:
                for gid in expr_data[genr][id][ap.params["rid"]]:
                    yvals = []
                    for t in tarr:
                        yvals.append( float(expr_data[genr][id][this_rid][gid][t]) )
                    print "\nGenr:", genr, "ID:", id, "Reg.Problem:", this_rid, "Gene:", gid
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


#########################################################################
#
# main
#
resolve_cli_options()
if ap.params["mode"] == "cran":
    plot_cran()
elif ap.params["mode"] == "cli":
    plot_cli()
        
        
