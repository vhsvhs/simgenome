#
# Input: the captured STDOUT stream from a sim-reg run.
# Output: an R plot of time versus expression level.
#
import os, re, sys

from argparser import *
ap = ArgParser(sys.argv)

fin = open(sys.argv[1], "r")

generation = ap.getOptionalArg("--gen")
if generation != False:
    generation = int(generation)


# data is  key = gene id, value = hash: key = time, value = expression level

tarr = []
color = {}

data = {}
gene_mode = {}

def color_for_run(x):
    if x not in color:
        this_color = (color.__len__()+1).__str__()
        color[x] = this_color
    return color[x]

def get_data():
    """Parses a text file containing the STDOUT from SimGenome."""
    for l in fin.xreadlines():
        if l.startswith("r:"):
            tokens = l.split()
            rid = int( tokens[1] )
            genr = int( tokens[3] )
            if generation != False and generation != genr:
                continue
            t = int( tokens[5] )
            id = int( tokens[7] )
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
    fin.close()
    print data
    return data

def plot_data( data, gen, id):
    foutpath = sys.argv[1] + "gen" + gen.__str__() + ".id" + id.__str__()
    fout = open(foutpath + ".rscript", "w")
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
    
    miny = None
    maxy = 5.0
    
    var_str = {}
    for r in data[gen][id]:
        for g in data[gen][id][r]:
            var = "Gene" + g.__str__() + gene_mode[g]
            str = var + " <-c("
    
            for t in tarr:
                val = data[gen][id][r][g][t]
                if miny == None:
                    miny = val
                elif miny > val:
                    miny = val
                if maxy < val:
                    maxy = val
                str += val.__str__() + ","
            str = re.sub(",$", "", str)
            str += ");\n"
            var_str[var] = str
    
    fout.write("pdf('" + foutpath + ".pdf', width=7, height=3.5);\n")
    fout.write(ts)
    
    fout.write("x <-c(" + mint.__str__() + "," + maxt.__str__() + ");\n" )
    fout.write("y <-c(" + miny.__str__() + "," + maxy.__str__() + ");\n" )
    fout.write("plot(x,y,type='n',log='y');\n")
    for var in var_str:
        fout.write( var_str[var] )
        fout.write("points(ts," + var + ",type='l',col='" + color_for_run(var) + "');\n")
    
    lstr = "legend(\"topleft\", c("
    varkeys = var_str.keys()
    varkeys.sort()
    for var in varkeys:
        lstr += "'" + var + "',"
    lstr = re.sub(",$", "", lstr)
    lstr += "), "
    lstr += "col=c("
    for var in varkeys:
        lstr += color_for_run(var) + ","
    lstr = re.sub(",$", "", lstr)
    lstr += "), "
    lstr += "pch=c("
    for var in varkeys:
        lstr += "21,"
    lstr = re.sub(",$", "", lstr)
    lstr += ")"
    lstr += ");\n"
    fout.write(lstr);
    fout.write("dev.off();\n")
    fout.close()
    
    os.system("r --no-save < " + foutpath + ".rscript")
    
    
data = get_data()
for genr in data:
    for id in data[genr]:
        plot_data(data, genr, id)