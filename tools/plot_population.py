"""
plot_population.py

Plots the entire population, generations versus fitness.
The coalescent process is drawn.
"""


import os, re, sys
from argparser import *
ap = ArgParser(sys.argv)
outputdir = ap.getArg("--outputdir") #outputdir is the folder into which a SimGenome run placed output.

if False == os.path.exists(outputdir + "/PLOTS"):
    os.system("mkdir " + outputdir + "/PLOTS")

# All output files written by this script will be prefaced with
# the contents of the variable output_filename_seed
output_filename_seed = outputdir + "/PLOTS/popgraph"

title = "Population Fitness"
xlab = "generations"
ylab = "fitness"

#
# Read LOG/generations.txt
#
def build_plot_basics():
    """Reads LOGS/generations.txt in order to gather basic information about the population history,
    including min, max, and mean fitness at each generation."""
    maxx = None
    maxy = None
    genstring = ""
    genstring += "g <-c("
    maxfstring = "maxf <-c("
    minfstring = "minf <-c("
    meanfstring = "meanf <-c("

    fin = open(outputdir + "/generations.txt", "r")
    lines = fin.readlines()
    fin.close()
    
    count_gen = 0
    for l in lines:
        if l.__len__() > 1:
            tokens = l.split()
            gen = int(tokens[1])
            if maxx < gen:
                maxx = gen
            maxf = float(tokens[3])
            if maxy < maxf:
                maxy = maxf
            minf = float(tokens[5])
            meanf = float(tokens[7])
            
            genstring += gen.__str__() + ","   
            maxfstring += maxf.__str__() + ","
            minfstring += minf.__str__() + ","
            meanfstring += meanf.__str__() + ","
            count_gen += 1
    
    genstring = re.sub(",$", "", genstring)
    genstring += ");\n"
    maxfstring = re.sub(",$", "", maxfstring)
    maxfstring += ");\n"
    minfstring = re.sub(",$", "", minfstring)
    minfstring += ");\n"
    meanfstring = re.sub(",$", "", meanfstring)
    meanfstring += ");\n"
    
    plotstring = "plot(c(0.0," + maxx.__str__() + "), c(0," + maxy.__str__() + "), type='n', main='" + title + "', xlab='" + xlab + "', ylab='" + ylab + "');\n"
    
    pointsstring = "points(g,maxf,type='l', col='black', lwd=1);\n"
    pointsstring += "points(g,minf,type='l', col='black', lwd=1);\n"
    pointsstring += "points(g,meanf,type='l', col='black', lwd=3);\n"
    
    r = genstring + "\n"
    r += maxfstring + "\n"
    r += minfstring + "\n"
    r += meanfstring + "\n"
    r += plotstring + "\n"
    r += pointsstring + "\n"
    return [r, count_gen]

def build_poptrace(ngen):
    [paths,marks] = build_coalescent_paths(ngen)
    #print paths
    return plot_coalescent_paths(paths,marks)

def build_coalescent_paths(ngen):
    """Reads MATING/mating.gen... files in order to build a coalescent history
    of the population.
    And also reads FITNESS/fitness.gen... files in order to learn the fitness of
    every individual. """
    paths = [] # array of [ (x1,y1), (x2,y2) ] lists
    marks = [] # array of [ (x,y), char ]
    
    parent_fitness = {}
    for gen in range(0, ngen):
        my_parents = {}

        # put values in my_parents
        had_babies = []
        mpath = None
        if gen > 0:
            mpath = outputdir + "MATING/mating.gen" + (gen - 1).__str__() + ".txt"
            for line in open(mpath, "r").readlines():
                tokens = line.split()
                if line.__contains__("clone"):
                    p1 = tokens[4]
                    child = tokens[7]
                    my_parents[child] = (p1,p1)
                    if p1 not in had_babies:
                        had_babies.append(p1)
                elif line.__contains__("child"):
                    p1 = tokens[3]
                    p2 = tokens[6]
                    child = tokens[9]
                    my_parents[child] = (p1,p2)
                    if p1 not in had_babies:
                        had_babies.append(p1)
                    if p2 not in had_babies:
                        had_babies.append(p2)
                        

        children = []
        child_fitness = {}
        fpath = outputdir + "FITNESS/fitness.gen" + gen.__str__() + ".txt"
        for line in open(fpath, "r").readlines()[1:]:
            tokens = line.split()
            id = tokens[0]
            f = float(tokens[1])
            #print id, f
            child_fitness[id] = f
            children.append(id)
        
        if gen == 0:
            for child in children:
                my_parents[child] = (child,child)   
                
        #print "109:", gen    
        #print "110:", children
        #print "111:", child_fitness
        #print "112:", my_parents
                    
        if gen > 0:  
            for child in children:
                # add lines
                x1 = gen-1
                #print my_parents[child][0]
                y1 = parent_fitness[ my_parents[child][0] ]
                
                x2 = gen-1
                y2 = parent_fitness[ my_parents[child][1] ]
                
                x3 = gen
                y3 = child_fitness[child]

                paths.append( [(x1,y1),(x3,y3)] )
                paths.append( [(x2,y2),(x3,y3)] )

                if child not in had_babies:
                    if my_parents[child][0] == my_parents[child][1]: # clonal:
                        marks.append( [(gen-1, parent_fitness[child]), "w"] )
                    else:
                        marks.append( [(gen-1, parent_fitness[child]), "x"] )
                else:
                    if my_parents[child][0] == my_parents[child][1]: # clonal:
                        marks.append( [(gen-1, parent_fitness[child]), "c"] )
                    else:
                        marks.append( [(gen-1, parent_fitness[child]), "+"] )

        for child in children:
            parent_fitness[child] = child_fitness[child]
            
    return [paths, marks]         

def plot_coalescent_paths(paths, marks):
    ret = ""
    for path in paths:
        print path
        s =  "x <-c(" + path[0][0].__str__() + "," + path[1][0].__str__() + ");\n"
        s += "y <-c(" + path[0][1].__str__() + "," + path[1][1].__str__() + ");\n"
        
        if path[0][1] > path[1][1]:
            color = "indianred1"
            lwd = 0.05
        else:
            color = "royalblue1"  
            lwd = 0.15
        s += "points(x,y,type='l',col='" + color + "',lwd=" + lwd.__str__() + ");\n"
        ret += s
    for mark in marks:
        s = "points(c(" + mark[0][0].__str__() + "), c(" + mark[0][1].__str__() + "),"
        if mark[1] == "x":
            s +=" type='p',pch=4,col='red3',cex=0.35,ps=5,font=1);\n" #sexual death
        elif mark[1] == "w":
            s +=" type='p',pch=6,col='red1',cex=0.35,ps=5,font=1);\n" #clonal death
        elif mark[1] == "+":
            s +=" type='p',pch=1,col='royalblue3',cex=0.35,ps=5,font=1);\n" # mating
        else:
            s +=" type='p',pch=2,col='deepskyblue',cex=0.35,ps=5,font=1);\n" # clonal
        ret += s
    return ret

def get_legend():
    l = "legend(\"bottomright\", c('clonal birth', 'sexual birth', 'clonal death', 'sexual death'), col=c('deepskyblue','royalblue3','red1', 'red3'), pch=c(2,1,6,4), cex=0.6);\n"
    return l



def execute_cran_string(cstring):    
    fout_cran = open(output_filename_seed + ".cran", "w")
    fout_cran.write("pdf('" + output_filename_seed + ".pdf', width=6, height=5);\n")
    fout_cran.write(cstring)
    fout_cran.write("dev.off();\n")
    fout_cran.close()
    os.system("r --no-save < " + output_filename_seed + ".cran")


############################################
#
# main. . .
#
[s,ngen] = build_plot_basics()
p = build_poptrace(ngen)
l = get_legend()
final = s + "\n" + p + "\n" + l
execute_cran_string(final)



