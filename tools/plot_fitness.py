"""
    (This script is depricated, replaced by plot_population.py)

    USAGE: 

    >$ python plot_fitness.py PATH

    where PATH is a file path to a generations.txt log file,
    found in the LOG directory.
"""

import os, re, sys

output_filename_seed = sys.argv[1]
title = "Fitness (" + output_filename_seed + ")"
xlab = "generations"
ylab = "fitness"

fin = open(output_filename_seed, "r")
lines = fin.readlines()
fin.close()

maxx = None
maxy = None
genstring = ""
genstring += "g <-c("
maxfstring = "maxf <-c("
minfstring = "minf <-c("
meanfstring = "meanf <-c("
for l in lines:
    if l.__len__() > 1:
        print l
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

genstring = re.sub(",$", "", genstring)
genstring += ");\n"
maxfstring = re.sub(",$", "", maxfstring)
maxfstring += ");\n"
minfstring = re.sub(",$", "", minfstring)
minfstring += ");\n"
meanfstring = re.sub(",$", "", meanfstring)
meanfstring += ");\n"

plotstring = "plot(c(0.0," + maxx.__str__() + "), c(0," + maxy.__str__() + "), type='n', main='" + title + "', xlab='" + xlab + "', ylab='" + ylab + "');\n"

pointsstring = "points(g,maxf,type='l', col='black', lwd=1.3);\n"
pointsstring += "points(g,minf,type='l', col='black', lwd=1.3);\n"
pointsstring += "points(g,meanf,type='l', col='blue', lwd=3);\n"

fout_cran = open(output_filename_seed + ".cran", "w")
fout_cran.write("pdf('" + output_filename_seed + ".pdf', width=7, height=4);\n")
fout_cran.write(genstring)
fout_cran.write(maxfstring)
fout_cran.write(minfstring)
fout_cran.write(meanfstring)
fout_cran.write(plotstring)
fout_cran.write(pointsstring)
fout_cran.write("dev.off();\n")
fout_cran.close()

os.system("r --no-save < " + output_filename_seed + ".cran")


