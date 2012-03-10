from configuration import *

#
# Write CRAN scripts and plot PDFs.
#
def plot_expression(genome, output_filename_seed, title, xlab, ylab):
    """plots the time vs. expression for all genes in one individual genome"""
    """genome must contain valid data in genome.gene_expr"""
    maxx = None
    maxy = None
    miny = None

    series = genome.gene_expr.keys()
    series.sort()

    # the ingredients for our CRAN script:
    plotstring = ""
    pointsstring = ""
    legendstring = "legend(\"bottomright\", c("
    for gid in series:
        legendstring += "\"" + gid.__str__() + "\","
    legendstring = re.sub(",$", "", legendstring)
    legendstring += "), cex=0.8, col=c("
    for gid in series:
        legendstring += "\"" + (gid+1).__str__().__str__() + "\","
    legendstring = re.sub(",$", "", legendstring)        
    legendstring += "), pch=15);\n"

    maxx = None
    string = ""
    string += "t <-c("
    for t in range(1, MAX_TIME):
        if maxx == None:
            maxx = t
        elif t > maxx:
            maxx = t
        string += t.__str__() + ","   
    string = re.sub(",$", "", string)
    string += ");\n"
    pointsstring = string

    miny = None
    maxy = None
    for gid in series:
        string = "y" + gid.__str__() + "<-c("
        for t in range(1, MAX_TIME):
            this_expr = genome.gene_expr[gid][t]
            if miny == None:
                miny = this_expr
            elif this_expr < miny:
                miny = this_expr
            if maxy == None:
                maxy = this_expr
            elif this_expr > maxy:
                maxy = this_expr
            string += this_expr.__str__() + ","
        string = re.sub(",$", "", string)
        string += ");\n"
        string += "points(t,y" + gid.__str__() + ", type='l', col='" + (gid+1).__str__() + "', lwd=3);\n"
        pointsstring += string
    
    plotstring = "xlimits <-c(0.0," + maxx.__str__() + ");\n"
    plotstring += "ylimits <-c(" + miny.__str__() + "," + maxy.__str__() + ");\n"
    plotstring += "plot(xlimits, ylimits, type='n', main='" + title + "', xlab='" + xlab + "', ylab='" + ylab + "');\n"
    
    fout_cran = open(output_filename_seed + ".cran", "w")
    fout_cran.write("pdf('" + output_filename_seed + ".pdf', width=7, height=4);\n")
    fout_cran.write(plotstring)
    fout_cran.write(pointsstring)
    fout_cran.write(legendstring)
    fout_cran.write("dev.off();\n")
    fout_cran.close()
    
    os.system("r --no-save < " + output_filename_seed + ".cran")