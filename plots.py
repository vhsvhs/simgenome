from configuration import *

"""Write CRAN scripts and plot PDFs."""
def plot_expression(genome, output_filename_seed, title, xlab, ylab, ap):
    """plots the time vs. expression for all genes in one individual genome"""
    """genome must contain valid data in genome.gene_expr"""
    maxx = None
    maxy = None
    miny = None

    series = genome.gene_expr.keys()
    series.sort()

    # the ingredients for our CRAN script:
    plotstring = ""
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
    tstring = ""
    tstring += "t <-c("
    for t in range(1, ap.params["maxtime"]):
        if maxx == None:
            maxx = t
        elif t > maxx:
            maxx = t
        tstring += t.__str__() + ","   
    tstring = re.sub(",$", "", tstring)
    tstring += ");\n"

    miny = None
    maxy = None
    tr_pointsstring = ""
    rep_pointsstring = ""
    for gid in series:
        string = "y" + gid.__str__() + "<-c("
        for t in range(1, ap.params["maxtime"]):
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
        if gid < ap.params["numtr"]:
            tr_pointsstring += string
        else:
            rep_pointsstring += string
    
    tr_colostring = "palette( rainbow(" + (ap.params["numtr"] + ap.params["numreporter"]).__str__() + ", start=0, end=4/6));\n"
    rep_colostring = "palette( rainbow(" + ap.params["numreporter"].__str__() + ", start=2/6, end=4/6));\n"
    
    plotstring = "xlimits <-c(0.0," + maxx.__str__() + ");\n"
    plotstring += "ylimits <-c(" + miny.__str__() + "," + maxy.__str__() + ");\n"
    plotstring += "plot(xlimits, ylimits, type='n', main='" + title + "', xlab='" + xlab + "', ylab='" + ylab + "');\n"
    
    fout_cran = open(output_filename_seed + ".cran", "w")
    fout_cran.write("pdf('" + output_filename_seed + ".pdf', width=7, height=4);\n")
    fout_cran.write(tstring)
    fout_cran.write(plotstring)
    fout_cran.write(tr_colostring)
    fout_cran.write(tr_pointsstring)
    #fout_cran.write(rep_colostring)
    fout_cran.write(rep_pointsstring)
    fout_cran.write(legendstring)
    fout_cran.write("dev.off();\n")
    fout_cran.close()
    
    if ap.getOptionalArg("--verbose") > 4:
        os.system("r --no-save < " + output_filename_seed + ".cran")