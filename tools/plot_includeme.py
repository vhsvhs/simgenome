import os, re, sys
from test_common import *


def getr_binding(path, timeslice, rid):
    """Returns an rscript that can be embedded in other r scripts,
    with code to plot the site vs. prob(expr) peaks for all transcription
    factors binding all sites in the gene specified by path,
    for the timeslice and regulatory problem (rid)."""
    
    #
    # Find the lines that are relevant to this timeslice and gene
    #
    foundit = False
    fin = open(path, "r")
    relevant_lines = []
    last_time = 0
    last_rid = 0
    for l in fin.readlines():
        if l.startswith(". time"):
            tokens = l.split()
            last_time = int(tokens[2])
            last_rid = int(tokens[4])
        if last_time == timeslice and last_rid == rid and l.startswith("site"):
            relevant_lines.append(l)
    fin.close()
    
    
    #
    # Fill site_tf_expre
    #
    site_tf_expr = {}
    tf_site_expr = {}
    for l in relevant_lines:
        tokens = l.split()
        this_site = int(tokens[1])
        tokens = tokens[3:]
        i = 0
        while i < tokens.__len__():
            this_tf = int(tokens[i])
            this_p = float(tokens[i+2])
            if this_tf not in tf_site_expr:
                tf_site_expr[this_tf] = {}
            tf_site_expr[this_tf][this_site] = this_p
            
            if this_site not in site_tf_expr:
                site_tf_expr[this_site] = {}
            site_tf_expr[this_site][this_tf] = this_p
            i += 3
            #i += 6
    
    #
    # for testing:
    #
    #for site in site_tf_expr:
    #    print site, site_tf_expr[site]
    
    #for tf in tf_site_expr:
    #    print tf, tf_site_expr[tf]
    
    maxy = 1.0
    miny = 0.0
    maxx = 0 # N sites
    minx = 0
    lout = "" # line out
    
    # scan for the max timepoint
    for tf in tf_site_expr:
        for site in tf_site_expr[tf]:
            if site > maxx:
                maxx = site
    
    for tf in tf_site_expr:
        l = "tf" + tf.__str__() + " <-c("    
        lastval = None
        for site in range(1, maxx+1):
            if site in tf_site_expr[tf]:
                lastval = tf_site_expr[tf][site]
                l += "%.3f"%tf_site_expr[tf][site]
            else:
                lastval = lastval / 2 # decay function
                l += "%.3f"%lastval
            l += ","
        l = re.sub(",$", "", l)
        l += ");\n"
        lout += l
    
    l = "x <-c("
    for site in range(1, maxx+1):
        l += site.__str__() + ","
    l = re.sub(",$", "", l)
    l += ");\n"
    lout += l
    
    title = ""
    xlab = "Sites"
    ylab = "P"
    lout += "par(xpd=TRUE);\n"
    lout += "plot(c(1," + maxx.__str__() + "), c(0," + maxy.__str__() + "), type='n',xaxt='n',las=1, main='" + title + "', xlab='', ylab='');\n"
    lout += "axis(side=1, at=" + minx.__str__() + ":" + maxx.__str__() + ",cex.axis=0.5);\n"
    lout += "mtext(\"sites\", side=1, line=3);\n"
    lout += "mtext(\"P\", side=2, las=1, line=3);\n"
    
    
    for tf in tf_site_expr:
        lout += "points(x, tf" + tf.__str__() + ",type='l',lwd='" + lwd_for_gene(tf) + "', col='" + gene_color[tf] + "');\n"
    
    lstr = "legend(0,1.5, cex=0.7, ncol=5, c("
    for tf in tf_site_expr:
        lstr += "'" + tf.__str__() + "',"
    lstr = re.sub(",$", "", lstr)
    lstr += "), "
    lstr += "col=c("
    for tf in tf_site_expr:
        lstr += gene_color[tf] + ","
    lstr = re.sub(",$", "", lstr)
    lstr += "), "
    lstr += "pch=c("
    for tf in tf_site_expr:
        lstr += "15,"
    lstr = re.sub(",$", "", lstr)
    lstr += ")"
    lstr += ");\n"
    
    #lout += lstr
    
    lout += "mtext(\"" + path + " - t=" + timeslice.__str__() + " - r=" + rid.__str__() + "\", side=3,cex=0.7);\n"
    
    return lout