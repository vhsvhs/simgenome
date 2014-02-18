#
# Python script functions for dealing with SimReg population files.
#

import os, sys, re

def get_individual(poppath, id):
    """Returns array of strings from pop.save.txt. array[0] is the header string,
    and array[1] is the body string for the pop.save.txt"""
    if os.path.exists(poppath) == False:
        print "Oops. I can't find the poppath file ", poppath
        exit()
    out = []        
    header = []
    body = []
    fin = open(poppath, "r")
    ingenome = False
    foundbody = False
    for l in fin.xreadlines():
        l = l.strip()
        if l.__contains__("total genes"):
            this_id = int( l.split()[1] )
            if this_id != id:
                ingenome = False
                continue
            else:
                ingenome = True
                header.append( l )
        elif l.__contains__("Genome " + id.__str__() + " :"):
            foundbody = True
            body.append(l)
            ingenome = True
        elif ingenome and foundbody == False:
            header.append( l )
        elif ingenome and foundbody == True:
            body.append( l )
        # reset on empty lines
        if l.__len__() <= 1:
            ingenome = False
    return [ header, body ]

def clone_individual(lines, id, n):
    out = []
    for i in range(0, n):
        for l in lines[0]: # for each header line:
            outl = l
            outl = re.sub("Genome: " + id.__str__(), "Genome: " + i.__str__(), outl)
            out.append(outl)
    for i in range(0, n):
        for l in lines[1]: # for each body line:
            outl = l
            outl = re.sub("Genome " + id.__str__() + " :", "Genome " + i.__str__() + " :", outl)
            out.append(outl)            
    return out
            