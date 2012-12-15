"""USAGE: python plot_expression FOLDERPATH
...where FOLDERPATH is the EXPR_HISTORY folder for a project"""

import re, sys, os
import matplotlib.pyplot as plt
import math
import networkx as nx
from argparser import *

######################################################################
ap = ArgParser(sys.argv)
fpath = ap.getArg("--expdir")
generation = ap.getOptionalArg("--gen")
individual = ap.getOptionalArg("--gid")
MAX_EDGE_WIDTH = float( ap.getOptionalArg("--max_edge_width") )
if MAX_EDGE_WIDTH == False:
    MAX_EDGE_WIDTH = 15    
EXP_MODIFIER = float( ap.getOptionalArg("--exp_modifier") )
if EXP_MODIFIER == False:
    EXP_MODIFIER = 0.06
WEIGHT_CUTOFF = float( ap.getOptionalArg("--weight_cutoff") )
if WEIGHT_CUTOFF == False:
    WEIGHT_CUTOFF = 0.05
##########################################################################

def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i+lv/3], 16) for i in range(0, lv, lv/3))

def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % rgb

def float_to_rgb(f, max, one=None, two=200, three=200):
    if one == None:
        scalar = float(max * 1.0 / 200)
        return (int(f/scalar)%200 + 55,(int(f/scalar)%200 + 55),int(f/scalar)%200 + 55)

def normalize_edge_weights(G):
    """First, find the maximum weight"""
    max_w = 0.0
    for e in G.edges():
        w = G[e[0]][e[1]]['weight']
        if w > max_w:
            max_w = w
    """Next, normalize by the maximum."""
    for e in G.edges():
        G[e[0]][e[1]]['weight'] = MAX_EDGE_WIDTH * G[e[0]][e[1]]['weight'] #/ max_w
        #print G[e[0]][e[1]]['weight']
    return G

def build_graph(epath):
    """Build a directional graph, using the Network X package."""
    G = nx.MultiDiGraph()

    """Fill the graph."""
    fin = open(epath)
    this_time = None
    this_gene = None
    for l in fin.readlines():
        #print l
        if l.__len__() < 2:
            """Skip empty lines."""
            continue
        if l.startswith("."):
            """New time or gene information...."""
            tokens = l.split()
           # print tokens
            this_time = int(tokens[2])
            this_gene = int(tokens[4])
            if this_gene in G:
                G.add_node(this_gene)
        elif l.startswith("site") and this_time != None and this_gene != None:
            print l
            tokens = l.split()
            #print tokens
            if tokens.__len__() > 3:
                i = 3
                while (i < tokens.__len__()):
                    this_tf = int(tokens[i])
                    t = tokens[i+1]
                    p = float(tokens[i+2].split("=")[1])
                    b = float(tokens[i+3].split("=")[1])
                    #print p*b
                    this_w = math.exp(EXP_MODIFIER*p*b) - 1
                    if False == G.nodes().__contains__(this_tf):
                        G.add_node(this_tf, type=t)
                        #G[this_tf]['type'] = t
                    if this_w > WEIGHT_CUTOFF:
                        if False == G.has_edge(this_tf, this_gene):
                            #print "Adding edge", this_tf, this_gene
                            G.add_edge(this_tf, this_gene)
                            G[this_tf][this_gene]['weight'] = this_w
                            if t.__contains__("-"):
                                G[this_tf][this_gene]['type'] = "-"
                            elif t.__contains__("+"):
                                G[this_tf][this_gene]['type'] = "+"
                        else:
                            #print "+=", this_w
                            G[this_tf][this_gene]['weight'] += this_w
                    i += 4     
    fin.close()
    return G



def draw_network(epath):
    """epath is the path to a file printing the configuration."""
    outpath = re.sub("txt", "png", epath)
    print outpath
    
    """Fill the graph with data"""
    G = build_graph(epath)
        
    """Normalize edge weights."""
    #G = normalize_edge_weights(G)

    """Draw the graph."""
    pos=nx.shell_layout(G)
    
    for n in G.nodes():
        edge_sum = 0.0
        for e in G.edges():
            if e[1] == n:
                if G[e[0]][e[1]]['type'] == "-":
                    edge_sum -= G[e[0]][e[1]]["weight"]
                else:
                    edge_sum += G[e[0]][e[1]]["weight"]
        if edge_sum < 0:
            edge_sum = 0
        #color = rgb_to_hex( float_to_rgb(edge_sum, 20) )
        color = "white"
        #print edge_sum, "color=", color
        nx.draw_networkx_nodes(G,pos, nodelist=[ n ], node_color=color, node_size=1000)
    for e in G.edges():
        ecolor = "royalblue"
        if G[e[0]][e[1]]['type'] == "-":
            ecolor = "orangered"
        nx.draw_networkx_edges(G, pos, edgelist=[ e ], edge_color=ecolor, width=G[e[0]][e[1]]['weight'], arrows=True)
    nx.draw_networkx_labels(G,pos,font_size=14,font_family='Helvetica')
    plt.axis('off')
    plt.savefig(outpath) # save as png
    plt.close()


##########################################################
#
# main
#
##########################################################
files = os.listdir(fpath)
for f in files:
    if generation != False:
        """Skip this file if it's not for the generation we want."""
        if f.__contains__("gen" + generation.__str__()) == False:
            continue
    if individual != False:
        """SKip this file if its not for the individual we want."""
        if f.__contains__("gid" + individual.__str__()) == False:
            continue
    if f.__contains__("cran"):
        os.system("r --no-save < " + fpath + "/" + f)
    elif f.__contains__("txt"):
        draw_network(fpath + "/" + f)
