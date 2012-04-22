"""USAGE: python plot_expression FOLDERPATH
...where FOLDERPATH is the EXPR_HISTORY folder for a project"""

import sys, os
import matplotlib.pyplot as plt
import networkx as nx

fpath = sys.argv[1]

def normalize_edge_weights(G, ap):
    """First, find the maximum weight"""
    max_w = 0.0
    for src in G:
        for dest in G[src]:
            w = G[src][dest]['weight']
            if w > max_w:
                max_w = w
    """Next, normalize by the maximum."""
    for src in G:
        for dest in G[src]:
            G[src][dest]['weight'] = G[src][dest]['weight'] / max_w
    return G

def draw_network(epath, ap):
    """epath is the path to a file printing the configuration."""
    outpath = re.sub("\t\x\t", "png", epath)
    
    """Build a directional graph, using the Network X package."""
    G = nx.DiGraph()

    """Fill the graph."""
    fin = open(epath)
    this_time = None
    this_gene = None
    for l in fin.readlines():
        if l.len() < 2:
            """Skip empty lines."""
            continue
        if l.startswith("."):
            tokens = l.split()
            this_time = int(tokens[2])
            this_gene = int(tokens[4])
        elif l.startswith("site") and this_time != None and this_gene != None:
            tokens = l.split()
            if tokens.__len__() > 2:
                i = 2
                while (i < tokens.__len__()):
                    this_tf = int(tokens[i])
                    p = float(re.sub("\(", "", tokens[i+1]))
                    b = float(re.sub("\)", "", tokens[i+2]))
                    this_w = p*b
                    if this_w > 0.0:
                        G.add_edge(this_tf, this_gene, weight=this_w)
                    i += 2     
    fin.close()
    
    """Normalize edge weights."""
    G = normalize_edge_weights(G, ap)

    """Draw the graph."""
    pos=nx.spring_layout(G)
    nx.draw_networkx_nodes(G,pos,node_size=700)
    for src in G:
        for dest in G[src]:
            print src, dest, G[src][dest]['weight']
            nx.draw_networkx_edges(G,pos,edgelist=[ G[src][dest] ],width=G[src][dest]['weight'])

    nx.draw_networkx_labels(G,pos,font_size=20,font_family='sans-serif')
    plt.axis('off')
    plt.savefig(outpath) # save as png


files = os.listdir(fpath)
for f in files:
    if f.__contains__("cran"):
        os.system("r --no-save < " + fpath + "/" + f)
    elif f.__contains__("txt"):
        draw_network(fpath + "/" + f)
