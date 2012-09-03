############################
#
#
# probe_landscape.py
#
#
# INPUT: 
# - a pickled population
# - a fitness landscape
# - an individual in the population (presumably the best-fit individual)
#
# METHOD:
# - for each site in each cis sequence, mutate the site the other possible states,
#     and then recalculate (i) the individual's fitness, (ii) if there was a network
#    re-wiring event
#
# OUTPUT:
# - overall stats:
#    (a) proportion of 1-hop mutations that increase/decrease/unaffect fitness
#    (b) proportion of 1-hop mutations that rewire the network
#    (c) the intersection of data in (a) and (b)
#
# - a big table listing:
#    site, x->y, fitness change, rewiring change
#
# - (maybe) an graphic showing rewiring hotspots.
#
############################


from argparser import *
from configuration import *
from genetic_algorithm import *
from population import *
from version import *
from debug_tools import *
from cli import *

# read the pickled population

# extract the target genome

gene_site_state_result = {}
# for each gene:
#    for each site in the gene's URS:
#        for each of the possible alternate states:
#            make the mutation, calculate fitness
#            record the fitness in gene_site....
#            record the size of network shift in gene_site....

# Write a table expressing gene_site...


