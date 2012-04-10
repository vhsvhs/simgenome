from argparser import *
from configuration import *
from genetic_algorithm import *
from population import *
from version import *


#
# x is an array [population, landscape]
#
def get_ptable(x):
    tf_expr_level = {}
    for tf in range(0, ap.params["numtr"]):
        tf_expr_level[ x[0].genomes[0].genes[tf].id ] = 0.0
    ptables = ProbTable( ap.params["numtr"], MAX_GD, x[0].genomes[0].genes[0].urs.__len__() )
    ptables = x[1].calc_prob_tables(x[0].genomes[0], x[0].genomes[0].genes[0], tf_expr_level, ptables)
    return ptables
    
    