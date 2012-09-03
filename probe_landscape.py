############################
#
#
# probe_landscape.py
#
#
# INPUT: 
# - a pickled population (--poppath)
# - a fitness landscape (--patternpath)
# - an individual in the population (presumably the best-fit individual) (--genomeid)
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
from cli import *

#
ap = ArgParser(sys.argv)
read_cli(ap)
ap.params["generation"] = 0

# read the pickled population
poppath = ap.getArg("--poppath")
print "\n. Loading the population in the file", poppath
popfile = open(poppath, "r")
p0 = pickle.load( popfile )
popfile.close()
pop = Population()
pop.uncollapse(p0)

# extract the target genome
genomeid = int( ap.getArg("--genomeid"))
genm = pop.genomes[genomeid]

# build the fitness landscape
landscape = Landscape(ap)
landscape.init(ap, genome = genm)

gene_site_state_result = {}
# for each gene:
for gene in genm.genes:
    print "\n. Evaluating 1 bp mutations for gene", gene.id
    #    for each site in the gene's URS:
    for site in range( 0, gene.urs.__len__() ):
        print ". Evaluating site", site
        #        for each of the possible alternate states:
        org_state = gene.urs[site]
        for mut_state in ALPHABET:
            if mut_state != org_state:
                #            make the mutation, calculate fitness
                new_urs = gene.urs[0:site] + mut_state + gene.urs[(site+1):gene.urs.__len__()]
                genm.genes[ gene.id ].urs = new_urs
                
                # (this next block is modified from genetic_algorithm "runsim_slave" method...
                new_fitness = landscape.get_fitness( pop.genomes[genomeid], ap)
                print org_state, "->", mut_state, "f= ", new_fitness 
                #if ap.params["verbosity"] >= 2:
                    #"""and then plot the expression of all genes in this genome"""
                    #filenameseed = ap.getArg("--runid") + "/" + EXPR_PLOTS + "/expr.gen" + i.__str__() + ".gid" + gid.__str__() 
                    #plot_expression( genm, filenameseed, filenameseed, "time", "expression", ap)
                
                #            record the fitness in gene_site....
                gene_site_state_result[gene][site][mut_state] = new_fitness
                
                #            record the size of network shift in gene_site....
                # to-do
                #
        # reset the URS to its original sequence
        genm.gene[ gene.id ].urs = gene.urs[0:site] + org_state + gene.urs[(site+1):gene.urs.__len__()]
                
# Write a table expressing gene_site...
