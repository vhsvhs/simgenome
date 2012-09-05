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


def check_workspace(ap):
    if False == os.path.exists(WORKSPACE + "/" + ap.getArg("--runid")):
        os.system("mkdir " + WORKSPACE + "/" + ap.getArg("--runid"))
    dirs = ["PROBE"]
    for d in dirs:
        if False == os.path.exists(WORKSPACE + "/" + ap.getArg("--runid") + "/" + d):
            os.system("mkdir " + WORKSPACE + "/" + ap.getArg("--runid") + "/" + d)

def splash():
    print "======================================="
    print "."
    print ". Sim Genome - probe_landscape.py"
    print "."
    print ". " + VERSIONSTRING
    print "."
    print "======================================="
    

#
ap = ArgParser(sys.argv)
read_cli(ap)
ap.params["generation"] = 0
splash()
check_workspace(ap)


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


# calculate the original fitness
org_fitness = landscape.get_fitness( pop.genomes[genomeid], ap)

gene_site_state_result = {}
# for each gene:
for gene in genm.genes:
    print "\n. Evaluating all SNP mutations in the URS for gene", gene.id
    gene_site_state_result[gene.id] = {}
    # for each site in the gene's URS:
    for site in range( 0, gene.urs.__len__() ):
        print ". Evaluating site", site
        gene_site_state_result[gene.id][site] = {}
        # for each of the possible alternate states:
        org_state = gene.urs[site]
        for mut_state in ALPHABET:
            if mut_state != org_state:
                # make the mutation, calculate fitness
                new_urs = gene.urs[0:site] + mut_state + gene.urs[(site+1):gene.urs.__len__()]
                genm.genes[ gene.id ].urs = new_urs
                
                # (this next block is modified from genetic_algorithm "runsim_slave" method...
                new_fitness = landscape.get_fitness( pop.genomes[genomeid], ap)
                print org_state, "to", mut_state, "\tf_delta: %.4f"%(new_fitness - org_fitness) 
                
                # record the fitness in gene_site....
                gene_site_state_result[gene.id][site][mut_state] = new_fitness
                
                # record the size of network shift in gene_site....
                # to-do
                #
        # reset the URS to its original sequence
        genm.genes[ gene.id ].urs = gene.urs[0:site] + org_state + gene.urs[(site+1):gene.urs.__len__()]

# calculate properties of the fitness landscape
count_total_mu = 0  # all mutations
count_up_mu = 0     # adaptive mutations
count_down_mu = 0   # deleterious mutations
count_neu_mu = 0    # neutral mutations
for gene in genm.genes:
    for site in range( 0, gene.urs.__len__() ):
        org_state = gene.urs[site]
        for mut_state in ALPHABET:
            if mut_state != org_state:
                count_total_mu += 1
                if gene_site_state_result[gene.id][site][mut_state] > org_fitness:
                    count_up_mu += 1
                if gene_site_state_result[gene.id][site][mut_state] < org_fitness:
                    count_down_mu += 1
                if gene_site_state_result[gene.id][site][mut_state] == org_fitness:
                    count_neu_mu += 1

# Write a table expressing gene_site...
fout = open(WORKSPACE + "/" + ap.getArg("--runid") + "/PROBE/probe_results.txt", "w")
for gene in genm.genes:
    for site in range( 0, gene.urs.__len__() ):
        org_state = gene.urs[site]
        for mut_state in ALPHABET:
            if mut_state != org_state:
                line = "gene: " + gene.id.__str__()
                line += "\tsite: " + site.__str__()
                line += "\t" + org_state + " to " + mut_state
                line += "\tf_delta: %.4f"%(gene_site_state_result[gene.id][site][mut_state] - org_fitness)
                print line
                fout.write(line + "\n")
line = "\n=================================================\n"
line += count_total_mu.__str__() + " possible cis mutations:\n"
line += "%.3f are adaptive.\n"%(float(count_up_mu) / count_total_mu)
line += "%.3f are deleterious.\n"%(float(count_down_mu) / count_total_mu)
line += "%.3f are neutral.\n"%(float(count_neu_mu) / count_total_mu)
print line
fout.write(line)
fout.close()