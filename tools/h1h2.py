#
# Measure H1 versus H2. See my lab notes for more info.
#
# Compare the K logs for two individuals from different generations.
#
from h1h2_include import *

ap.params["indir"] = ap.getArg("--dir") #outputdir is the folder into which a SimGenome run placed output.
ap.params["gen1"] = int( ap.getArg("--gen1") )
ap.params["gen2"] = int( ap.getArg("--gen2") )
ap.params["rids"] =  ap.getList("--rids") # desired regulatory problems
for i in range(0, ap.params["rids"].__len__()):
    ap.params["rids"][i] = int(ap.params["rids"][i])
ap.params["times"] =  ap.getList("--times") # desired time slices to compare
for i in range(0, ap.params["times"].__len__()):
    ap.params["times"][i] = int(ap.params["times"][i])
ap.params["reporters"] = ap.getList("--reporters")
for i in range(0, ap.params["reporters"].__len__()):
    ap.params["reporters"][i] = int(ap.params["reporters"][i])
ap.params["regulators"] = ap.getList("--regulators")
for i in range(0, ap.params["regulators"].__len__()):
    ap.params["regulators"][i] = int(ap.params["regulators"][i])

ap.params["id1"] = get_maxfit_id(ap.params["indir"], ap.params["gen1"])
ap.params["id2"] = get_maxfit_id(ap.params["indir"], ap.params["gen2"])

verify_files_exist()
      
print "\nDifference in Affinity"
d_mat = get_diff_matrix(d_op=operator.sub)
#print_mat(d_mat)
#h1h2_print(d_mat)
print h1h2_test(d_mat)[0]

print "\nAbsolute Difference in Affinity"
d_mat = get_diff_matrix(d_op=operator.sub, abs=True)
#print_mat(d_mat)
#h1h2_print(d_mat)
print h1h2_test(d_mat)[0]

print "\nRelative Difference Factor:"
d_mat = get_diff_matrix(d_op=operator.div)
#print_mat(d_mat)
#h1h2_print(d_mat)
print h1h2_test(d_mat)[0]

    