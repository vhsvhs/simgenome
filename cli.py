from argparser import *
from configuration import *
from genetic_algorithm import *
from population import *
from version import *
from rules import *

def read_cli(ap):
    """This method reads the command-line interface arguments and sets global parameters"""
    x = ap.getArg("--workspace")
    ap.params["workspace"] = x

    x = ap.getArg("--runid")
    ap.params["runid"] = x

    x = ap.getOptionalArg("--verbose")
    if x != False:
        ap.params["verbosity"] = int(x)
    else:
        ap.params["verbosity"] = 2

    x = ap.getOptionalToggle("--keep") # Erase old data? (False= old data is left, it may be overwritten, but won't be erased). Default is to enable clearcache
    if x != False:
        ap.params["clearcache"] = False 
    else:
        ap.params["clearcache"] = True
    
    x = ap.getOptionalArg("--growth_rate")
    if x != False:
        ap.params["growth_rate"] = float(x)
    else:
        ap.params["growth_rate"] = MAX_TRANSCRIPTION_RATE

    x = ap.getOptionalArg("--iid_samples")
    if x != False:
        ap.params["iid_samples"] = int(x)
    else:
        ap.params["iid_samples"] = IID_SAMPLES
        
    x = ap.getOptionalArg("--decay_rate")
    if x != False:
        ap.params["decay_rate"] = float(x)
    else:
        ap.params["decay_rate"] = MAX_DECAY_RATE
    
    x = ap.getOptionalArg("--popsize")
    if x != False:
        ap.params["popsize"] = int(x)
    else:
        ap.params["popsize"] = N_GENOMES

    x = ap.getOptionalArg("--ngoals")
    if x != False:
        ap.params["ngoals"] = int(x)
    else:
        ap.params["ngoals"] = GOAL_COMPLEXITY    
        
    x = ap.getOptionalArg("--maxtime")
    if x != False:
        ap.params["maxtime"] = int(x)
    else:
        ap.params["maxtime"] = MAX_TIME

    x = ap.getOptionalArg("--maxgenerations")
    if x != False:
        ap.params["maxgens"] = int(x)
    else:
        ap.params["maxgens"] = MAX_GA_GENS

    x = ap.getOptionalArg("--numtr")
    if x != False:
        ap.params["numtr"] = int(x)
    else:
        ap.params["numtr"] = N_TR
    
    x = ap.getOptionalArg("--numreporter")
    if x != False:
        ap.params["numreporter"] = int(x)
    else:
        ap.params["numreporter"] = N_REPORTER 

    x = ap.getOptionalArg("--maxgd")
    if x != False:
        ap.params["maxgd"] = int(x)
        if ap.params["maxgd"] <= 0:
            ap.params["maxgd"] = 1
    else:
        ap.params["maxgd"] = MAX_GD
    # Here we precompte the range of GD values, because we'll use this range very often in the code. 
    ap.params["rangegd"] = []
    for d in range(0, ap.params["maxgd"]):
        ap.params["rangegd"].append( d )
    
    # Is the TF coop matrix
    ap.params["coopinit"] = None
    x = ap.getOptionalArg("--tfcoop")
    if x == False:
        ap.params["coopinit"] = "zeros"
    elif x == "random":
        ap.params["coopinit"] = "random"
    
    x = ap.getOptionalArg("--mu")
    if x != False:
        ap.params["mu"] = float(x)
    else:
        ap.params["mu"] = MU 
        
    x = ap.getOptionalArg("--elitemu")
    if x != False:
        ap.params["elitemu"] = float(x)
    else:
        ap.params["elitemu"] = ELITE_MU
        
    x = ap.getOptionalArg("--dbdmu") # prob. of mutating a DBD
    if x != False:
        ap.params["dbdmu"] = float(x)
    else:
        ap.params["dbdmu"] = DBD_MU

    x = ap.getOptionalArg("--pwmmu") # how far in bits can a PWM site be mutated?
    if x != False:
        ap.params["pwmmu"] = float(x)
    else:
        ap.params["pwmmu"] = PWM_MU

    x = ap.getOptionalArg("--pwmlenmu") # once selected for mutation, this is the prob. of inserting/deleting content to a DBD's PWM
    if x != False:
        ap.params["pwmlenmu"] = float(x)
    else:
        ap.params["pwmlenmu"] = PWM_LEN_MU

    x = ap.getOptionalArg("--pwmmulenmax") # the max length of a PWM indel.
    if x != False:
        ap.params["pwmmulenmax"] = float(x)
    else:
        ap.params["pwmmulenmax"] = PWM_MU_LEN_MAX

    x = ap.getOptionalArg("--urslenmu")
    if x != False:
        ap.params["urslenmu"] = float(x)
    else:
        ap.params["urslenmu"] = URS_LEN_MU

    x = ap.getOptionalArg("--cismu")
    if x != False:
        ap.params["cismu"] = float(x)
    else:
        ap.params["cismu"] = URS_LEN_MU

    x = ap.getOptionalArg("--p2pmu") 
    if x != False:
        ap.params["p2pmu"] = float(x)
    else:
        ap.params["p2pmu"] = P2P_MU

    x = ap.getOptionalArg("--p2pmudelta")
    if x != False:
        ap.params["p2pmudelta"] = float(x)
    else:
        ap.params["p2pmudelta"] = P2P_MU_DELTA

    x = ap.getOptionalArg("--eliteproportion")
    if x != False:
        ap.params["eliteproportion"] = float(x)
    else:
        ap.params["eliteproportion"] = ELITE_PROPORTION
        
    x = ap.getOptionalArg("--init_pwm_len")
    if x != False:
        ap.params["init_pwm_len"] = int(x)
    else:
        ap.params["init_pwm_len"] = INIT_PWM_LEN
        
    x = ap.getOptionalArg("--init_urs_len")
    if x != False:
        ap.params["init_urs_len"] = int(x)
    else:
        ap.params["init_urs_len"] = INIT_URS_LEN

    x = ap.getOptionalArg("--pe_scalar")
    if x != False:
        ap.params["pe_scalar"] = float(x)
    else:
        ap.params["pe_scalar"] = PE_SCALAR

    x = ap.getOptionalArg("--start_generation")
    if x != False:
        ap.params["generation"] = int(x)
    else:
        ap.params["generation"] = INIT_GEN

    x = ap.getOptionalArg("--sexual_ratio")
    if x != False:
        ap.params["sexual_ratio"] = float(x)
    else:
        ap.params["sexual_ratio"] = SEXUAL_RATIO

    ap.params["enable_epigenetics"] = False
    x = ap.getOptionalArg("--enable_heritable_expression")
    if x != False:
        if int(x) == 1:
            ap.params["enable_epigenetics"] = True
    
    x = ap.getOptionalArg("--ko")
    if x != False:
        ap.params["doko"] = True
        ap.params["kogenome"] = int(x)
    else:
        ap.params["doko"] = False
        ap.params["kogenome"] = -1
    
    ap.params["rangetrs"] = []
    """Here we precompute the range of TR numbers "rangertrs", and the range + 1 (for no dimerization) "rangetrs+".
    We'll use these ranges very often in the code, so this precomputation step saves time later."""
    for i in range(0, ap.params["numtr"]):
        ap.params["rangetrs"].append(i)
    ap.params["rangetrs+"] = ap.params["rangetrs"] + [ap.params["numtr"]]
    # . . . this produces an array [0,1,2,3,4,...,n] where 1 through n-1 correspond to TR indices, and n corresponds to the empty case. 


def check_world_consistency(ap, population, landscape):
    """Check for out-of-bounds genes"""
    for rid in landscape.rulecollections:
        for rule in landscape.rulecollections[rid].rules:
            if rule.reporter_id > population.genomes[0].genes.__len__() - 1:
                print "Ooops. One of your rules uses reporter gene", rule.reporter_id, "but your genome only contains", (population.genomes[0].genes.__len__() -1 ), "genes."
                print "I am quitting."
                exit(1)

def correct_maxtime_error(implied_max_time, ap):
        if implied_max_time > ap.params["maxtime"]:
            if comm.Get_rank() == 0 and ap.params["verbosity"] >= 1:
                print "\n. You specified the maximum timepoint to be", ap.params["maxtime"]
                print "  However, your configuration of inputs and fitness rules implies a maximum timepoint of", implied_max_time
                print "--> I am setting maximum time to", implied_max_time
            ap.params["maxtime"] = implied_max_time

def get_input_rules_from_file(ap):
    """See the examples, included with the source code, for information about the required file format."""
    if False == ap.getOptionalArg("--patternpath"):
        return None
    else:
        ret = {}  # see Landscape.rulecollection
        patternpath = ap.getOptionalArg("--patternpath")
        
        fin = open(patternpath, "r")
        lines = fin.readlines()
        fin.close()
        
        # 1, scan to find the maximum timepoint implied by the INPUT setup. . .
        for l in lines: 
            if l.startswith("RULE") or l.startswith("INPUT"):
                tokens = l.split()
                start_time = int(tokens[2])
                if start_time > ap.params["maxtime"]:
                    #print "cli.py 231", start_time, ap.params["maxtime"]
                    correct_maxtime_error(start_time, ap)
                end_time = int(tokens[3])
                if end_time > ap.params["maxtime"]:
                    #print "cli.py 235", end_time, ap.params["maxtime"]
                    correct_maxtime_error(end_time, ap)
                    
        # 2. scan to find the maximum timepoint implied by the RULE setup. . .
        max_time = 0
        for l in lines:
            if l.startswith("RULE"):
                tokens = l.split()
                if tokens.__len__() >= 6:
                    this_rc_id = int(tokens[1])
                    this_basal_gene_id = int(tokens[2])
                    this_timepoint = int(tokens[3])
                    this_expr_level = float(tokens[4])
                    if this_expr_level < MINIMUM_ACTIVITY_LEVEL:
                        this_expr_level = MINIMUM_ACTIVITY_LEVEL
                    if this_expr_level > MAXIMUM_ACTIVITY_LEVEL:
                        this_expr_level = MAXIMUM_ACTIVITY_LEVEL
                    this_reporter_gene_id = int(tokens[5])
                    this_rule_type = tokens[6]
                    if this_rule_type == "ge":
                        this_rule_type = operator.ge
                    elif this_rule_type == "eq":
                        this_rule_type = operator.eq
                    elif this_rule_type == "le":
                        this_rule_type = operator.le
                    else:
                        this_rule_type = this_rule_type = operator.ge
                    this_rule = Fitness_Rule(this_timepoint, this_expr_level, this_reporter_gene_id, this_rule_type)
                    if this_rc_id not in ret:
                        ret[this_rc_id] = Rule_Collection(this_rc_id)
                    ret[ this_rc_id ].rules.append( this_rule )
            elif l.startswith("INPUT"):
                tokens = l.split()
                if tokens.__len__() >= 4:
                    # first parse the parameters of the time pattern
                    this_rc_id = int(tokens[1])

                    this_basal_gene_id = int(tokens[2])
                    this_timepoint_start = int(tokens[3])
                    if this_timepoint_start > max_time:
                        max_time = this_timepoint_start
                    this_timepoint_stop = int(tokens[4])
                    if this_timepoint_stop > max_time:
                        max_time = this_timepoint_stop
                    this_expr_level = float(tokens[5])
                    if this_expr_level < MINIMUM_ACTIVITY_LEVEL:
                        this_expr_level = MINIMUM_ACTIVITY_LEVEL
                    if this_expr_level > MAXIMUM_ACTIVITY_LEVEL:
                        this_expr_level = MAXIMUM_ACTIVITY_LEVEL
                    # build the input pattern
                    if this_rc_id not in ret:
                        ret[this_rc_id] = Rule_Collection(this_rc_id)
                    for t in range(this_timepoint_start, this_timepoint_stop+1):
                        if t not in ret[this_rc_id].inputs:
                            ret[this_rc_id].inputs[t] = []
                        ret[this_rc_id].inputs[t].append( [this_basal_gene_id, this_expr_level] )
        
        
        # 3. Ensure we don't have rules for timepoints beyond the max timepoint
        correct_maxtime_error(max_time, ap)
        #print "cli.py 295 -", ap.params["maxtime"]
        return ret

def get_genes_from_file(ap):
    """Returns either a list of genes read from a file,
    or returns None if the user did not specify to use genes from a file."""
    
    if comm.Get_rank() == 0 and ap.params["verbosity"] >= 1:
        if False == ap.getOptionalArg("--urspath"): 
            print "\n. I found no value for --urspath."
            print "--> I will build random URSs."
        else:
            print "--> Reading the URSs described in", ap.getOptionalArg("--urspath")
                
        if False == ap.getOptionalArg("--pwmpath"): 
            print "\n. I found no value for --pwmpath."
            print "--> I will build random PWMs."
        else:
            print "--> Reading the PWMs described in", ap.getOptionalArg("--pwmpath")
    
    if False == ap.getOptionalArg("--urspath"):
        return None
    
    if False != ap.getOptionalArg("--urspath"):
        ret_genes = []
        urspath = ap.getOptionalArg("--urspath")
        fin = open(urspath, "r")
        lines = fin.readlines()
        
        #
        #
        #
#        count_tf = 0
#        count_reporter = 0
#        for l in lines:
#            if l.startswith("#"):
#                continue            
#            else:
#                tokens = l.split()
#                if tokens.__len__() >= 4:
#                    this_has_dbd = int(tokens[1])
#                    if this_has_dbd == 1:
#                        count_tf += 1
#                    else:
#                        count_reporter += 1
#        print count_tf, count_reporter
#        exit()
#        ap.params["numtr"] = count_tf
#        print "\n. I found " + count_tf.__str__() + " regulators in your gene file."
#        ap.params["numreporter"] = count_reporter
#        print "\n. I found " + count_reporter.__str__() + " reporters in your gene file."
#        
#        ap.params["rangetrs"] = []
#        for i in range(0, ap.params["numtr"]):
#            ap.params["rangetrs"].append(i)
#        ap.params["rangetrs+"] = ap.params["rangetrs"] + [ap.params["numtr"]]
        
        
        #
        #
        #
        for l in lines:
            if l.startswith("#"):
                continue
            else:
                tokens = l.split()
                if tokens.__len__() >= 4:
                    this_id = int(tokens[0])
                    this_pwm = None
                    if ap.getOptionalArg("--pwmpath"):
                        pwmpath = ap.getOptionalArg("--pwmpath")
                        this_pwm = PWM()
                        this_pwm.read_from_file( pwmpath, this_id )
                    this_has_dbd = int(tokens[1])
                    this_repressor = int(tokens[2])
                    if this_repressor == 0:
                        this_repressor = False
                    elif this_repressor == 1:
                        this_repressor = True
                    this_urs = tokens[3]
                    this_gene = Gene(this_id, this_urs.__len__(), urs=this_urs, has_dbd=this_has_dbd, repressor=this_repressor, pwm=this_pwm, ap=ap) 
                    ret_genes.append(this_gene)
        fin.close()
        
        return ret_genes
