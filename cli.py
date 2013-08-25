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

    x = ap.getOptionalArg("--maxgd") # How close can two proteins be?
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
        ap.params["mu"] = int(x)
    else:
        ap.params["mu"] = MU 
        
    x = ap.getOptionalArg("--elitemu")
    if x != False:
        ap.params["elitemu"] = float(x)
    else:
        ap.params["elitemu"] = ELITE_MU
        
    x = ap.getOptionalArg("--dbdmu") 
    """There will be dbdmu * n_tr mutation events, where each event includes
    changing one PWM site's preference by pwmmu amount, and an indel event with probability pwmlenmu.
    """
    if x != False:
        ap.params["dbdmu"] = float(x)
    else:
        ap.params["dbdmu"] = DBD_MU

    x = ap.getOptionalArg("--pwmdeltamax") # how far in bits can a PWM site be mutated?
    if x != False:
        ap.params["pwmdeltamax"] = float(x)
    else:
        ap.params["pwmdeltamax"] = PWM_MU

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
    """Equals proportion of genes in the genome that will have length mutation event.""" 
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
    
    ap.params["stopconvergence"] = False
    x = ap.getOptionalArg("--stop_early") # followed by the desired mean pop. fitness goal
    if x != False:
        ap.params["stopconvergence"] = True 
        ap.params["fgoal"] = float(x)
        if ap.params["fgoal"] > 1.0:
            ap.params["fgoal"] = 1.0
        if ap.params["fgoal"] < 0.0:
            ap.params["fgoal"] = 0.0
    
    x = ap.getOptionalArg("--ko")
    if x != False:
        ap.params["doko"] = True
        ap.params["kogenome"] = int(x)
    else:
        ap.params["doko"] = False
        ap.params["kogenome"] = -1
    
    ap.params["trlist"] = []
    """Here we precompute the range of TR numbers "rangertrs", and the range + 1 (for no dimerization) "trlist+".
    We'll use these ranges very often in the code, so this precomputation step saves time later."""
    for i in range(0, ap.params["numtr"]):
        ap.params["trlist"].append(i)
    ap.params["trlist+"] = ap.params["trlist"] + [ap.params["numtr"]]
    # . . . this produces an array [0,1,2,3,4,...,n] where 1 through n-1 correspond to TR indices, and n corresponds to the empty case. 

    #ap.params["gene_names"] = {} # key = gene ID (integer), value = name of gene (as a string)

def check_world_consistency(ap, population, landscape):
    """Check for out-of-bounds genes"""
    for rid in landscape.rulecollections:
        for rule in landscape.rulecollections[rid].rules:
            if rule.reporter_id > population.genomes[0].genes.__len__()-1:
                print "Ooops. One of your rules uses reporter gene", rule.reporter_id, "but your genome contains", (population.genomes[0].genes.__len__() ), "genes."
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
                time = int(tokens[3])
                if time > ap.params["maxtime"]:
                    #print "cli.py 231", start_time, ap.params["maxtime"]
                    correct_maxtime_error(time, ap)
            if l.startswith("INPUT"):
                end_time = int(tokens[4])
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
                    #this_basal_gene_id = int(tokens[2])
                    this_timepoint = int(tokens[3])
                    this_expr_level = float(tokens[4])
                    if this_expr_level < MINIMUM_ACTIVITY_LEVEL:
                        this_expr_level = MINIMUM_ACTIVITY_LEVEL
                    if this_expr_level > MAXIMUM_ACTIVITY_LEVEL:
                        this_expr_level = MAXIMUM_ACTIVITY_LEVEL
                    this_reporter_gene_id = int(tokens[2])
                    this_rule_type = tokens[5]
                    if this_rule_type == "ge":
                        this_rule_type = operator.ge
                    elif this_rule_type == "eq":
                        this_rule_type = operator.eq
                    elif this_rule_type == "le":
                        this_rule_type = operator.le
                    else:
                        this_rule_type = this_rule_type = operator.ge
                    this_multiplier = 1.0
                    if tokens.__len__() > 7:
                        this_multiplier = float(tokens[6])
                    if this_multiplier < 0.0:
                        this_multiplier = 0.0
                    this_rule = Fitness_Rule(this_timepoint, this_expr_level, this_reporter_gene_id, this_rule_type,m=this_multiplier)
                    if this_rc_id not in ret:
                        ret[this_rc_id] = Rule_Collection(this_rc_id)
                    ret[ this_rc_id ].rules.append( this_rule )
            elif l.startswith("INPUT"):
                tokens = l.split()
                if tokens.__len__() >= 4:
                    # first parse the parameters of the time pattern
                    this_rc_id = int(tokens[1]) # Rule Collection ID

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
        if False == ap.getOptionalArg("--pwmpath"): 
            print "\n. I found no value for --pwmpath."
            print "--> I will build random PWMs."
        else:
            print "--> Reading the PWMs described in", ap.getOptionalArg("--pwmpath")
    if False == ap.getOptionalArg("--urspath"):
        return None
    gene_pwms = {} # key = gene_id, value = PWM for this gene
    gene_regmode = {}
    fin = open(ap.getOptionalArg("--pwmpath"), "r")
    prev = None
    countsites = -1
    for l in fin.xreadlines():
        if l.__len__() < 1:
            continue
        if l.__contains__("pwm"):
            tokens = l.split()
            gene_id = int(tokens[1])
            prev = gene_id
            reg_mode = int( tokens[2] )
            gene_pwms[ gene_id ] = PWM()
            gene_regmode[ gene_id ] = reg_mode
            countsites = -1
        elif prev != None and l.__len__() > 2 and False == l.__contains__("#"):
            countsites += 1
            tokens = l.split()
            cc = 0
            gene_pwms[ prev ].P.append( {} )
            gene_pwms[ prev ].rangesites.append(countsites)
            for c in ALPHABET:
                gene_pwms[ prev ].P[countsites][c] = float(tokens[cc])
                cc += 1
    fin.close()
    
    
    if comm.Get_rank() == 0 and ap.params["verbosity"] >= 1:
        if False == ap.getOptionalArg("--urspath"): 
            print "\n. I found no value for --urspath."
            print "--> I will build random URSs."
        else:
            print "--> Reading the URSs described in", ap.getOptionalArg("--urspath")
    if False != ap.getOptionalArg("--urspath"):
        ret_genes = []
        urspath = ap.getOptionalArg("--urspath")
        fin = open(urspath, "r")
        lines = fin.readlines()

        for i in range(0, lines.__len__()):
            l = lines[i]
            if l.startswith("#"):
                continue
            elif l.startswith(">"):
                l = l.strip()
                tokens = l.split()
                gene_id = int( re.sub(">", "", tokens[0]) )
                gene_name = re.sub(tokens[0] + " ", "", l)
                if gene_name.__len__() < 1:
                    gene_name = gene_id.__str__()
                
                j = i+1
                this_urs = ""
                while j < lines.__len__():
                    if False == lines[j].startswith(">") and False == lines[j].startswith("#"):
                        this_urs += lines[j].strip()
                        j += 1
                    else:
                        j = lines.__len__()
                if gene_id in gene_pwms:
                    this_gene = Gene(gene_id, this_urs.__len__(), urs=this_urs, has_dbd=True, repressor=gene_regmode[ gene_id ], pwm=gene_pwms[ gene_id ], ap=ap, name=gene_name)
                else:
                    this_gene = Gene(gene_id, this_urs.__len__(), urs=this_urs, has_dbd=False, ap=ap, name=gene_name) 
                ret_genes.append(this_gene)
        fin.close()

    # for debugging:
    if ap.params["verbosity"] >= 99:
        print ". Gene Summary:"
        for gene in ret_genes:
            name = "(" + gene.name + ")"
            if gene.has_dbd:
                rep = "Activator"
                if gene.is_repressor:
                    rep = "Repressor"
                print " ", rep, gene.id, name, ": URS length:", gene.urs.__len__()
            else:
                print "  Reporter", gene.id, name, ": URS length:", gene.urs.__len__()
    
            
    return ret_genes


