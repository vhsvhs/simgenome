from argparser import *
from configuration import *
from genetic_algorithm import *
from population import *
from version import *
from debug_tools import *

def read_cli(ap):
    """Read the CLI args and set global parameters"""
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
        
    x = ap.getOptionalArg("--pwmmu")
    if x != False:
        ap.params["pwmmu"] = float(x)
    else:
        ap.params["pwmmu"] = PWM_MU

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


def check_world_consistency(ap, population, landscape):
    """Check for out-of-bounds genes"""
    for tp in landscape.timepatterns:
        if tp.basal_gene_id > population.genomes[0].genes.__len__() - 1:
            print "Ooops. One of your rules uses basal gene", tp.basal_gene_id, "but your genome only contains genes 0 through", (population.genomes[0].genes.__len__() - 1)
            exit(1)
        for rule in tp.rules:
            if rule.reporter_id > population.genomes[0].genes.__len__() - 1:
                print "Ooops. One of your rules uses reporter gene", rule.reporter_id, "but your genome only contains", population.genomes[0].genes.__len__(), "genes."
                exit(1)

def get_timepatterns_from_file(ap):
    if ap.getOptionalArg("--patternpath"):
        ret_timepatterns = []
        patternpath = ap.getOptionalArg("--patternpath")
        fin = open(patternpath, "r")
        for l in fin.readlines():
            if l.startswith("#"):
                continue
            else:
                tokens = l.split()
                if tokens.__len__() >= 6:
                    this_timepattern_id = int(tokens[0])
                    this_basal_gene_id = int(tokens[1])
                    this_timepoint = int(tokens[2])
                    this_expr_level = float(tokens[3])
                    this_reporter_gene_id = int(tokens[4])
                    this_rule_type = tokens[5]
                    if this_rule_type == "ge":
                        this_rule_type = operator.ge
                    elif this_rule_type == "eq":
                        this_rule_type = operator.eq
                    elif this_rule_type == "le":
                        this_rule_type = operator.le
                    else:
                        this_rule_type = this_rule_type = operator.ge

                    if ret_timepatterns.__len__() <= this_timepattern_id:
                        this_timepattern = Time_Pattern(this_basal_gene_id)
                        ret_timepatterns.append(this_timepattern)
                
                    if this_timepoint > ap.params["maxtime"]:
                        ap.params["maxtime"] = this_timepoint
                                            
                    this_rule = Fitness_Rule(this_timepoint, this_expr_level, this_reporter_gene_id, this_rule_type)
                    ret_timepatterns[ this_timepattern_id ].rules.append( this_rule )
        return ret_timepatterns
        fin.close()

"""Returns either a list of genes read from a file, OR returns None if the user
did not specify to use genes from a file."""
def get_genes_from_file(ap):
    #print "Getting genes from file..."
    if ap.getOptionalArg("--urspath"): 
        ret_genes = []
        urspath = ap.getOptionalArg("--urspath")
        fin = open(urspath, "r")
        for l in fin.readlines():
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
                    this_gene = Gene(this_id, this_urs.__len__(), urs=this_urs, has_dbd=this_has_dbd, repressor=this_repressor, pwm=this_pwm) 
                    ret_genes.append(this_gene)
        fin.close()
        
        """Now check what we build versus any command-line parameters...."""
        count_tf = 0
        count_reporter = 0
        for gid in ret_genes:
            if gid.has_dbd:
                count_tf += 1
            else:
                count_reporter += 1
        ap.params["numtr"] = count_tf
        ap.params["numreporter"] = count_reporter
        return ret_genes
    else:
        return None