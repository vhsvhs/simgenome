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
        
    x = ap.getOptionalArg("--init_pwm_len")
    if x != False:
        ap.params["init_pwm_len"] = int(x)
    else:
        ap.params["init_pwm_len"] = INIT_PWM_LEN
        
    x = ap.getOptionalArg("--urs_len")
    if x != False:
        ap.params["init_urs_len"] = int(x)
    else:
        ap.params["init_urs_len"] = URS_LEN