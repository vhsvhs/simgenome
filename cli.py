from argparser import *
from configuration import *
from genetic_algorithm import *
from population import *
from version import *
from debug_tools import *

def read_cli(ap):
    x = ap.getOptionalArg("--growth_rate")
    if x != False:
        ap.params["growth_rate"] = float(x)
    else:
        ap.params["growth_rate"] = GROWTH_FACTOR
    x = ap.getOptionalArg("--decay_rate")
    if x != False:
        ap.params["decay_rate"] = float(x)
    else:
        ap.params["decay_rate"] = DECAY_FACTOR
    