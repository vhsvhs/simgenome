try:
   import cPickle as pickle
except:
   import pickle
from datetime import datetime
import math
from mpilibs import *
from numpy import *
import operator
from operator import mul, add
import os
from progressbar import *
import random
import re
import sys
#from tools import *

POPPICKLES = "POP_HISTORY"
EXPR_PLOTS = "EXPR_HISTORY"

"""The values in this configuration file will be used as defaults.
Most of these parameter values can be overrided by command-line values."""

"""URS_LEN = nucleotide length of upstream regulatory regions.
   If the mutation model allows indels, then URS_LENGTH will be
   the mean expected length, but the URS length will vary
   for each individual gene."""
URS_LEN = 20

"""What is the initial length (# sites) in a new PWM?"""
INIT_PWM_LEN = 3

ALPHABET = ["A", "C", "G", "T"]

"""How many genomes (i.e. individuals) in the population?"""
N_GENOMES = 10

"""How many transcription factors per genome?"""
N_TR = 2

"""How many reporter genes per genome?"""
N_REPORTER = 1

"""Mean mutation rate"""
MU = 0.1
ELITE_MU = 0.01
ELITE_PROPORTION = 0.5
PWM_MU = MU
"""How many generations to run the genetic algorithm?"""
MAX_GA_GENS = 5000

"""When sampling configurations, how many i.i.d. samples should be drawn from the cumulative distribution?"""
IID_SAMPLES = 2000

"""Each random fitness goal will involve GOAL_COMPLEXITY number of expected unique 
gene expression pulses in time."""
GOAL_COMPLEXITY = 1

"""Each random fitness goal will include ideal expression pulses that are
smaller than or equal to GOAL_PULSE_MAX number of time ticks apart from each other.""" 
GOAL_PULSE_MAX = 2

"""maximum distance between TFs to include gamma term"""
MAX_GD = 1

"""Fitness of temporal expression will be evaluated on timeslices zero through MAX_TIME."""
MAX_TIME = 3

"""Minimum distance between transcription factors bound to the same upstream regulatory sequence."""
MIN_TF_SEPARATION = 0

"""Inversely proportion to URS length. Higher values make k_act higher, which makes it easer for genes to be expressed."""
PE_SCALAR = 1.0

"""Genes are never truly 'off'.  The minimum expression level is MINIMUM_ACTIVITY_LEVEL."""
MINIMUM_ACTIVITY_LEVEL = 0.001
MAXIMUM_ACTIVITY_LEVEL = 1.0

"""Larger values make the fitness function sharper around optimal expression, smaller values spread the function's hills further out over poor expression."""
FITNESS_SCALAR = -3.0

"""If RNA polymerase is active, then the expression level at time t+1 will be GROWTH_FACTOR times the expr. level at current time t."""
MAX_TRANSCRIPTION_RATE = 2.2
MAX_DECAY_RATE = 1.4

"""This rate gets used in the function coopfunc (in Landscape) to
control the rate of decay of the cooperative or competitive interactions
between TFs."""
V_RATE_OF_COOP_DECAY = 80

#"""The scalar used to convert information to Kd.... k = 1/(e^(aI)), where a is alpha and I is the binding strength (in bits)"""
#INFO_ALPHA = 0.01


