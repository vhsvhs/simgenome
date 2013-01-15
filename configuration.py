try:
   import cPickle as pickle
except:
   import pickle
try:
    import copy
    from datetime import datetime
    import math
    from mpilibs import *
    import numpy
    from numpy import *
    import operator
    from operator import mul, add
    import os
    import random
    import re
    # we also import scipy.stats in the population.do_mutations method
    import sys
    import cProfile
except ImportError:
    print "\n. Sorry, I had a problem importing one (or more) Python libraries."
    exit()
#from tools import *

POPPICKLES = "POP_HISTORY"
EXPR_PLOTS = "EXPR_HISTORY"

"""The values in this configuration file will be used as defaults.
Most of these parameter values can be overrided by command-line values."""

"""URS_LEN = nucleotide length of upstream regulatory regions.
   If the mutation model allows indels, then INIT_URS_LEN will be
   the mean expected length, but the URS length will vary
   for each individual gene."""
INIT_URS_LEN = 20
MIN_URS_LEN = 10

"""What is the initial length (# sites) in a new PWM?"""
INIT_PWM_LEN = 5

ALPHABET = ["A", "C", "G", "T"]

"""How many genomes (i.e. individuals) in the population?"""
N_GENOMES = 10

"""How many transcription factors per genome?"""
N_TR = 2

"""How many reporter genes per genome?"""
N_REPORTER = 1

"""Mean mutation rate"""
MU = 0.05
ELITE_MU = 0.01
ELITE_PROPORTION = 0.3
DBD_MU = 0.5 # how mutations will be made per DBD?
PWM_MU = 0.5 # how much change occurs to a particular PWM once selected for mutation?
PWM_LEN_MU = 0.3
PWM_MU_LEN_MAX = 2
URS_LEN_MU = 0.05
URS_LEN_INDEL_MAX = 10
P2P_MU = 0.1
P2P_MU_DELTA = 1.0

NORM_SD = 0.05 # the sd value for normal distributions, created with SciPy

"""When reproduction occurs, what proportion of reproductive events are sexual? 1.0 means only sexual reproduction,
wherease 0.0 means all clonal reproduction."""
SEXUAL_RATIO = 1.0

"""How many generations to run the genetic algorithm?"""
MAX_GA_GENS = 5000
INIT_GEN = 0

"""When sampling configurations, how many i.i.d. samples should be drawn from the cumulative distribution?"""
IID_SAMPLES = 4000

"""Each random fitness goal will involve GOAL_COMPLEXITY number of expected unique 
gene expression pulses in time.  This variable is only valid if random fitness goals are being generated.
i.e., you don't specify --patternpath"""
GOAL_COMPLEXITY = 1

"""Each random fitness goal will include ideal expression pulses that are
smaller than or equal to GOAL_PULSE_MAX number of time ticks apart from each other.
This variable is only valid if random fitness goals are being generated.
i.e., you don't specify --patternpath""" 
GOAL_PULSE_MAX = 2

"""maximum distance between TFs to include gamma term"""
MAX_GD = 3

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
between TFs.  Big values facilitate cooperative binding across long distances,
whereas small values make the coopfunc dropoff quickly."""
V_RATE_OF_COOP_DECAY = 3



