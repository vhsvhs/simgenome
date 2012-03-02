from datetime import datetime
import math
from numpy import *
import operator
from operator import mul, add
import os
from progressbar import *
import random
import re
import sys
from tools import *

"""The values in this configuration file will be used as defaults.
Most of these parameter values can be overrided by command-line values."""

"""URS_LEN = nucleotide length of upstream regulatory regions.
   If the mutation model allows indels, then URS_LENGTH will be
   the mean expected length, but the URS length will vary
   for each individual gene."""
URS_LEN = 100

"""What is the initial length (# sites) in a new PWM?"""
INIT_PWM_LEN = 3

ALPHABET = ["A", "C", "G", "T"]

"""How many genomes (i.e. individuals) in the population?"""
N_GENOMES = 2

"""How many transcription factors per genome?"""
N_TR = 3

"""How many reporter genes per genome?"""
N_REPORTER = 2

"""Mean mutation rate"""
MU = 0.001

"""How many generations to run the genetic algorithm?"""
MAX_GA_GENS = 10

"""When sampling configurations, how many i.i.d. samples should be drawn from the cumulative distribution?"""
IID_SAMPLES = 1000

"""Each random fitness goal will involve GOAL_COMPLEXITY number of expected unique 
gene expression pulses in time."""
GOAL_COMPLEXITY = 1

"""Each random fitness goal will include ideal expression pulses that are
smaller than or equal to GOAL_PULSE_MAX number of time ticks apart from each other.""" 
GOAL_PULSE_MAX = 2

"""maximum distance between TFs to include gamma term"""
MAX_GD = 1

"""Fitness of temporal expression will be evaluated on timeslices zero through MAX_TIME."""
MAX_TIME = 5

"""Minimum distance between transcription factors bound to the same upstream regulatory sequence."""
MIN_TF_SEPARATION = 0

"""If the method get_expression returns a value greater than ACTIVATION_THRESHOLD then activation will begin on the downstream 
gene.  Else, the expression level of the downstream gene will begin to decay."""
ACTIVATION_THRESHOLD = 0.1

"""Genes are never truly 'off'.  The minimum expression level is MINIMUM_ACTIVITY_LEVEL."""
MINIMUM_ACTIVITY_LEVEL = 0.001
MAXIMUM_ACTIVITY_LEVEL = 1.0

"""If RNA polymerase is active, then the expression level at time t+1 will be GROWTH_FACTOR times the expr. level at current time t."""
GROWTH_FACTOR = 1.3
DECAY_FACTOR = 1.2


"""This rate gets used in the function coopfunc (in Landscape) to
control the rate of decay of the cooperative or competitive interactions
between TFs."""
V_RATE_OF_COOP_DECAY = 80

"""The scalar used to convert information to Kd.... k = 1/(e^(aI)), where a is alpha and I is the binding strength (in bits)"""
INFO_ALPHA = 0.01

"""Inversely proportion to URS length."""
DELTA_G_SCALAR = 0.1
