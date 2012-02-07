import math
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
URS_LEN = 1000

"""What is the initial length (# sites) in a new PWM?"""
INIT_PWM_LEN = 3

ALPHABET = ["A", "C", "G", "T"]

"""How many genomes (i.e. individuals) in the population?"""
N_GENOMES = 1000

"""How many transcription factors per genome?"""
N_TR = 200

"""How many reporter genes per genome?"""
N_REPORTER = 2000

"""Mean mutation rate"""
MU = 0.001

"""How many generations to run the genetic algorithm?"""
MAX_GA_GENS = 1000

"""Each random fitness goal will involve GOAL_COMPLEXITY number of expected unique 
gene expression pulses in time."""
GOAL_COMPLEXITY = 20

"""Each random fitness goal will include ideal expression pulses that are
smaller than or equal to GOAL_PULSE_MAX number of time ticks apart from each other.""" 
GOAL_PULSE_MAX = 10

"""maximum distance between TFs to include gamma term"""
MAX_GD = 80