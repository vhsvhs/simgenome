import os
import random
import re
import sys

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
N_GENOMES = 3

"""How many transcription factors per genome?"""
N_TR = 10

"""How many reporter genes per genome?"""
N_REPORTER = 10

"""Mean mutation rate"""
MU = 0.001

"""How many generations to run the genetic algorithm?"""
MAX_GA_GENS = 1000