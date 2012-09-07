from argparser import *
from configuration import *

def log_generation(ap, line):
    fout = open(ap.getArg("--workspace") + "/" + ap.getArg("--runid") + "/LOGS/generations.txt", "a")
    fout.write(line + "\n")
    fout.close()