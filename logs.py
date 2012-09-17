from argparser import *
from configuration import *

def log_generation(ap, line):
    fout = open(ap.params["workspace"] + "/" + ap.params["runid"] + "/LOGS/generations.txt", "a")
    fout.write(line + "\n")
    fout.close()