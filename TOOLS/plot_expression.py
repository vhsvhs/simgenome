"""USAGE: python plot_expression FOLDERPATH
...where FOLDERPATH is the EXPR_HISTORY folder for a project"""

import sys, os

fpath = sys.argv[1]

files = os.listdir(fpath)
for f in files:
    if f.__contains__("cran"):
        os.system("r --no-save < " + fpath + "/" + f)