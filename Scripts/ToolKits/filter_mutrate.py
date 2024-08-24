#!/home/users/astar/gis/zhangyu2/.conda/envs/Env_py3/bin/python
__doc__="""
    the script is used for filter the mutrate file with coverage threshold
    criteria:
        1) absolute coverage > 100 reads in one position
        2) <optional> filter by normalized coverage ( eg: 10 / 1M )
    input:
        mutrate.txt.gz
    """
__version__="v1.0"
__author__="noahzy"
__last_modify__="10-oct-2019"

import gzip
import io
from sys import stdout
import argparse
import re

def filter_mutrate_by_cov(inputfile,outputfile,cov_threshold,filter_type):
    if inputfile.split(".")[-1] in ["gz","gzip"]:
        input = io.TextIOWrapper(io.BufferedReader(gzip.open(inputfile,'rb')))
    else:
        input = open(inputfile)
    if outputfile == "-":
        output = stdout
    else:
        output = gzip.open(outputfile,'wt')
    ## choose the thresholds for cov or normalized cov
    if filter_type == "ncov":
        idx = "normalized_cov"
    else:
        idx = "coverage"
    ## filter out each line
    names = ["gene","pos","end","id","mutrate","strand","coverage","mutant","normalized_cov","ntref","detail"]
    for l in input:
        i = l.strip('\n').split('\t')
        cols = dict(zip(names,i))
        test = float(cols[idx])
        if test >= cov_threshold:
            #print (cols[idx],cov_threshold)
            output.write(l)
    output.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i","--input",help="inputfile of mutrate.txt.gz file")
    parser.add_argument("-o","--output",default="-",help="outputfile of mutrate.filter.gz file")
    parser.add_argument("-c","--cov_threshold",default=100,help="the threshold of coverage, > c reads")
    parser.add_argument("-t","--filter_type",default="cov", help="use cov or normlized cov, default: cov")
    args = parser.parse_args()
    cov_threshold = int(args.cov_threshold)

    filter_mutrate_by_cov(args.input, args.output, cov_threshold, args.filter_type)
