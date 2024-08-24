#!python
__doc__="""\n
    the script is used to calculate the Rsq between mutant and coverage as 
    INPUT:
        wide mutant matrix
        wide coverage matrix
    OUTPUT:
        Rsq of each nucleotide (windows)
    """
import argparse
import glob
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import multiprocessing as mp
import matplotlib

from scSHAPE_Toolkits import *

def Calculate_Rsq(mutant_df, coverage_df):
    ref_rsquare = heterogeneity_by_rsquare(mutant_df, coverage_df)
    return ref_rsquare

def main(mutant_matrix_file, coverage_matrix_file, outputfile):
    if "hd5" in mutant_matrix_file:
        mutant_df = pd.read_hdf(mutant_matrix_file, key="df")
    elif "csv" in mutant_matrix_file:
        mutant_df = pd.read_csv(mutant_matrix_file, header=[0,1,2,3], index_col=[0,1])
        
    if "hd5" in coverage_matrix_file:
        coverage_df = pd.read_hdf(coverage_matrix_file, key="df")
    elif "csv" in coverage_matrix_file:
        coverage_df = pd.read_csv(coverage_matrix_file, header=[0,1,2,3], index_col=[0,1])
        
    
    rsq_df = Calculate_Rsq(mutant_df, coverage_df)
    #rsq_df = Calculate_Rsq(mutant_df.iloc[0:1000], coverage_df.iloc[0:1000])
    
    
    if "hd5" in outputfile:
        rsq_df.to_hdf(outputfile, key="df")
    elif "csv" in outputfile:
        rsq_df.to_csv(outputfile)
    else:
        print (rsq_df)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-m","--mutant",
                        required=True,
                        help="input mutant matrix")
    parser.add_argument("-c","--coverage",
                        required=True,
                        help="input coverage matrix")
    parser.add_argument("-o","--output",
                        required=False,
                        default="stdout",
                        help="output Rsq results")
    args = parser.parse_args()
    main(args.mutant, args.coverage, args.output)