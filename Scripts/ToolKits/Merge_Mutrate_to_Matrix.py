#!python
__doc__="""\n
    the script is used to merge mutrate dataframe of each sample together
    INPUT:
        a directory contaning all mutrate dataframes
        optional: the metasheet of samples
    OUTPUT:
        mutrate matrix
        reactivity matrix
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

def main(lib_table_file, matrix_type, outputfile=None, win_size=1, cov_threshold=60, genes=None):
    sample_lst = load_large_dataFrame(lib_table_file, sep="\t", header=[0], index_col=None)
    sample_lst = sample_lst.set_index('sample_id')
    print (sample_lst)
    
    #pattern  = re.compile(r'T\/(R\S+?)\.transcripts.mutrate.filter.gz')
    #mutrate_files = glob.glob("R*/R*transcripts.mutrate.filter.gz")
    pattern  = re.compile(r'T\/(R\S+?)\.transcripts.mutrate.filter.cov_30.gz')
    mutrate_files = glob.glob("R*/R*transcripts.mutrate.filter.cov_30.gz")
    
    concat_lst = []
    for mutrate_file in mutrate_files:
        print (mutrate_file)
        sample_id = pattern.findall(mutrate_file)[0]
        print (sample_id)
    
        if matrix_type == "mutrate":
            out_df = read_mutrate_file((sample_id,mutrate_file), win_size=win_size, cov_threshold=cov_threshold, genes=None)
        elif matrix_type == "mutant":
            out_df = read_mutant_file((sample_id,mutrate_file), win_size=win_size, cov_threshold=cov_threshold, genes=None)
        elif matrix_type == "coverage":
            out_df = read_coverage_file((sample_id,mutrate_file), win_size=win_size, cov_threshold=cov_threshold, genes=None)
        concat_lst.append(out_df)
        del out_df
        
    out_df = pd.concat(concat_lst, axis=1, join="outer")
    col_df = out_df.columns.to_frame(name="id")
    col_df = pd.concat([col_df, sample_lst[['transfer', 'treatment', 'ref']]], axis=1, join="inner")
    col_df.sort_index(inplace=True)
    out_df.sort_index(axis=1,inplace=True)
    out_df.columns = pd.MultiIndex.from_frame(col_df)
    
    if outputfile == None:
        print (out_df)
    elif "hd5" in outputfile:
        out_df.to_hdf(outputfile, key="df")
    elif "csv" in outputfile:
        out_df.to_csv(outputfile)
    else:
        print (out_df)
        
        
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-l","--inputlist",
                        required=True,
                        help="input list of sample bin files")
    parser.add_argument("-t","--type",
                        required=True,
                        default="mutrate", 
                        help="output matrix type, choosing from mutrate, mutant, coverage, reactivity")
    parser.add_argument("-c","--cov",
                        default=100, 
                        help="coverage threshold of each nucleotide")
    parser.add_argument("-w","--win",
                        default=1, 
                        help="coverage threshold of each nucleotide")
    parser.add_argument("-o","--output",
                        required=False,
                        default=None,
                        help="output results")
    args = parser.parse_args()
    main(args.inputlist, args.type, args.output, win_size=int(args.win), cov_threshold=int(args.cov))
