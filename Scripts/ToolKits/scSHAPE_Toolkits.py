## loading package 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
#from matplotlib import pylab
#import seaborn as sns 
from sklearn.metrics.pairwise import cosine_similarity
from statannot import add_stat_annotation
import matplotlib as mpl
from scipy import stats, cluster
import glob
import re
from sklearn.metrics import mean_squared_error, r2_score
from statsmodels.stats import multitest
#from matplotlib.gridspec import GridSpec
#from matplotlib_venn import venn2, venn3
#from matplotlib import pyplot as plt

import warnings 
warnings.simplefilter('ignore')



## load_large_dataFrame
def load_large_dataFrame(input_file, sep=",", header=0, names=None, index_col=0, chunksize=100000, compressed=False, comment="#"):
    if compressed:
        TextFileReader = pd.read_csv(input_file, names=names, chunksize=chunksize, sep=sep, header=header,index_col=index_col, compression='gzip', comment=comment)
    else:
        TextFileReader = pd.read_csv(input_file, names=names, chunksize=chunksize, sep=sep, header=header,index_col=index_col, comment=comment)
    dfList=[]
    for df in TextFileReader:
        dfList.append(df)
    final_df = pd.concat(dfList,sort=False)
    return final_df


#### modules for merge mutrate matrix and filter by the coverage | average coverage
def make_mutrate_of_N_size_windows(mutrate_df, win_size=10, cov_threshold=60):
    win_size = float(win_size)
    mutrate_df['pos'] = (mutrate_df['pos']//win_size)*win_size
    mutrate_df['pos'] = mutrate_df['pos'].astype('int64')
    win_mutrate_df = mutrate_df.groupby(['gene','pos'], axis=0).sum()
    win_mutrate_df['mutrate'] = np.divide(win_mutrate_df['mutant'], win_mutrate_df['coverage'])
    # filter by coverage in each window, 60 (reads) * 50 (nt) as default
    win_mutrate_df = win_mutrate_df[win_mutrate_df['coverage'] >= win_size*cov_threshold]
    win_mutrate_df = win_mutrate_df.loc[:, ['mutrate']]
    return win_mutrate_df
def make_mutant_of_N_size_windows(mutrate_df, win_size=10, cov_threshold=60):
    win_size = float(win_size)
    mutrate_df['pos'] = (mutrate_df['pos']//win_size)*win_size
    mutrate_df['pos'] = mutrate_df['pos'].astype('int64')
    win_mutrate_df = mutrate_df.groupby(['gene','pos'], axis=0).sum()
    win_mutrate_df['mutrate'] = np.divide(win_mutrate_df['mutant'], win_mutrate_df['coverage'])
    # filter by coverage in each window, 60 (reads) * 50 (nt) as default
    win_mutrate_df = win_mutrate_df[win_mutrate_df['coverage'] >= win_size*cov_threshold]
    win_mutrate_df = win_mutrate_df.loc[:, ['mutant']]
    return win_mutrate_df
def make_coverage_of_N_size_windows(mutrate_df, win_size=10, cov_threshold=60):
    win_size = float(win_size)
    mutrate_df['pos'] = (mutrate_df['pos']//win_size)*win_size
    mutrate_df['pos'] = mutrate_df['pos'].astype('int64')
    win_mutrate_df = mutrate_df.groupby(['gene','pos'], axis=0).sum()
    win_mutrate_df = win_mutrate_df[win_mutrate_df['coverage'] >= win_size*cov_threshold]
    win_mutrate_df = win_mutrate_df.loc[:, ['coverage']]
    return win_mutrate_df

## process mutrate file into matrix
def read_mutrate_file(params, win_size=10, cov_threshold=50, genes=None):
    label,sample_file = params
    print (label)
    col_names = ["gene","pos","end","id","mutrate","strand","coverage","mutant","normalized_cov","gene_expr","ntref","detail"]
    #input_df = pd.read_csv(sample_file,compression="gzip", sep="\t", names=col_names, index_col="id")
    input_df = load_large_dataFrame(sample_file, compressed=False, sep="\t", names=col_names, index_col="id")
    if genes:
        input_df['gene'] = input_df['gene'].astype('str')
        input_df = input_df.loc[input_df['gene'].isin(genes)]
    if win_size == 1:
        sample_mutrate_df = input_df.loc[:, ['gene','pos','mutrate','coverage']]
        sample_mutrate_df = sample_mutrate_df.loc[sample_mutrate_df['coverage']>=cov_threshold]
        sample_mutrate_df = sample_mutrate_df[['gene','pos','mutrate']]
        sample_mutrate_df.rename(columns={"mutrate":label},inplace=True)
        sample_mutrate_df.set_index(['gene','pos'], inplace=True)
    else:
        sample_mut_cov_df = input_df.loc[:, ['gene','pos','mutant','coverage']]
        sample_mutrate_df = make_mutrate_of_N_size_windows(sample_mut_cov_df, win_size, cov_threshold)
        sample_mutrate_df.rename(columns={"mutrate":label},inplace=True)
    return sample_mutrate_df

def read_mutant_file(params, win_size=10, cov_threshold=50, genes=None):
    label,sample_file = params
    print (label)
    col_names = ["gene","pos","end","id","mutrate","strand","coverage","mutant","normalized_cov","gene_expr","ntref","detail"]
    #input_df = pd.read_csv(sample_file,compression="gzip", sep="\t", names=col_names, index_col="id")
    input_df = load_large_dataFrame(sample_file, compressed=False, sep="\t", names=col_names, index_col="id")
    if genes:
        input_df['gene'] = input_df['gene'].astype('str')
        input_df = input_df.loc[input_df['gene'].isin(genes)]
    if win_size == 1:
        sample_mutrate_df = input_df.loc[:, ['gene','pos','mutant','coverage']]
        sample_mutrate_df = sample_mutrate_df.loc[sample_mutrate_df['coverage']>=cov_threshold]
        sample_mutrate_df = sample_mutrate_df[['gene','pos','mutant']]
        sample_mutrate_df.rename(columns={"mutant":label},inplace=True)
        sample_mutrate_df.set_index(['gene','pos'], inplace=True)
    else:
        sample_mut_cov_df = input_df.loc[:, ['gene','pos','mutant','coverage']]
        sample_mutrate_df = make_mutant_of_N_size_windows(sample_mut_cov_df, win_size, cov_threshold)
        sample_mutrate_df.rename(columns={"mutant":label},inplace=True)
    return sample_mutrate_df

def read_coverage_file(params, win_size=10, cov_threshold=50, genes=None):
    label,sample_file = params
    print (label)
    col_names = ["gene","pos","end","id","mutrate","strand","coverage","mutant","normalized_cov","gene_expr","ntref","detail"]
    #input_df = pd.read_csv(sample_file,compression="gzip", sep="\t", names=col_names, index_col="id")
    input_df = load_large_dataFrame(sample_file, compressed=False, sep="\t", names=col_names, index_col="id")
    if genes:
        input_df['gene'] = input_df['gene'].astype('str')
        input_df = input_df.loc[input_df['gene'].isin(genes)]
    if win_size == 1:
        sample_mutrate_df = input_df.loc[:, ['gene','pos','mutrate','coverage']]
        sample_mutrate_df = sample_mutrate_df.loc[sample_mutrate_df['coverage']>=cov_threshold]
        sample_coverage_df = sample_mutrate_df[['gene','pos','coverage']]
        sample_coverage_df.rename(columns={"coverage":label},inplace=True)
        sample_coverage_df.set_index(['gene','pos'], inplace=True)
    else:
        sample_mut_cov_df = input_df.loc[:, ['gene','pos','mutant','coverage']]
        sample_coverage_df = make_coverage_of_N_size_windows(sample_mut_cov_df, win_size, cov_threshold)
        sample_coverage_df.rename(columns={"coverage":label},inplace=True)
    return sample_coverage_df
#### module END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#### calculate Rsq of mutant to coverage as heterogeneity
def calculate_Rsq_of_coverage_to_mutant(mutant_df, coverage_df):
    mutant_df = mutant_df.dropna(how="all")
    coverage_df = coverage_df.dropna(how="all")
    intersect_index = mutant_df.index.intersection(coverage_df.index)
    concat_list = []
    for idx in intersect_index:
        mutant_series = mutant_df.loc[idx]
        coverage_series = coverage_df.loc[idx]
        mutant_series = mutant_series.dropna()
        coverage_series = coverage_series.dropna()
        intersect_index = mutant_series.index.intersection(coverage_series.index)
        mutant_series = mutant_series.loc[intersect_index]
        coverage_series = coverage_series.loc[intersect_index]
        
        if (mutant_series.size <= 5) or (coverage_series.size <= 5):
            continue
        # print (mutant_series)
        # print (coverage_series)
        try:
            slope, intercept, r, p, stderr = stats.linregress(coverage_series, mutant_series)
        except:
            slope, intercept, r, p, stderr = np.nan, np.nan, np.nan, np.nan, np.nan
        out_series = pd.Series([slope, intercept, r, p, stderr], index=['slope','intercept','r','p','stderr'], name=idx)
        out_series['rsq'] = np.power(out_series['r'], 2)
        #print (out_series)
        #out_series.name = pd.MultiIndex.from_tuples(idx, names=['gene','pos'])
        concat_list.append(out_series)
    out = pd.concat(concat_list, axis=1)
    out = out.T
    out.index.names = ['gene', 'pos']
    return out

def heterogeneity_by_rsquare(mutant_df, coverage_df):
    dmso_mutant = mutant_df.loc[:, mutant_df.columns.get_level_values('treatment')=="DMSO"]
    nain3_mutant = mutant_df.loc[:, mutant_df.columns.get_level_values('treatment')=="NAI_N3"]
    dmso_coverage = coverage_df.loc[:, coverage_df.columns.get_level_values('treatment')=="DMSO"]
    nain3_coverage = coverage_df.loc[:, coverage_df.columns.get_level_values('treatment')=="NAI_N3"]
    
    #dmso_rsquare = calculate_Rsq_of_coverage_to_mutant(dmso_mutant, dmso_coverage)
    nai_rsquare = calculate_Rsq_of_coverage_to_mutant(nain3_mutant, nain3_coverage)
    #return dmso_rsquare, nai_rsquare
    return nai_rsquare  
#### module END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<