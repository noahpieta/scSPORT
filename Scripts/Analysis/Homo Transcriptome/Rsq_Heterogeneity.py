## loading package 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from matplotlib import pylab
import seaborn as sns 
from sklearn.metrics.pairwise import cosine_similarity
from statannot import add_stat_annotation
import matplotlib as mpl
from scipy import stats, cluster
import glob
import re
from sklearn.metrics import mean_squared_error, r2_score
from statsmodels.stats import multitest
from matplotlib.gridspec import GridSpec

from matplotlib_venn import venn2, venn3
from matplotlib import pyplot as plt

import warnings 
warnings.simplefilter('ignore')

## load_large_dataFrame
def load_large_dataFrame(input_file, sep=",", header=0, index_col=0, chunksize=100000, compressed=False):
    if compressed:
        TextFileReader = pd.read_csv(input_file, chunksize=chunksize, sep=sep, header=header,index_col=index_col, compression='gzip')
    else:
        TextFileReader = pd.read_csv(input_file, chunksize=chunksize, sep=sep, header=header,index_col=index_col)
    dfList=[]
    for df in TextFileReader:
        dfList.append(df)
    final_df = pd.concat(dfList,sort=False)
    return final_df


def calculate_Rsq_of_coverage_to_mutant(mutant_df, coverage_df):
    mutant_df = mutant_df.dropna(how="all")
    coverage_df = coverage_df.dropna(how="all")
    intersect_index = mutant_df.index.intersection(coverage_df.index)
    n = 0
    concat_list = []
    for idx in intersect_index:
        mutant_series = mutant_df.loc[idx]
        coverage_series = coverage_df.loc[idx]
        mutant_series = mutant_series.dropna()
        coverage_series = coverage_series.dropna()
        intersect_index = mutant_series.index.intersection(coverage_series.index)
        mutant_series = mutant_series.loc[intersect_index]
        coverage_series = coverage_series.loc[intersect_index]
        
        if (mutant_series.size <= 10) or (coverage_series.size <= 10):
            continue
        # print (mutant_series)
        # print (coverage_series)
        try:
            slope, intercept, r, p, stderr = stats.linregress(coverage_series, mutant_series)
        except:
            slope, intercept, r, p, stderr = np.nan, np.nan, np.nan, np.nan, np.nan
        out_series = pd.Series([slope, intercept, r, p, stderr], index=['slope','intercept','r','p','stderr'], name=idx)
        #print (out_series)
        #out_series.name = pd.MultiIndex.from_tuples(idx, names=['gene','pos'])
        concat_list.append(out_series)
        
        n += 1
        if n % 100000 == 0:
            print (n)
    out = pd.concat(concat_list, axis=1)
    out = out.T
    out.index.names = ['gene', 'pos']
    return out

def heterogeneity_by_rsquare(mutant, coverage):
    #dmso_mutant = mutant.xs('dmso', level="agent", axis=1)
    nai_mutant = mutant.xs('nai-n3', level="agent", axis=1)
    #dmso_coverage = coverage.xs('dmso', level="agent", axis=1)
    nai_coverage = coverage.xs('nai-n3', level="agent", axis=1)

    #dmso_rsquare = calculate_Rsq_of_coverage_to_mutant(dmso_mutant, dmso_coverage)
    nai_rsquare = calculate_Rsq_of_coverage_to_mutant(nai_mutant, nai_coverage)
    #return dmso_rsquare, nai_rsquare
    return nai_rsquare   


## KEY Module: iterate consectutive segments by (gene, pos) in index >>>>>>>>>>>
def positions_to_continuous_segments(lst, gap=5):
    from itertools import groupby, chain
    fun = lambda x: x[1]-x[0]
    segs = []
    junctions = []
    for k, g in groupby(enumerate(lst), fun):
        #print (k)
        #print (g)
        l1 = [ j for i, j in g ]
        #print (l1)
        junctions.append(k)
        segs.append([l1[0],l1[-1]])
    last_pos = max(lst)
    junctions = junctions + [last_pos]
    #print(segs)
    return segs, junctions
def merge_adjecent_intervals(intervals, min_gap=5):
    ## merge adjecent intervals
    def iterator_to_merge_flanking_ranges(ranges, min_gap=2):
        """
        Merge overlapping and adjacent ranges and yield the merged ranges
        in order. The argument must be an iterable of pairs (start, stop).

        >>> list(merge_ranges([(5,7), (3,5), (-1,3)]))
        [(-1, 7)]
        >>> list(merge_ranges([(5,6), (3,4), (1,2)]))
        [(1, 2), (3, 4), (5, 6)]
        >>> list(merge_ranges([]))
        []
        """
        ranges = iter(sorted(ranges))
        current_start, current_stop = next(ranges)

        for start, stop in ranges:
            if start - current_stop > min_gap:
                # Gap between segments: output current segment and start a new one.
                yield current_start, current_stop
                current_start, current_stop = start, stop
            else:
                # Segments adjacent or overlapping: merge.
                current_stop = max(current_stop, stop)
        yield current_start, current_stop
    concat_list = []
    for s,e in iterator_to_merge_flanking_ranges(intervals, min_gap):
        #print (s,e)
        concat_list.append([s,e])
    return concat_list
def iterate_segments(wide_shape_reac):
    ## iterate continuous intervals
    [gene_col, pos_col] = wide_shape_reac.index.names
    idx_df = wide_shape_reac.index.to_frame()
    idx_df.index = list(range(idx_df.index.size))
    idx_df
    tx_concat_list = []
    for gene, tx_idx in idx_df.groupby([gene_col]):
        seg_concat_list = []
        seg_lst, junction_lst = positions_to_continuous_segments(tx_idx[pos_col])
        seg_lst =  merge_adjecent_intervals(seg_lst, min_gap=5)

        #print (junction_lst)
        for seg in seg_lst:
            #print (seg)
            selected_idx = pd.MultiIndex.from_frame(tx_idx.loc[tx_idx[pos_col].isin(list(range(seg[0], seg[1]+1)))])
            yield gene, "-".join(list(map(str,seg))), wide_shape_reac.loc[selected_idx]
## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

def smooth_by_avg_reac_of_rolling_wins(crude_reac, rolling_wins=3):
    concat_lst = []
    for gene, seg_id, seg_reac in iterate_segments(crude_reac):
        seg_rolling_reac = seg_reac.rolling(window=rolling_wins, min_periods=rolling_wins, center=True).mean(numeric_only=True)
        concat_lst.append(seg_rolling_reac)
    rolling_reac = pd.concat(concat_lst)
    return rolling_reac