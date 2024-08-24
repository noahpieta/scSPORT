#!/home/users/astar/gis/zhangyu2/.conda/envs/Env_py3/bin/python
__doc__="""
    the script is for calculate mutrate and coverage of each line (nucleotide)
    iteratively instead of make a directory

    it's similar but not same with previous get_mutrate.v1.*.py. the previous one used dictionary but
    here the iterator of each line was used. therefore, it would not able to 1) select specified gene,
    2) not all nucleotide would be reserved in the final outpufile,

    changes to release v2.2
    the previous one is to calculated gene_expr normlized coveragee
        normalized_cov = reads_count(nt)/read_count(gene)
    in the v2.2 we add a choice to normlized by total reads counts of the sample
        normalized_cov = reads_count(nt)/total_reads_count(gene, by million)

    changes in release 3.4
    1) don't output each gene mutrate as space limit

    """
__version__="v3.3"
__author__="noahzy"
__last_modify__="10-Oct-2019"

import gzip
import io
from sys import stdout
import argparse
import os
import shutil
from copy import deepcopy
from multiprocessing import Pool

def mk_idx(reffile):
    seqs={}
    genes=[]
    with open(reffile) as ref:
        for l in ref:
            if l[0] == ">":
                genes.append(l[1:].strip("\n"))
                seqs[genes[-1]]=""
            else:
                seqs[genes[-1]]+=l.strip('\n').upper()
    ref_idx={}
    for g in genes:
        ref_idx[g]=[0 for n in range(len(seqs[g]))]
    return ref_idx

def iterator_mutrate(inputfile, fcfile, out_diretory, out_prefix, win_threshold, cov_threshold=-1, gene_output_cov_threshold=10):
    ## make output file
    if out_diretory[-1] == "/":
        pass
    else:
        out_diretory = out_diretory+"/"
    if out_prefix != "-":
        output = gzip.open(out_prefix,'wt')
    else:
        output = stdout
    ## make a dictionary for read count in each gene_expr
    gexpr, ref_idx, total_reads_count = make_gene_cov(fcfile)
    total_reads_count = float(total_reads_count)/1000000## million reads
    ## make variances for the function
    #cov_map = deepcopy(ref_idx)
    #mutrate_map = deepcopy(ref_idx)
    ## input the bam-readcount file (gzipped or not)
    if inputfile.split(".")[-1] in ["gz","gzip"]:
        #myopen = gzip.open
        input = io.TextIOWrapper(io.BufferedReader(gzip.open(inputfile,'rb')))
    else:
        input = open(inputfile)
    ## iterate the bam-readcount file
    ## set up win_cov, for a x nt window for estimate the position

    #    process_line(l,ref_idx, gexpr, output,cov_threshold=1)
    ## make a process_line funtion
    #def process_line(l,ref_idx, gexpr, output,cov_threshold=1):
    def process_line(l):
        ## set nonlocal variables
        nonlocal gexpr
        nonlocal total_reads_count ## for the normalization by total reads of the sample (million reads)
        #nonlocal output
        nonlocal cov_threshold
        ## deconstruct each line
        i = l.strip('\n').split()
        gene=i[0]
        pos=int(i[1])-1
        refnt=i[2]
        coverage=int(i[3])
        ## calculate normalized coverage
        g_len, g_readcount = gexpr[gene]
        ## normalized_cov = round(float(coverage)/g_readcount,8)
        normalized_cov = round(float(coverage)/total_reads_count,8)
        ## if there is any filter for each line, just add here
        if cov_threshold == -1:
            pass
        elif coverage < cov_threshold:
            #continue
            return
        ## make detail column to show propotion of each nt(ATCGN)
        ## the format is A:numa;T:numt...
        detail = []
        ins, delete = 0, 0
        for item in i[4:]:
            ## sum number of  insertions or deleltion together
            if item[0] == "+":
                ins += int(item.split(":")[1])
            elif item[0] == "-":
                delete += int(item.split(":")[1])
            else:
                nt = item.split(":")[0]
                num = item.split(":")[1]
                detail.append(nt+":"+num)
        if ins != 0:
            nt = "+"
            num = str(ins)
            detail.append(nt+":"+num)
        if delete != 0:
            nt = "-"
            num = str(delete)
            detail.append(nt+":"+num)
        ## the mutant read counts without modification
        ## insertions and deletions are also counted
        ident=int(i[4].split(":")[1])
        mutant=coverage-ident  ## all mutants are counted,
        if coverage == 0:
            mutrate = "NA"
        else:
            mutrate = round(float(mutant)/coverage,8)
        ## put the coverage and mutrate into pseudo_matrix (map)
        #cov_map[gene][pos]=coverage
        #mutrate_map[gene][pos]=mutrate
        ## output bed-like format
        output.write("\t".join([gene,i[1],str(pos+2),gene+"."+i[1],str(mutrate),".",
                                str(coverage),str(mutant),str(normalized_cov),
                                str(g_readcount),refnt,";".join(detail)])+"\n")
        return gene,"\t".join([gene,i[1],str(pos+2),gene+"."+i[1],str(mutrate),".",
                    str(coverage),str(mutant),str(normalized_cov),
                    str(g_readcount),refnt,";".join(detail)]), coverage, normalized_cov, pos


    ## try multiprocessing failed,
    ## single threshod were used
    ## output to files for each gene which reach the threshold
    def parse_single_gene(win_len=100, normalized_cov_threshold=10,len_prop_threshold=None):
        ## notes for thresholds:
        ## normalized_cov_threshold is x reads / total reads, x: cov at position, total reads is measured by million
        ## win_threshold is the threshold for positions number which > normalized_cov_threshold
        nonlocal input
        output_temp = []
        gcov_meet_condition = []
        previous_gene = ""
        gene_lens = 0
        for l in input:
            gene, line, cov, ncov, pos = process_line(l)
            if gene == previous_gene:
                output_temp.append(line)
                gene_lens += 1
                if ncov >= normalized_cov_threshold:
                    gcov_meet_condition.append(ncov)
                #OUT.write(line)
            else:
                #OUT.write(line)
                ## measure if output the gene
                ## add a option for propotional threshold of win_len
                if len_prop_threshold == None:
                    pass
                else:
                    win_len = float(len_prop_threshold)*gene_lens
                if len(gcov_meet_condition) >= win_len:
                    ## skip output each genes as limit of space
                #    OUT = gzip.open(out_diretory+out_prefix+"."+gene+".mutrate.txt.gz",'wt')
                #    OUT.write("\n".join(output_temp)+"\n")
                #    OUT.close()
                    pass
                else:
                    pass
                ## reset varibles for the next cycle
                previous_gene = gene
                output_temp = [line]
                gcov_meet_condition = []
                gene_lens = 1
                if ncov >= normalized_cov_threshold:
                    gcov_meet_condition.append(ncov)
        ## run the last time for the last one
        ## add a option for propotional threshold of win_len
        if len_prop_threshold == None:
            pass
        else:
            win_len = float(len_prop_threshold)*gene_lens
        ## >>> last one start >>>
        if len(gcov_meet_condition) >= win_len:
        #    OUT = gzip.open(out_diretory+out_prefix+"."+gene+".mutrate.gene.gz",'wt')
        #    OUT.write("\n".join(output_temp)+"\n")
        #    OUT.close()
            pass  ## skip output each genes as limit of storage
        else:
            pass
        ## <<< last one end <<<

    # parse_single_gene(len_prop_threshold=0.5)
    if win_threshold < 1:
        parse_single_gene(normalized_cov_threshold=gene_output_cov_threshold,len_prop_threshold=win_threshold)
    else:
        parse_single_gene(win_len=win_threshold,normalized_cov_threshold=gene_output_cov_threshold)
    ## run process_line only
    #for l in input:
    #    gene, line, cov, ncov, pos = process_line(l)

def avg(test,digits=8):
    test = list(map(float,test))
    if len(test) == 0:
        return 0
    else:
        average = round(sum(test)/len(test),digits)
        return average

def make_gene_cov(fcfile):
    d = {}
    ref_idx={}
    total_reads_count = 0
    with open(fcfile) as fc:
        for l in fc:
            if l[0] == "#": continue
            i = l.strip('\n').split('\t')
            if i[0] == "Geneid":    continue
            [gene_len, gene_cov] = i[5:7]
            total_reads_count += int(gene_cov)
            d[i[1]] = (int(gene_len), int(gene_cov))
            ref_idx[i[1]]=[0 for n in range(int(gene_len))]
    return d,ref_idx,total_reads_count

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i","--input",required=True, help="input bam-readcount output file")
    #parser.add_argument("-f","--ref",required=True, help="input ref fasta file")
    parser.add_argument("-o","--outprefix",default="-")
    parser.add_argument("-d","--directory",default="./temp")
    parser.add_argument("-c","--cov_threshold",default=-1, help="cov threshold for outputing each record")
    parser.add_argument("-gc","--gene_output_cov_threshold",default=10,help="cov threshold for outputing each genes")
    parser.add_argument("-w","--win_threshold",default=50,help="win threshold for outputing each genes")
    parser.add_argument("-e","--gene_readcounts",help="input the featureCounts file")
    args = parser.parse_args()

    #ref_idx=mk_idx(args.ref)
    if os.path.exists(args.directory):
        #os.remove(args.directory)
        shutil.rmtree(args.directory, ignore_errors=True) # remove unempty directory
        #os.makedirs(args.directory)
        folder = args.directory
        folder.strip('/')
    else:
        #os.makedirs(args.directory)
        folder = args.directory
        folder.strip('/')

    #ref_idx = mk_idx(args.ref)
    iterator_mutrate(args.input, args.gene_readcounts, folder, args.outprefix, win_threshold=float(args.win_threshold),cov_threshold=float(args.cov_threshold),gene_output_cov_threshold=float(args.gene_output_cov_threshold))
