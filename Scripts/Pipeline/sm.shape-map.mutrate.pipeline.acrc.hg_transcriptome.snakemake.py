## __doc__="""
##    this sctrip is using for Shape-Map which seem no difference between
##    SC or BULK
##
##    the update of version 1.5
##    cutadapt twice to remove the internal adapter which ligated which reverse transcript
##
##    the change in version 2.5:
##    1) change to zy made reference
##    2) bowtie2 used two params: -k 20, -R 4
##    """
##
## __author__="noahzy"
## __version__="3.0"
## __last_modify__="17-Aug-2022"

shell.executable("/bin/bash")
## unofficial bash strict mode
shell.prefix("source ~/.bashrc; set -euo pipefail;")

#ADAPTER = {"a":"TGTCTCTTATACACATCT", "A":"TGTCTCTTATACACATCT"}
ADAPTER = {"a":"CCGAGCCCACGAGAC", "A":"GACGCTGCCGACGA"} ## i5, i7
INTERADAPTER = {"g":"AAGCAGTGGTATCAACGCAGAGTACATGGG","a":"GTACTCTGCGTTGATACCAGTGCTT","G":"AAGCAGTGGTATCAACGCAGAGTACATGGG","A":"GTACTCTGCGTTGATACCAGTGCTT"}

## Riboswithes and rRNAs
RiboSwitch_WT_REF = "/home/users/astar/gis/zhangyu2/scratch/Projects/Genome_Reference/RiboSwitches_Ref/RiboSwitches_WT/RiboSwitch_Ref_WT.bowtie2/RiboSwitch_Ref_WT"
RiboSwitch_MT_REF = "/home/users/astar/gis/zhangyu2/scratch/Projects/Genome_Reference/RiboSwitches_Ref/RiboSwitches_MT/RiboSwitch_Ref_MT.bowtie2/RiboSwitch_Ref_MT"
Riboswitch_WT_FASTA = "/home/users/astar/gis/zhangyu2/scratch/Projects/Genome_Reference/RiboSwitches_Ref/RiboSwitches_WT/RiboSwitch_Ref_WT.fa"
Riboswitch_MT_FASTA = "/home/users/astar/gis/zhangyu2/scratch/Projects/Genome_Reference/RiboSwitches_Ref/RiboSwitches_MT/RiboSwitch_Ref_MT.fa"
Riboswitch_WT_SAF = "/home/users/astar/gis/zhangyu2/scratch/Projects/Genome_Reference/RiboSwitches_Ref/RiboSwitches_WT/RiboSwitch_Ref_WT.saf"
Riboswitch_MT_SAF = "/home/users/astar/gis/zhangyu2/scratch/Projects/Genome_Reference/RiboSwitches_Ref/RiboSwitches_MT/RiboSwitch_Ref_MT.saf"
RiboSwitches_WT_BED = "/home/users/astar/gis/zhangyu2/scratch/Projects/Genome_Reference/RiboSwitches_Ref/RiboSwitches_WT/RiboSwitch_Ref_WT.bed"
RiboSwitches_MT_BED = "/home/users/astar/gis/zhangyu2/scratch/Projects/Genome_Reference/RiboSwitches_Ref/RiboSwitches_MT/RiboSwitch_Ref_MT.bed"

## Reference
HOMO_BOWTIE2_INDEX="/home/users/astar/gis/zhangyu2/scratch-LORNAGS/zhangyu2/Genome_Ref/Gencode_GRCh38.longest_transcriptome.modified/gencode.v36.longestCDS.rRNA.fa"
HOMO_BOWTIE2_FA="/home/users/astar/gis/zhangyu2/scratch-LORNAGS/zhangyu2/Genome_Ref/Gencode_GRCh38.longest_transcriptome.modified/gencode.v36.longestCDS.rRNA.fa"
HOMO_SAF = "/home/users/astar/gis/zhangyu2/scratch-LORNAGS/zhangyu2/Genome_Ref/Gencode_GRCh38.longest_transcriptome.modified/gencode.v36.longestCDS.rRNA.saf"
HOMO_BED = "/home/users/astar/gis/zhangyu2/scratch-LORNAGS/zhangyu2/Genome_Ref/Gencode_GRCh38.longest_transcriptome.modified/gencode.v36.longestCDS.rRNA.bed"

## Tools
CUTADAPT = "/home/users/astar/gis/zhangyu2/.conda/envs/Env_py3/bin/cutadapt"
BOWTIE2 = "/home/users/astar/gis/zhangyu2/.conda/envs/Env_py3/bin/bowtie2"
BEDTOOLS = "/home/users/astar/gis/zhangyu2/.conda/envs/Env_py3/bin/bedtools"
SAMTOOLS = "/home/users/astar/gis/zhangyu2/local_bin/samtools"
FEATURECOUNT = "/home/users/astar/gis/zhangyu2/.conda/envs/Env_py3/bin/featureCounts"
BAMREADCOUNT = "/home/users/astar/gis/zhangyu2/local_bin/bam-readcount"
GETMUTRATE = "/home/users/astar/gis/zhangyu2/scratch/Projects/scSHAPE_Map/scSHAPE_scripts/get_mutrate.v3.5.py"
FILTERMUTRATE = "/home/users/astar/gis/zhangyu2/scratch/Projects/scSHAPE_Map/scSHAPE_scripts/filter_mutrate.py"

def get_id():
    import os
    import glob
    fqi = glob.glob("merged*fq.gz")[0]
    id = os.path.abspath(fqi).split('/')[-2]
    return id
ID = get_id()
print (ID)
if "WT" in ID:
    BOWTIEIDX = RiboSwitch_WT_REF
    ANNOTATED_SAF = Riboswitch_WT_SAF
    TX_EBD = RiboSwitches_WT_BED
    REFFA = Riboswitch_WT_FASTA
elif "MT" in ID:
    BOWTIEIDX = RiboSwitch_MT_REF
    ANNOTATED_SAF = Riboswitch_MT_SAF
    TX_EBD = RiboSwitches_MT_BED
    REFFA = Riboswitch_MT_FASTA
ID = ID.split("_")[0]
print (ID)

def get_fq(r):
    import glob
	#print (glob.glob("R*_R{}_001.fastq.gz".format(r)))
    #return glob.glob("R*_R{}.fastq.gz".format(r))
    #print (glob.glob("merged*{}.fq.gz".format(r)))
    return glob.glob("merged.R*.{}.fq.gz".format(r))

rule all:
    input:
        bai1 = ID+".ribosxitch.bam.bai",
        bai2 = ID+".ribosxitch.md.bam.bai",
        bai3 = ID+".transcripts.bam.bai",
        bai4 = ID+".transcripts.md.bam.bai",
        outbed_ribosxitch = ID+".ribosxitch.mutrate.filter_pos.merge.bed",
        outbed_transcripts = ID+".transcripts.mutrate.filter_pos.merge.bed",
        #out_gzip = ID+".mutrate.txt.gz",
        #filtered_mutrate = ID+".mutrate.filter.gz",

rule cutadapt:
    input:
        fqi = get_fq(1),
        fqj = get_fq(2),
    output:
        trim_fqi = temp("{ID}.R1.cutadapt.fastq.gz"),
        trim_fqj = temp("{ID}.R2.cutadapt.fastq.gz"),
    params: min_len=30
    shell:"""
        {CUTADAPT} -m {params.min_len} -a {ADAPTER[a]} -A {ADAPTER[A]} \
        -o {output.trim_fqi} -p {output.trim_fqj} {input.fqi} {input.fqj}
        """

rule cutadapt2:
    input:
        fqi = rules.cutadapt.output.trim_fqi,
        fqj = rules.cutadapt.output.trim_fqj,
    output:
        trim_fqi = "{ID}.R1.cutadapt2.fastq.gz",
        trim_fqj = "{ID}.R2.cutadapt2.fastq.gz",
    params: min_len=30
    shell:"""
        {CUTADAPT} -m {params.min_len} -g {INTERADAPTER[g]} -a {INTERADAPTER[a]} -G {INTERADAPTER[G]} -A {INTERADAPTER[A]} \
        -o {output.trim_fqi} -p {output.trim_fqj} {input.fqi} {input.fqj}
        """

rule bowtie_map_to_ribosxitch:
    input:
        fqi = rules.cutadapt2.output.trim_fqi,
        fqj = rules.cutadapt2.output.trim_fqj,
        bowtie_index = BOWTIEIDX
    output:
        bam="{ID}.ribosxitch.bam",
        unmapped_prefix = "{ID}_non_ribosxitch",
    params: k=20, R=4,
    shell:"""
        {BOWTIE2} --quiet -k {params.k} -R {params.R} -p 4 -x {input.bowtie_index} \
        -1 {input.fqi} -2 {input.fqj} --un-conc-gz {output.unmapped_prefix} | \
        {SAMTOOLS} view -@ 4 -bS | {SAMTOOLS} sort -@ 4 -o {output.bam} &&
        ls {output.unmapped_prefix}* > {output.unmapped_prefix}
        """

rule bam_idx_ribosxitch:
    input:
        bam1 = rules.bowtie_map_to_ribosxitch.output.bam,
    output:
        bai1 = "{ID}.ribosxitch.bam.bai",
    shell:"""
        {SAMTOOLS} index {input.bam1} {output.bai1};
        """

rule featurecounts_ribosxitch:
    input:
        bam = rules.bowtie_map_to_ribosxitch.output.bam,
        annotation=ANNOTATED_SAF,
    output:
        fcout= "{ID}.ribosxtich.featurecount"
    params:
        T=1,
        F='SAF',
    shell:"""
        {FEATURECOUNT} -M -d 20 -D 100000 -T {params.T} -F {params.F} \
        -a {input.annotation} -o {output.fcout} {input.bam}
    """

rule calculate_md_ribosxitch:
    input:
        bam = rules.bowtie_map_to_ribosxitch.output.bam,
        reffa=REFFA
    output:
        mdbam = temp("{ID}.ribosxitch.md.bam")
    params:
        T=4
    shell:"""
        {SAMTOOLS} fillmd -@ {params.T} -euArE {input.bam} {input.reffa} > {output.mdbam}
        """

rule bam_idx_2_ribosxitch:
    input:
        bam1 = rules.calculate_md_ribosxitch.output.mdbam,
    output:
        bai1 = "{ID}.ribosxitch.md.bam.bai",
    shell:"""
        {SAMTOOLS} index {input.bam1} {output.bai1};
        """

rule bam_readcount_ribosxitch:
    input:
        bam = rules.calculate_md_ribosxitch.output.mdbam,
        reffa=REFFA,
        txbed=TX_EBD,
        bai = rules.bam_idx_2_ribosxitch.output.bai1
    output:
        pos_count = "{ID}.ribosxitch.bam_readcount.txt.gz"
    shell:"""
        {BAMREADCOUNT} -w 0 -f {input.reffa} -l {input.txbed} {input.bam} | gzip -c > {output.pos_count}
        """

rule calculate_mutrate_ribosxitch:
    input:
        pos_count = rules.bam_readcount_ribosxitch.output.pos_count,
        fc = rules.featurecounts_ribosxitch.output.fcout
    output:
        out_gzip = "{ID}.ribosxitch.mutrate.txt.gz",
    params:
        win_len = 0, ## window length must be > 50% of the gene length
        gene_output_cov_threshold = 10 ## coverage threshold is 10/1M reads
    shell:"""
        {GETMUTRATE} -i {input.pos_count} -e  {input.fc} -o {output.out_gzip} \
        -d each_gene -w {params.win_len} -gc {params.gene_output_cov_threshold}
        """

rule filter_by_cov_ribosxitch:
    input:
        mutrate = rules.calculate_mutrate_ribosxitch.output.out_gzip,
    output:
        filtered_mutrate = "{ID}.ribosxitch.mutrate.filter.gz",
    params:
        cov_threshold=100,
        type="ncov"
    shell:"""
        {FILTERMUTRATE} -i {input.mutrate} -o {output.filtered_mutrate} -c {params.cov_threshold} -t {params.type}
        """

rule merge_to_filter_pos_bed_ribosxitch:
    input:
        filter_mutrate = rules.filter_by_cov_ribosxitch.output.filtered_mutrate
    output:
        filter_pos_bed = "{ID}.ribosxitch.mutrate.filter_pos.merge.bed"
    shell:"""
        {BEDTOOLS} merge -i {input.filter_mutrate} > {output.filter_pos_bed}
        """

## Transcriptome

rule change_fastq_gzip:
    input:
        fq_prefix = rules.bowtie_map_to_ribosxitch.output.unmapped_prefix,
    output:
        fqi = rules.bowtie_map_to_ribosxitch.output.unmapped_prefix + ".1.fq.gz",
        fqj = rules.bowtie_map_to_ribosxitch.output.unmapped_prefix + ".2.fq.gz",
    shell:
        """
        mv {input.fq_prefix}.1 {output.fqi} &&
        mv {input.fq_prefix}.2 {output.fqj}
        """

rule bowtie_map_to_transcripts:
    input:
        fqi = rules.change_fastq_gzip.output.fqi,
        fqj = rules.change_fastq_gzip.output.fqj,
        bowtie_index = HOMO_BOWTIE2_INDEX
    output:
        bam="{ID}.transcripts.bam",
    params: k=20, R=4,
    shell:"""
        {BOWTIE2} --quiet -k {params.k} -R {params.R} -p 4 -x {input.bowtie_index} \
        -1 {input.fqi} -2 {input.fqj} | \
        {SAMTOOLS} view -@ 2 -bS | {SAMTOOLS} sort -@ 2 -o {output.bam}
        """

rule bam_idx_transcripts:
    input:
        bam1 = rules.bowtie_map_to_transcripts.output.bam,
    output:
        bai1 = "{ID}.transcripts.bam.bai",
    shell:"""
        {SAMTOOLS} index {input.bam1} {output.bai1};
        """

rule featurecounts_transcripts:
    input:
        bam = rules.bowtie_map_to_transcripts.output.bam,
        annotation=HOMO_SAF,
    output:
        fcout= "{ID}.transcripts.featurecount"
    params:
        T=1,
        F='SAF',
    shell:"""
        {FEATURECOUNT} -M -d 20 -D 100000 -T {params.T} -F {params.F} \
        -a {input.annotation} -o {output.fcout} {input.bam}
    """

rule calculate_md_transcripts:
    input:
        bam = rules.bowtie_map_to_transcripts.output.bam,
        reffa=HOMO_BOWTIE2_FA
    output:
        mdbam = temp("{ID}.transcripts.md.bam")
    params:
        T=2
    shell:"""
        {SAMTOOLS} fillmd -@ {params.T} -euArE {input.bam} {input.reffa} > {output.mdbam}
        """

rule bam_idx_2_transcripts:
    input:
        bam1 = rules.calculate_md_transcripts.output.mdbam,
    output:
        bai1 = "{ID}.transcripts.md.bam.bai",
    shell:"""
        {SAMTOOLS} index {input.bam1} {output.bai1};
        """

rule bam_readcount_transcripts:
    input:
        bam = rules.calculate_md_transcripts.output.mdbam,
        reffa=HOMO_BOWTIE2_FA,
        txbed=HOMO_BED,
        bai = rules.bam_idx_2_transcripts.output.bai1
    output:
        pos_count = "{ID}.transcripts.bam_readcount.txt.gz"
    shell:"""
        {BAMREADCOUNT} -w 0 -f {input.reffa} -l {input.txbed} {input.bam} | gzip -c > {output.pos_count}
        """

rule calculate_mutrate_transcripts:
    input:
        pos_count = rules.bam_readcount_transcripts.output.pos_count,
        fc = rules.featurecounts_transcripts.output.fcout
    output:
        out_gzip = "{ID}.transcripts.mutrate.txt.gz",
    params:
        win_len = 0, ## window length must be > 50% of the gene length
        gene_output_cov_threshold = 10 ## coverage threshold is 10/1M reads
    shell:"""
        {GETMUTRATE} -i {input.pos_count} -e  {input.fc} -o {output.out_gzip} \
        -d each_gene -w {params.win_len} -gc {params.gene_output_cov_threshold}
        """

rule filter_by_cov_transcripts:
    input:
        mutrate = rules.calculate_mutrate_transcripts.output.out_gzip,
    output:
        filtered_mutrate = "{ID}.transcripts.mutrate.filter.gz",
    params:
        cov_threshold=100,
        type="ncov"
    shell:"""
        {FILTERMUTRATE} -i {input.mutrate} -o {output.filtered_mutrate} -c {params.cov_threshold} -t {params.type}
        """

rule merge_to_filter_pos_bed:
    input:
        filter_mutrate = rules.filter_by_cov_transcripts.output.filtered_mutrate
    output:
        filter_pos_bed = "{ID}.transcripts.mutrate.filter_pos.merge.bed"
    shell:"""
        {BEDTOOLS} merge -i {input.filter_mutrate} > {output.filter_pos_bed}
        """
