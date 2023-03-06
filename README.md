# scSPORT
*Copyright (c) 2023 Noah Yu Zhang*. This project is licensed under the terms of the MIT license.

__scSPORT__: a Single-cell RNA Structure Probing Method

## Introduction

RNA structure is critical for multiple steps in gene regulation. However, how structures of transcripts differ both within and between individual cells is unknown. Here, we develop __Single Cell Structure Probing Of RNA Transcripts__ (__scSPORT__), a method enabling the simultaneous determination of transcript structure and abundance at single-cell resolution, that opens the door to the systematic characterization of RNA structure-function relationships at single-cell resolution. 

## Usage
### Input:
The single-cell structure probing data were sequenced with illumina sequencing platform. The reads from each cell were separated by the barcodes. The for each cell, the read file was put into a individual folder. The hierarchical tree is as below:

/-----MUX0001  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-----MUX0001.R1.fastq   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-----MUX0001.R2.fastq  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|----MUX0002  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-----MUX0002.R1.fastq  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-----MUX0002.R2.fastq  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|----MUX0003   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-----MUX0003.R1.fastq  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-----MUX0003.R2.fastq  
......
&nbsp;&nbsp;&nbsp;&nbsp;
### Pipeline for each cell
The scSPORT used the snakemake script to cut the adapter, map to transcriptome reference and calculate mutant rate of each nucleotide for each single cell. You may change the reference and the path of dependant softwares according to your env. we suggest you to build a env with conda.   
&nbsp;&nbsp;&nbsp;&nbsp;  
The snakemake scripts are under the PATH ```scSPORTs_New/Scripts/Pipeline```  
>+ *```sm.shape-map.mutrate.pipeline.acrc.RiboSxitch.snakemake``` is used for the wildtype and mutated well-known riboswitches() which are dopped in;*  
>+ *```sm.shape-map.mutrate.pipeline.acrc.hg_transcriptome.snakemake``` is used for structural probing of human transcriptome in cells of neural differentiation.*  

After running the pipeline of each cell, the output is a matrix for each single cell containing the mutant rate of each nucleotide or 10nt-windown. Then the ```Merge_Mutrate_to_Matrix.py``` was used to merge mutant rates from cells together into a matrix. 
&nbsp;&nbsp;&nbsp;&nbsp;
### Calculate Reactivity and Heterogeneity
We used jupyter lab to do the next analysis. The jupyter notebooks are under path ```Scripts/Analysis/```. The jupyter notebooks in ```Homo Transcriptome``` folder are for cells from human neural dfferentiation; The jupyter notebooks in ```RiboSwitch``` are for the dop-in benchmarks. Here we used the Riboswitches as example:
&nbsp;&nbsp;&nbsp;&nbsp;  
<img src="/Figures/Reactivity_Heatmaps_of_Riboswitches/reac_heatmap.gNorm.Tetrahymena.svg"  width="1200" height="600"><img src="/Figures/Reactivity_Heatmaps_of_Riboswitches/reac_heatmap.gNorm.Tetrahymena.svg"  width="1200" height="600">
## Authors

*Dr. Noah Yu Zhang (zhang_yu@gis.a-star.edu.sg)*  
*Dr. Tong Zhang (zhangtong516@gmail.com)*  
*from Genome Institute of Singapore(GIS), Agency of Science, Techonology and Research(A\*STAR)*   
*under the supervision of Dr. Yue Wan (wany@gis.a-star.edu.sg)*


## License

This program is free software, and can be redistribute and/or modified under the terms of the MIT License.

## Documentation

For any information, please refer to the documentation: <http://kkkkk/>
