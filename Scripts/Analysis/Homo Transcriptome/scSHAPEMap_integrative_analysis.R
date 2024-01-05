library(data.table)
library(stringr)
library(scales) 
library(Rfast)
library(ggplot2)
library(MOFA2)
## give path to python executable 
reticulate::use_python(python = "/Users/zhangtong/Documents/002_software/anaconda3/bin/python")

args = commandArgs(trailingOnly=TRUE)
cell_metadata_file = args[1]
heterogeneity_data_file = args[2]
expression_data_file = args[3]
alternative_splicing_data_file = args[4]

## usage: 
## Rscript scSHAPEMap_integrative_analysis.R <cell metadata> <gene heterogeneity matrix> <expression matrix> <alternative splicing PSI score matrix>


myCol = hue_pal()(4)
stageList = c("h9","day0", "day1", "day7")
names(myCol) =stageList

## set current directory as workding directory 
baseDir = getwd()
setwd(baseDir)
## redirect all output to the folder
outDir = str_c(baseDir, "/integrative_analysis/")
if(!dir.exists(outDir)) {dir.create(outDir)}

## load cell meta 
cells_meta = read.table(cell_metadata_file,sep="\t")
cells = as.character(rownames(cells_meta))

## loading structural heterogeneity per gene
gene_distance_to_pseudo = fread(heterogeneity_data_file)
names(gene_distance_to_pseudo)= c("gene", "cell", "d")
gene_distance_to_pseudo = gene_distance_to_pseudo[!is.na(d)]

gene_distance_to_pseudo_wide = dcast(gene_distance_to_pseudo, gene~ cell, value.var = "d")
dim(gene_distance_to_pseudo_wide)
dim(na.omit(gene_distance_to_pseudo_wide))
gene_distance_to_pseudo_wide_nona = na.omit(gene_distance_to_pseudo_wide)
hetero = as.data.frame(na.omit(gene_distance_to_pseudo_wide_nona[,-1]))
rownames(hetero) = str_c("hetero_",gene_distance_to_pseudo_wide_nona$gene)
#### qualtile normalization
hetero_qnorm = as.data.frame(preprocessCore::normalize.quantiles(as.matrix(hetero)))
colnames(hetero_qnorm) = colnames(hetero)
rownames(hetero_qnorm) = rownames(hetero)

## loading expression matrix: pre-compiled data by SCTransform from Seurat 
exp_dt = fread(expression_data_file)
exp_all = as.data.frame(exp_dt[,-c("gene")])
rownames(exp_all) = str_c("exp_", exp_dt$gene)
old_names = colnames(exp_all)
setkey(cell_cycle, cell)
new_names = cell_cycle[old_names, new_cell]
colnames(exp_all) = new_names
# #### qualtile normalization
# exp_all_qnorm = as.data.frame(preprocessCore::normalize.quantiles(as.matrix(exp_all)))
# colnames(exp_all_qnorm) = colnames(exp_all)
# rownames(exp_all_qnorm) = rownames(exp_all)
exp_vars = data.table("gene" = rownames(exp_all), "var" = rowVars(as.matrix(exp_all)))
#### select 1000 most variable genes for analysis: 
most_variable_genes = exp_vars[order(var, decreasing = T),][1:1000,]$gene
exp=exp_all[most_variable_genes,colnames(hetero)]

## loading Alternative splicing 
as_dt = fread(alternative_splicing_data_file)
as_all = as.data.frame(as_dt[,-c("txn")])
rownames(as_all) = str_c("AS_", as_dt$txn)
#### qualtile normalization
as = as.data.frame(preprocessCore::normalize.quantiles(as.matrix(as_all)))
colnames(as) = colnames(as_all)
rownames(as) = rownames(as_all)
as= as[,colnames(hetero)]


############################################
## MOFA+ : build  multiple models
############################################
exp_hetero_as_model_data= list("hetero" = as.matrix(hetero),"exp" = as.matrix(exp), "as" = as.matrix(as))
exp_as_model_data = list("exp" = as.matrix(exp), "as" = as.matrix(as))
exp_hetero_model_data = list("hetero" = as.matrix(hetero),"exp" = as.matrix(exp))
exp_model_data = list("exp" = as.matrix(exp)) 
#
exp_hetero_as_model = create_mofa(exp_hetero_as_model_data)
exp_as_model = create_mofa(exp_as_model_data)
exp_hetero_model = create_mofa(exp_hetero_model_data)
exp_model = create_mofa(exp_model_data)
##### define parameters for model training 
### define parameters:
get_parameters = function(MOFAobject, num_factors=6 ) {
    data_opts <- get_default_data_options(MOFAobject)
    model_opts <- get_default_model_options(MOFAobject)
    model_opts$num_factors = num_factors
    train_opts <- get_default_training_options(MOFAobject)
    # set stochastic to TRUE
    train_opts$stochastic <- T
    train_opts$gpu_mode <- T
    stochastic_opts <- get_default_stochastic_options(MOFAobject)
    return(list(
        "data_opts" = data_opts,
        "model_opts" = model_opts,
        "train_opts" = train_opts,
        "stochastic_opts" = stochastic_opts
    ))
    
} 
## 
exp_hetero_as_model_parameters = get_parameters(exp_hetero_as_model)
exp_hetero_as_model <- prepare_mofa(exp_hetero_as_model,
                        #    data_options = exp_hetero_as_model_parameters$data_opts,
                        #    training_options = exp_hetero_as_model_parameters$train_opts,
                           model_options = exp_hetero_as_model_parameters$model_opts)
                        #    stochastic_options = exp_hetero_as_model_parameters$stochastic_opt)

exp_as_model_parameters = get_parameters(exp_as_model)
exp_as_model <- prepare_mofa(exp_as_model,
                        #    data_options = exp_as_model_parameters$data_opts,
                        #    training_options = exp_as_model_parameters$train_opts,
                           model_options = exp_as_model_parameters$model_opts)
                        #    stochastic_options = exp_as_model_parameters$stochastic_opt)

exp_hetero_model_parameters = get_parameters(exp_hetero_model)
exp_hetero_model <- prepare_mofa(exp_hetero_model,
                        #    data_options = exp_hetero_model_parameters$data_opts,
                        #    training_options = exp_hetero_model_parameters$train_opts,
                           model_options = exp_hetero_model_parameters$model_opts)
                        #    stochastic_options = exp_hetero_model_parameters$stochastic_opt)

exp_model_parameters = get_parameters(exp_model)
exp_model <- prepare_mofa(exp_model,
                        #    data_options = exp_model_parameters$data_opts,
                        #    training_options = exp_model_parameters$train_opts,
                           model_options = exp_model_parameters$model_opts)
                        #    stochastic_options = exp_model_parameters$stochastic_opt)

##

if (file.exists("MOFA2_model_exp_hetero_as_stochasticVariantionalInference.hdf5")) {
    exp_hetero_as_model.svi = load_model("MOFA2_model_exp_hetero_as_stochasticVariantionalInference.hdf5")
} else {
    exp_hetero_as_model.svi  <- run_mofa(exp_hetero_as_model, file.path("MOFA2_model_exp_hetero_as_stochasticVariantionalInference.hdf5"))
}

if(file.exists("MOFA2_model_exp_as_stochasticVariantionalInference.hdf5")) {
    exp_as_model.svi = load_model("MOFA2_model_exp_as_stochasticVariantionalInference.hdf5")
} else {
    exp_as_model.svi  <- run_mofa(exp_as_model, file.path("MOFA2_model_exp_as_stochasticVariantionalInference.hdf5"))
}

if(file.exists("MOFA2_model_exp_hetero_stochasticVariantionalInference.hdf5")) {
    exp_hetero_model.svi = load_model("MOFA2_model_exp_hetero_stochasticVariantionalInference.hdf5")
} else {
    exp_hetero_model.svi  <- run_mofa(exp_hetero_model, file.path("MOFA2_model_exp_hetero_stochasticVariantionalInference.hdf5"))
} 

if(file.exists("MOFA2_model_exp_stochasticVariantionalInference.hdf5")) {
    exp_model.svi = load_model("MOFA2_model_exp_stochasticVariantionalInference.hdf5")
} else {
    exp_model.svi  <- run_mofa(exp_model, file.path("MOFA2_model_exp_stochasticVariantionalInference.hdf5"))
} 

## add sample metadata 
Nsamples = sum(exp_hetero_as_model@dimensions$N)
sample_metadata <- cells_meta[,c("stage","batch","cellCycle")]
sample_metadata = data.frame(
  sample = colnames(hetero),
  condition = cells_meta[colnames(hetero),]$stage, 
  batch = cells_meta[colnames(hetero),]$batch
)
samples_metadata(exp_hetero_as_model.svi) <- sample_metadata
samples_metadata(exp_as_model.svi) <- sample_metadata
samples_metadata(exp_hetero_model.svi) <- sample_metadata
samples_metadata(exp_model.svi) <- sample_metadata

## UMAP and check the accuracy (ARI)
prepare_for_ARI=function(raw_identity, cluster_res) {
  d1= data.table("stage"=as.factor(raw_identity))
  d1[,cell:=names(raw_identity)]
  
  d2=data.table("cluster" = cluster_res$cluster)
  d2[,cell:=names(cluster_res$cluster)]
  
  d= merge(d1,d2, by="cell")
  return(d)
}

get_ari=function(raw_identity, cluster_res) {
  data_for_ARI = prepare_for_ARI(raw_identity, cluster_res)
  ari_value = ARI(data_for_ARI$stage, data_for_ARI$cluster) 
  nmi_value = NMI(data_for_ARI$stage, data_for_ARI$cluster)
  message("ARI: ", ari_value)
  message("NMI: ", nmi_value)
  return(ari_value)
  
}


plot_tsne_cluster=function(tsne_res_dt, tsne_cluster, title="TSNE plot") {
  tsne_res_dt[,c("stage", "treat", "type", "libid", "xx", "batch"):=tstrsplit(sample,"_")]
  tsne_res_dt[,cluster:=as.factor(tsne_cluster$cluster[sample])]  
  
  myplot = ggplot(tsne_res_dt, aes(TSNE1, TSNE2, col=stage, shape=cluster)) +
    geom_point(size=3) +
    theme_bw() +
    ggtitle(title)
  # myplot 
  return(myplot)
}

plot_umap_cluster=function(umap_res_dt, umap_cluster, title="TSNE plot") {
  umap_res_dt[,c("stage", "treat", "type", "libid", "xx", "batch"):=tstrsplit(sample,"_")]
  umap_res_dt[,cluster:=as.factor(umap_cluster$cluster[sample])]  
  
  myplot = ggplot(umap_res_dt, aes(UMAP1, UMAP2, col=stage, shape=cluster)) +
    geom_point(size=3) +
    theme_bw() +
    ggtitle(title)
  return(myplot)
}

#
library('aricode')
raw_identity = as.factor(cells_meta$stage)
names(raw_identity) = rownames(cells_meta)

#### for each model: cluster based on TSNE/UMAP and calculate the ARI  
plot_tsne_allin1 = function(model, model_name) {
    require(aricode)
    raw_identity = as.factor(cells_meta$stage)
    names(raw_identity) = rownames(cells_meta)

    set.seed(123)
    model <- run_tsne(model) ##, min_dist = 0.01, n_neighbors = 10)  
    tsne_dt = data.table(model@dim_red$TSNE)

    model_cluster = cluster_samples(model, k=4, factors=1:model@dimensions$K)
    model_cluster_ari = round(get_ari(raw_identity, model_cluster),3)
    tsne_plot = plot_tsne_cluster(tsne_dt, model_cluster, title = str_c(model_name, ": ARI=", model_cluster_ari)) 
}

plot_umap_allin1 = function(model, model_name) {
    require(aricode)
    raw_identity = as.factor(cells_meta$stage)
    names(raw_identity) = rownames(cells_meta)
    set.seed(123)
    model <- run_umap(model, min_dist = 0.01, n_neighbors = 10)  
    umap_dt = data.table(model@dim_red$UMAP)

    model_cluster = cluster_samples(model, k=4, factors=1:model@dimensions$K)
    model_cluster_ari = round(get_ari(raw_identity, model_cluster),3)
    umap_plot = plot_umap_cluster(umap_dt, model_cluster, title = str_c(model_name, ": ARI=", model_cluster_ari)) 
}


exp_hetero_as_model_tsne = plot_tsne_allin1(exp_hetero_as_model.svi, "EXP+HETERO+AS")
exp_as_model_tsne = plot_tsne_allin1(exp_as_model.svi, "EXP+AS")
exp_hetero_model_tsne = plot_tsne_allin1(exp_hetero_model.svi, "EXP+HETERO")
exp_tsne = plot_tsne_allin1(exp_model.svi, "EXP")

pdf(str_c(outDir, "clustering_plot_TSNE.pdf"), height=12, width=10)
ggpubr::ggarrange(exp_hetero_as_model_tsne, exp_as_model_tsne, 
                exp_hetero_model_tsne, exp_tsne, 
                common.legend = T,
                ncol=2, nrow=2,
                labels = 'AUTO')

dev.off()


# exp_hetero_as_model_umap = plot_umap_allin1(exp_hetero_as_model.svi, "EXP+HETERO+AS")
# exp_as_model_umap = plot_umap_allin1(exp_as_model.svi, "EXP+AS")
# exp_hetero_model_umap = plot_umap_allin1(exp_hetero_model.svi, "EXP+HETERO")
# exp_umap = plot_umap_allin1(exp_model.svi, "EXP")

# ####
# pdf(str_c(outDir, "clustering_plot_UMAP.pdf"), height=12, width=10)
# ggpubr::ggarrange(exp_hetero_as_model_umap, exp_as_model_umap, 
#                 exp_hetero_model_umap, exp_umap, 
#                 common.legend = T,
#                 ncol=2, nrow=2,
#                 labels = 'AUTO')

# dev.off()
