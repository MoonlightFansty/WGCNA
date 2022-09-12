#加载seurat数据和包
# single-cell analysis package
setwd('e:/writing/benke/hdWGCNA/')
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# gene enrichment packages
#install.packages('enrichR')
library(enrichR)

#BiocManager::install('GeneOverlap',update = F,ask = F)
library(GeneOverlap)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# load the Zhou et al snRNA-seq dataset
seurat_obj <- readRDS('hdWGCNA_object.rds')

#GO富集分析
# enrichr databases to test
dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021')

# 富集分析，会逐个模块分析
seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test
)

# retrieve the output table
enrich_df <- GetEnrichrTable(seurat_obj)

# make GO term plots作图，在文件夹下生成！
EnrichrBarPlot(
  seurat_obj,
  outdir = "enrichr_plots", # name of output directory
  n_terms = 10, # number of enriched terms to show (sometimes more show if there are ties!!!)
  plot_size = c(5,7), # width, height of the output .pdfs
  logscale=TRUE # do you want to show the enrichment as a log scale?
)

#气泡图
# GO_Biological_Process_2021
EnrichrDotPlot(
  seurat_obj,
  mods = c("turquoise","black"), # use all modules (this is the default behavior)
  database = "GO_Biological_Process_2021", # this has to be one of the lists we used above!!!
  n_terms=6
  # number of terms for each module
)


#气泡图
# GO_Cellular_Component_2021
EnrichrDotPlot(
  seurat_obj,
  mods = c("turquoise","black",'blue'), # use all modules (this is the default behavior)
  database = "GO_Cellular_Component_2021", # this has to be one of the lists we used above!!!
  n_terms=6
  # number of terms for each module
)

#气泡图
# GO_Biological_Process_2021
EnrichrDotPlot(
  seurat_obj,
  mods = c("turquoise","black",'red'), # use all modules (this is the default behavior)
  database = "GO_Molecular_Function_2021", # this has to be one of the lists we used above!!!
  n_terms=6
  # number of terms for each module
)

#差异基因重叠分析
## 这个分析帮助我们看到，哪些模块可能是相似的
# compute cell-type marker genes with Seurat:
# 常规方法计算差异基因/特征基因
Idents(seurat_obj) <- seurat_obj$seurat_clusters
markers <- Seurat::FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  logfc.threshold=1
)

# compute marker gene overlaps
overlap_df <- OverlapModulesDEGs(
  seurat_obj,
  deg_df = markers,
  fc_cutoff = 1 # log fold change cutoff for overlap analysis
)

#条形图
# overlap barplot, produces a plot for each cell type
plot_list <- OverlapBarPlot(overlap_df)

# stitch plots with patchwork
wrap_plots(plot_list, ncol=4)

#气泡图
# plot odds ratio of the overlap as a dot plot
OverlapDotPlot(
  overlap_df,
  plot_var = 'odds_ratio') +
  ggtitle('Overlap of modules & cell-type markers')

#----------------------------
#网络可视化
# network analysis & visualization package:
# network analysis & visualization package:
library(igraph)

#可视化每个模块的网络图
ModuleNetworkPlot(seurat_obj)

#组合网络图，在文件夹下生成
# hubgene network
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 3, n_other=5,
  edge_prop = 0.75,
  mods = "all"
)

#UMAP可视化
## 利用hub基因，重新UMAP，如此可以获得分群明显的图
g <- HubGeneNetworkPlot(seurat_obj,  return_graph=TRUE)

seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)


# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

# plot with ggplot
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()

ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=2 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
  
)

#监督UMAP
g <- ModuleUMAPPlot(seurat_obj,  return_graph=TRUE)
# run supervised UMAP:
seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10,
  n_neighbors=15,
  min_dist=0.1,
  supervised=TRUE,
  target_weight=0.5
)

# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

# plot with ggplot
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()