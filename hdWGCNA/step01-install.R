# install BiocManager
install.packages("BiocManager")

# install Bioconductor core packages
BiocManager::install("Seurat",update=F,ask=F)
BiocManager::install("WGCNA",update=F,ask=F)

# install additional packages:
install.packages(c("igraph", "devtools"))


# 正式安装
library(devtools)
devtools::install_github('smorabit/hdWGCNA', ref='dev')
# 当然运行不畅记得用本地下载install_local

