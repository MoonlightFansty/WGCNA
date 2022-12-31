# 计算指定表型的相关分析
pheno = 'stage'
pheno_type = phenotype[ , pheno, drop=F]
colnames(pheno_type) = pheno
# (1) Gene significance,GS: 即比较样本某个基因与对应表型的相关性
GS = as.data.frame(cor(wgcna_matrix, pheno_type, use = "p"))
colnames(GS) = "GS"

GS.p = as.data.frame(corPvalueStudent(as.matrix(GS), nrow(wgcna_matrix)))


# (2) Module Membership: 模块内基因表达与模块特征值的相关性
modNames = substring(names(MEs), 3)
# 计算3600个基因与18个模块的相关性
MM = as.data.frame(cor(wgcna_matrix, MEs, use = "p"))
colnames(MM) = paste("MM", modNames, sep="")

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(MM), nrow(wgcna_matrix)))
colnames(MMPvalue) = paste("p.MM", modNames, sep="");
MMPvalue[1:4,1:4]

# # (3) 可视化blue模块基因特征
# module = "turquoise"
# moduleGenes = names(net$colors)[net$colors==module]
# sizeGrWindow(7, 7)
# par(mfrow = c(1,1))
# verboseScatterplot(abs(MM[moduleGenes, paste("MM", module, sep="")]),
#                    abs(GS[moduleGenes, 1]),
#                    xlab = paste("Module Membership in", module, "module"),
#                    ylab = paste("Gene significance for", pheno),
#                    main = paste("Module membership vs. gene significance\n"),
#                    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


module = 'lightcyan'
column = match(module, modNames)
moduleGenes = net$colors==module
table(moduleGenes)
moduleGenes
module <- as.data.frame(dimnames(data.frame(wgcna_matrix))[[2]][moduleGenes]) 
names(module)='genename'


MM_1 <- abs(MM[moduleGenes,column])
GS_1 <- abs(GS[moduleGenes, 1])
b <- as.data.frame(cbind(MM_1,GS_1))
rownames(b) = module$genename

hub <- abs(b$MM_1)>0.8 & abs(b$GS_1)>0.12
table(hub)

b$group <- hub

library(ggplot2)
library(ggpubr)

ggplot(data=b) +
  geom_point(size=1.5, aes(x=MM_1, y=GS_1, color=group)) +
  scale_colour_manual(values=c("grey60", "#DE6757")) + 
  theme_bw() +  
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +  
  labs(x="Module Membership in lightcyan module", y="Gene significance for stage",title = "Module membership vs. gene significance")+ 
  theme(axis.title.x =element_text(size=12), axis.title.y=element_text(size=12),axis.text = element_text(size = 12),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"),plot.title = element_text(hjust = 0.5,size = 12,face = "bold"),plot.margin = unit(rep(2,4),'lines')) +
  theme(legend.position = 'none') +
  geom_hline(aes(yintercept=0.12),colour="#5B9BD5",lwd=1,linetype=5) +
  geom_vline(aes(xintercept=0.8),colour="#5B9BD5",lwd=1,linetype=5) +
  stat_cor(method = 'pearson', aes(x=MM_1, y=GS_1), color="#DE6757")
