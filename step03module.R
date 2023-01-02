power = sft$powerEstimate

net = blockwiseModules(wgcna_matrix, power = 2,
                       corType = "pearson",
                       networkType="unsigned",
                       TOMType = "unsigned", 
                       minModuleSize = 30,
                       mergeCutHeight = 0.25,
                       numericLabels = , 
                       verbose = 3)

# 最终得到的网络模块(合并之后)
unique(net$colors)
table(net$colors)
# grey module是 unassigned gene

# 合并之前的网络模块
unique(net$unmergedColors)
table(net$unmergedColors)

# 样本对于模块的特征值
dim(net$MEs)
net$MEs[1:4,1:4]
net$MEsOK


moduleColors = net$colors

moduleColors <- ifelse(moduleColors == 'blue', '#0074b3',
                ifelse(moduleColors == 'cyan', '#abd9e9',
                ifelse(moduleColors == 'green', '#459943',
                ifelse(moduleColors == 'magenta', '#de77ae', 
                ifelse(moduleColors == 'midnightblue', '#313695', 
                ifelse(moduleColors == 'orange', '#f47720',
                ifelse(moduleColors == 'purple', '#762a83', 
                ifelse(moduleColors == 'red', '#982b2b',
                ifelse(moduleColors == 'yellow', '#e8c559',
                ifelse(moduleColors == 'black', 'black',
                ifelse(moduleColors == 'brown', 'brown',
                ifelse(moduleColors == 'darkgreen', 'darkgreen',
                ifelse(moduleColors == 'darkgrey', 'darkgrey', 
                ifelse(moduleColors == 'darkred', 'darkred', 
                ifelse(moduleColors == 'darkturquoise', 'darkturquoise',
                ifelse(moduleColors == 'grey60', 'grey60', 
                ifelse(moduleColors == 'lightcyan', 'lightcyan',
                ifelse(moduleColors == 'lightgreen', 'lightgreen',
                ifelse(moduleColors == 'lightyellow', 'lightyellow',
                ifelse(moduleColors == 'pink', 'pink',
                ifelse(moduleColors == 'royalblue', 'royalblue',
                ifelse(moduleColors == 'salmon', 'salmon', 
                ifelse(moduleColors == 'tan', 'tan', 
                ifelse(moduleColors == 'turquoise', 'turquoise',
                ifelse(moduleColors == 'greenyellow', 'greenyellow', 'grey')))))))))))))))))))))))))

# Plot the dendrogram and the module colors underneath
# 如果对结果不满意,还可以recutBlockwiseTrees
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

color <- as.data.frame(net$colors)
color <- cbind(color, row.names(color))
colnames(color) <- c('color', 'symbol')
gene_set <- color[gene_set, ]

brown <- subset(color, subset = color == 'brown')
turquoise <- subset(color, subset = color == 'turquoise')
lightcyan <- subset(color, subset = color == 'lightcyan')

save(color, brown, turquoise, lightcyan, gene_set, file = 'Gene_set.Rdata')

save(net, file = 'Net.Rdata')
