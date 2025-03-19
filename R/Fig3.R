#after CCA-----------------------
library(Seurat)
library(dplyr)
load("../data/afterCCA_seurat_integrated.Rdata")
seurat_integrated <- RunPCA(object = seurat_integrated)
# Determine percent of variation associated with each PC
pct <- seurat_integrated[["pca"]]@stdev / sum(seurat_integrated[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2

seuratintegated <- seurat_integrated %>%
  FindNeighbors(dims = 1:45) %>%
  FindClusters(resolution = 1) %>%
  RunTSNE(dims = 1:45) %>%
  RunUMAP(dims = 1:45)


########################################
###################plot#################
#########################################
#marker
VlnPlot(seuratintegated, features = c("Cd48","Slamf1","Flt3","Cd34"), stack = TRUE, sort = F,cols = my36colors,group.by = "sample_celltype") +
  theme(legend.position = "none") + ggtitle("CellType_Features")

#umap--------------------------
library(tidydr)
library(tidyverse)
library(Seurat)
library(ggplot2)
library(viridis)
library(RColorBrewer)
DimPlot(seuratintegated, label = F, group.by = "sample_celltype",repel = T,reduction = "umap",
        cols = c("#D9534F","#CBE86B","#0099CC"))+
  theme_dr(xlength = 0.2,
           ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(aspect.ratio = 1,
        panel.grid = element_blank())

DimPlot(seuratintegated,group.by = "group",cols =  c("#E15759","#9EB9F3"),reduction = "umap",pt.size = .1)+
  theme_dr(xlength = 0.2,
           ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(aspect.ratio = 1,
        panel.grid = element_blank())

#violin score----------------------
#from UCell
library(introdataviz)
forpathways <- grep("HALLMARK",rownames(female_count),value = T)
female_count$celltype <- factor(female_count$celltype ,levels = c("LTHSC","STHSC",
                                                                  "MPP"))
for (i in forpathways) {
  pathway=i
  cat(paste("Processing",pathway,"\n",sep = " "))
  female_countfig <- dplyr::select(female_count,c("celltype",pathway,"group"))
  colnames(female_countfig)[2] <- "selectecols" 
  ggfig <- ggplot(female_countfig,aes(x = celltype,y = selectecols,fill = group)) +
    # split violin
    geom_split_violin(alpha = .5, trim = F,width = 0.7,color=NA,scale = "area") + 
    scale_fill_manual(values = c('#d22523','#3f7bbe'),name = "group")+
    # mean point
    stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
    theme_bw(base_size = 16) +
    labs(y=pathway)+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    theme(axis.text.x = element_text(angle = 90,color = 'black',hjust = 1),
          axis.text.y = element_text(color = 'black'),
          legend.position = 'top') +
    # add p
    stat_compare_means(aes(group=group),method = "wilcox",
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "NS")),label = "p.signif",
                       size = 5)
  ggsave(ggfig,filename = paste("./results/paperFig3/phenotype/",pathway,".pdf",sep = ""), height = 4,width = 5)
  
}
