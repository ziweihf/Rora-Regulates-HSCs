#TFs score----

library(ggplot2)
library(patchwork)
library(grid)
library(ggtext)
library(ggsignif)
library(stringr)
bg_colors <- c("white","#E8DFF0","#BFA0CC")
fill_color = "#A184BC"
gradient_grob <- rasterGrob(colorRampPalette(bg_colors)(256), width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)
ggplot(female_count, #Ucell
       aes(x = group, y = IduceTF, fill = group)) +
  annotation_custom(gradient_grob,xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  geom_violin(trim = FALSE) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
               geom="crossbar", width=0.3, size = 0.3) +
  geom_signif(comparisons = list(c("Old","young")), 
              map_signif_level = function(p) {paste("italic(P) == ", sprintf("%.2g", p))                },
              y_position = 0.23,  
              textsize = 4, tip_length = 0,
              parse = TRUE) +
  scale_fill_manual(values = c( fill_color,"#9C9C9C")) +
  labs(title="IduceTF<br>GSE147729 test set", x = NULL, y = "Gene set score") +
  theme_classic(base_size = 12) +
  theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10),
        plot.title = element_textbox_simple(size = 15, color = "white", halign = 0.5,
                                            fill = fill_color, width = 1.2, 
                                            padding = margin(3, 0, 3, 0),
                                            margin = margin(0, 0, 10, 0)),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black"),
        axis.text.y = element_text(size = 11,color = "black"),
        legend.position = "none")

#----------resist---
ggplot(female_count, 
       aes(x = group, y = ResistTF, fill = group)) +
  annotation_custom(gradient_grob,xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  geom_violin(trim = FALSE) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
               geom="crossbar", width=0.3, size = 0.3) +
  geom_signif(comparisons = list(c("Old","young")), 
              map_signif_level = function(p) {paste("italic(P) == ", sprintf("%.2g", p))                },
              y_position = 0.19,  
              textsize = 4, tip_length = 0,
              parse = TRUE) +
  scale_fill_manual(values = c( fill_color,"#9C9C9C")) +
  labs(title="ResistTF<br>GSE147729 test set", x = NULL, y = "Gene set score") +
  theme_classic(base_size = 12) +
  theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10),
        plot.title = element_textbox_simple(size = 15, color = "white", halign = 0.5,
                                            fill = fill_color, width = 1.2, 
                                            padding = margin(3, 0, 3, 0),
                                            margin = margin(0, 0, 10, 0)),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black"),
        axis.text.y = element_text(size = 11,color = "black"),
        legend.position = "none")



#GSEA plot----
plotnes <- read.csv("../data/GSEA/YvsO_GSEA.csv")
plotnes %>% filter(abs(NES)>1,pvalue<0.05) %>% 
  ggplot(aes(reorder(ID,NES), NES)) +
  geom_col(aes(fill= NES)) +
  scale_fill_gradient2(low = "#333333",high = "#0066CC",mid = "#CCCCCC")+
  coord_flip() +
  labs(x="Pathways", y="Normalized Enrichment Score",title="Human vs Young")+theme_bw()+
  theme_classic()+
  theme(plot.title = element_text(size = 15,
                                  colour = "black",
                                  hjust = 0.5),
        axis.title.y = element_text(size = 15, 
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        axis.title.x = element_text(size = 15, 
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5),
        axis.text.y = element_text(size = 13,  
                                   color = "black",
                                   face = "bold"),
        axis.text.x = element_text(size = 13,  
                                   color = "black",
                                   face = "bold", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0),
        legend.title = element_text(color="black", 
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", 
                                   size = 10, 
                                   face = "bold"))
