install_github("jokergoo/ComplexHeatmap")
source("http://bioconductor.org/biocLite.R")
biocLite("limma")
library(limma)
library(ComplexHeatmap)
library(RColorBrewer)


PD1_Response <- read.table("~/Documents/PhDProjects/PD1_Resistance_Signature/Data/GSE78220_Response.txt",row.names=1, header = T)
View(PD1_Response)

load("~/Documents/PhDProjects/PD1_Resistance_Signature/Results/PD1_results.rda")

head(test)

design <- model.matrix(~ factor(PD1_Response[,1]))
colnames(design) <- c("ALL", "1vs2")
fit <- lmFit(test, design)
fit <- eBayes(fit)
DEgenes1 <- topTable(fit,coef=2, number=100)
DEgenes50 <- topTable(fit,coef=2, number=50)
View(DEgenes1)
DEgenes1[64,]
row.names(DEgenes1)
write.table(DEgenes1, "~/Documents/PhDProjects/PD1_Resistance_Signature/Results/Top100DEGenes_RvNS.txt", sep='\t',quote=F, col.names = T, row.names= T)


RnaSeq_DEgenes=(test[row.names(DEgenes1),])
RnaSeq_DEgenes_50=(test[row.names(DEgenes50),])

View(RnaSeq_DEgenes_50)
RnaSeq_DEgenes_t=t(RnaSeq_DEgenes)
RnaSeq_DEgenes_t_50=t(RnaSeq_DEgenes_50)
RnaSeq_DEgenes_t_scaled_50=scale(RnaSeq_DEgenes_t_50)
View(RnaSeq_DEgenes_t_50)
RnaSeq_DEgenes_t_

annos <- HeatmapAnnotation(PD1_Response[,"Response",drop=F],
                           name="EGFR",
                           col = list(Response = c("R" =  "deeppink","NR" =  "darkolivegreen1")),
                           which = "row",
                           width = unit(0.5, "cm"),
                           show_legend = T)

hm <- Heatmap(t(RnaSeq_DEgenes), cluster_rows = F,
                  cluster_columns = F, show_row_names = F, show_column_names = F,
                  row_title_gp = gpar(fontsize =10), combined_name_fun = NULL,
                  name="RNA-Seq",col=my_palette,
                  column_title = "RNA-seq Values of Top 100 Diff. Expressed Genes (R vs NS)",
                  show_heatmap_legend = T, column_title_gp = gpar(fontsize = 9, fontface = "bold"))

PD1_Response
summary(PD1_Response)
ha = HeatmapAnnotation(df = PD1_Response)

df = data.frame(Response = c(rep("Responders", 13), rep("NonResponders", 15)))
df
ha = HeatmapAnnotation(df = df, col = list(Response = c("Responders" =  "deeppink", "NonResponders" = "darkolivegreen1")))

hm_scaled <- Heatmap(t(RnaSeq_DEgenes_t_scaled), cluster_rows = T,
              cluster_columns = F, show_row_names = F, show_column_names = F,
              row_title_gp = gpar(fontsize =10), combined_name_fun = NULL,
              name="RNA-Seq Values",col=my_palette,top_annotation = ha,
              column_title = "Top 100 Diff. Expressed Genes (R vs NS)",
              show_heatmap_legend = T, column_title_gp = gpar(fontsize = 12, fontface = "bold"))

hm_scaled <- Heatmap(t(RnaSeq_DEgenes_t_50), cluster_rows = T,
                     cluster_columns = F, show_row_names = T, show_column_names = F,
                     row_names_gp = gpar(fontsize = 8),
                     row_title_gp = gpar(fontsize =10), combined_name_fun = NULL,
                     name="RNA-Seq Values",col=my_palette,top_annotation = ha,
                     column_title = "Top 50 Diff. Expressed Genes  ",
                     show_heatmap_legend = T, column_title_gp = gpar(fontsize = 12, fontface = "bold"))
pdf("~/Documents/PhDProjects/PD1_Resistance_Signature/Results/Top50DEGenes_PD1.pdf")
draw(hm_scaled, row_title = "Genes", column_title_side = "bottom",annotation_legend_side = "bottom" )
dev.off()
