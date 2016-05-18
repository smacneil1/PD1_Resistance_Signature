library(ComplexHeatmap)
library(RColorBrewer)
library(ggplots2)

baseDir<- "~/Documents/PhDProjects/PD1_Resistance_Signature/Results/ASSIGN/"
gatherFile




data<-gatherFile(baseDir)
View(baseDir)
setwd(baseDir)
filenames<-system("ls */*/pathway_activity_testset*", intern=TRUE)
filenames
colnames(data)<-gsub(pattern = "/pathway_activity_testset.csv",replacement = "",x = colnames(data))

data=NULL
for(i in 1:length(filenames)){
  ###reading in the filess one at a time
  f<-read.csv(filenames[i], header=1,row.names=1)
  #colnames(f)<-paste(filenames[i],colnames(f),sep='/')
  sig_name=gsub('\\_.*$','',filenames[i])
  colnames(f)<-paste(toupper(sig_name))
  if(i==1){
    data<-f
  }
  else{
    data<-cbind(data,f)
  }
}
row.names(data)

dim(data)
data=data[,-7]
row.names(GSE78220_Response)
GSE78220_Response
data_resp=as.matrix(merge_drop(data,GSE78220_Response))


data_resp

responders= subset(data_resp, subset = "R")
non-responders=
pdf("~/Documents/PhDProjects/PD1_Resistance_Signature/Results/HeatMaps_ASSIGN_PD1.pdf")

row_R <- HeatmapAnnotation(GSE78220_Response[,"Response",drop=F],
                                name="EGFR",
                                col = list(Response = c("R" =  "deeppink","NR" =  "darkolivegreen1")),
                                which = "row",
                                width = unit(0.5, "cm"),
                                show_legend = T)

PD1_hm <- Heatmap(scale(data), cluster_rows = T,
              cluster_columns = T, show_row_names = F, show_column_names = T,
              row_title_gp = gpar(fontsize =10), combined_name_fun = NULL,
              name="Scaled\nPathway\nActivity",col=my_palette,
              column_title = "Pathway Activation Across Patients Treated PD-1 Inhibitor (Scaled)",
              show_heatmap_legend = T, column_title_gp = gpar(fontsize = 9, fontface = "bold"))
draw(PD1_hm+ row_R)

PD1_hm <- Heatmap(data, cluster_rows = T,
                  cluster_columns = T, show_row_names = F, show_column_names = T,
                  row_title_gp = gpar(fontsize =10), combined_name_fun = NULL,
                  name="Non-Scaled\nPathway\nActivity",col=my_palette,
                  column_title = "Pathway Activation Across Patients Treated PD-1 Inhibitor (Not Scaled)",
                  show_heatmap_legend = T, column_title_gp = gpar(fontsize = 9, fontface = "bold"))



draw(PD1_hm+ row_R)
dev.off()
################################################

#try plotting the genes from the gene list

geneListBAD=as.matrix(c(bad_gene_list_DOWN, bad_gene_list_UP))
geneListBAD
View(geneListBAD)
geneListBAD_m=merge(geneListBAD,test, by.x=1, by.y=0)
rownames(geneListBAD_m)=geneListBAD_m[,1]
head(geneListBAD_m)
geneListBAD_m=geneListBAD_m[,-1]

geneListBAD_m_combat=merge(geneListBAD,dat, by.x=1, by.y=0)
dim(geneListBAD_m_combat)

geneListAKT=as.matrix(c(akt_gene_list_DOWN, akt_gene_list_UP))
geneListAKT_m_combat=merge(geneListAKT,dat, by.x=1, by.y=0)
dim(geneListAKT_m_combat)


geneListHER2=as.matrix(c(her2_gene_list_DOWN, her2_gene_list_UP))
dim(geneListHER2) #10
geneListHER2_m_combat=merge(geneListHER2,dat, by.x=1, by.y=0)
View(geneListHER2_m_combat)

geneListIGF1R=as.matrix(c(igf1r_gene_list_DOWN, igf1r_gene_list_UP))
dim(geneListIGF1R) #100
geneListIGF1R_m_combat=merge(geneListIGF1R,dat, by.x=1, by.y=0)
dim(geneListIGF1R_m_combat) #97

geneListKRASGV=as.matrix(c(krasgv_gene_list_DOWN, krasgv_gene_list_UP))
dim(geneListKRASGV) #175
geneListKRASGV_m_combat=merge(geneListKRASGV,dat, by.x=1, by.y=0)
View(geneListKRASGV_m_combat) #168

geneListRAF=as.matrix(c(raf_gene_list_DOWN, raf_gene_list_UP))
dim(geneListRAF)
geneListRAF_m_combat=merge(geneListRAF,dat, by.x=1, by.y=0)
dim(geneListRAF_m_combat) #339

geneListEGFR=as.matrix(c(egfr_gene_list_DOWN, egfr_gene_list_UP))
dim(geneListEGFR) #50
geneListEGFR_m_combat=merge(geneListEGFR,dat, by.x=1, by.y=0)
dim(geneListEGFR_m_combat) #50

geneListKRADQH=as.matrix(c(krasqh_gene_list_DOWN, krasqh_gene_list_UP))
dim(geneListKRADQH) #300
geneListKRADQH_m_combat=merge(geneListKRADQH,dat, by.x=1, by.y=0)
dim(geneListKRADQH_m_combat) #288

test_t=
list(geneListBAD)
as.matrix(test[,geneListBAD])
as.matrix(t(test[geneListBAD,]))
dim(test)



class(test)

row_RNAseq<- HeatmapAnnotation(GSE78220_Response[,"Response",drop=F],
                               name="Response",
                               col = list(Response = c("R" =  "deeppink","NR" =  "darkolivegreen1")),
                               which = "row",
                               width = unit(0.5, "cm"),
                               show_legend = T)



geneListBAD_m_t=t(geneListBAD_m)
geneListBAD_m_t
geneListBAD_m_t_no28=geneListBAD_m_t[-12,]

geneListBAD_m_t
geneListBAD_m_scaled_no28=scale(geneListBAD_m_t_no28)
geneListBAD_m_scaled=scale(geneListBAD_m_t)

View(geneListBAD_m_scaled)


RNAseq_PD_hm <- Heatmap(as.matrix(geneListBAD_m_scaled_no28), cluster_rows = T,
                  cluster_columns = T, show_row_names = F, show_column_names = F,
                  row_title_gp = gpar(fontsize =8), combined_name_fun = NULL,
                  name="RNA-seq values",col=my_palette,
                  column_title = "RNA-seq Values for BAD signature w/o 28",
                  show_heatmap_legend = T, column_title_gp = gpar(fontsize = 9, fontface = "bold"))

draw(RNAseq_PD_hm+row_RNAseq)

RNAseq_PD_hm <- Heatmap(as.matrix(geneListBAD_m_scaled), cluster_rows = T,
                        cluster_columns = T, show_row_names = F, show_column_names = F,
                        row_title_gp = gpar(fontsize =8), combined_name_fun = NULL,
                        name="RNA-seq values",col=my_palette,
                        column_title = "RNA-seq Values for BAD signature with 28",
                        show_heatmap_legend = T, column_title_gp = gpar(fontsize = 9, fontface = "bold"))


draw(RNAseq_PD_hm+row_RNAseq)
