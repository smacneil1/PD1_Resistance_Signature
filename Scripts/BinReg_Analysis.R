######
# Runs the t.tests on BinReg predicts for the PD1 data
# Shelley Macneil
# May-16-2016
#
######

library(ComplexHeatmap)
source('~/Documents/PhDProjects/GFRN_signatures-master/Key_ASSIGN_functions_balancedsig.R')
setwd("~/Documents/PhDProjects/PD1_Resistance_Signature/Results/BinReg/Probabilities/")

files=system("ls", intern=TRUE)
files
response <- read.table("~/Documents/PhDProjects/PD1_Resistance_Signature/Data/GSE78220_Response.txt",row.names=1, header = T)


# function to read in a list of filenames and store as a variable (add to list)
for(i in 1:length(files)){
  ###reading in the filess one at a time
  f<-read.table(files[i], header=1,sep='\t')
  var_name=gsub('\\..*$','',files[i])
  print(var_name)
  assign(var_name, f) 
}  
length(files)


probabilities_RAF_50=probabilities_RAF_50[37:64, c(2,5) ]
probabilities_RAF_100=probabilities_RAF_100[37:64,c(2,5) ]
probabilities_RAF_175=probabilities_RAF_175[37:64,c(2,5) ]
probabilities_RAF_300=probabilities_RAF_300[37:64,c(2,5) ]
probabilities_RAF_350=probabilities_RAF_350[37:64,c(2,5) ]
probabilities_RAF_500=probabilities_RAF_500[37:64, c(2,5)]

probabilities_AKT_4=probabilities_AKT_4[37:64,c(2,5)]
probabilities_AKT_10=probabilities_AKT_10[37:64,c(2,5)]
probabilities_AKT_100=probabilities_AKT_100[37:64,c(2,5)]
probabilities_AKT_20=probabilities_AKT_20[37:64,c(2,5)]
probabilities_AKT_50=probabilities_AKT_50[37:64,c(2,5)]

probabilities_IGF1R_100=probabilities_IGF1R_100[37:64,c(2,5)]
probabilities_IGF1R_300=probabilities_IGF1R_300[37:64,c(2,5)]
probabilities_IGF1R_50=probabilities_IGF1R_50[37:64,c(2,5)]

probabilities_EGFR_25=probabilities_EGFR_25[25:52,c(2,5)]  
probabilities_EGFR_50=probabilities_EGFR_50[25:52,c(2,5)]
probabilities_EGFR_100=probabilities_EGFR_100[25:52,c(2,5)]  

probabilities_KRAS_100=probabilities_KRAS_100[37:64,c(2,5) ]
probabilities_KRAS_250=probabilities_KRAS_250[37:64,c(2,5) ]
probabilities_KRAS_50=probabilities_KRAS_50[37:64,c(2,5) ]

probabilities_HER2_10=probabilities_HER2_10[35:62,c(2,5) ]
probabilities_HER2_25=probabilities_HER2_25[35:62,c(2,5) ]

probabilities_BAD_100=probabilities_BAD_100[37:64,c(2,5) ]
probabilities_BAD_250=probabilities_BAD_250[37:64,c(2,5) ]
probabilities_BAD_350=probabilities_BAD_350[37:64,c(2,5) ]
probabilities_BAD_50=probabilities_BAD_50[37:64,c(2,5) ]



# make the boxplot function

plotBinReg= function (probs, main){
  rownames(probs)=probs[,1]
  probs=merge_drop(probs, response, by=0)
  probs=probs[order(probs$Response),] 
  par(cex.lab=1.3)
  par(cex.axis=1.3)
  margins=c(7,8)
  boxplot(probs$Probability~probs$Response, main = main, ylab="Pathway Activity", col=c("darkolivegreen1","deeppink" ), names=c("Non-Responders", "Responders"))
  ttest=t.test(probs$Probability~probs$Response)
  title(sub=paste("P-value = ", format(ttest$p.value, digits=3)), adj=1, line=3, font=1)
  return(probs)
}

pdf("~/Documents/PhDProjects/PD1_Resistance_Signature/Results/PD1_Boxplots_Andrea.pdf")

plotBinReg(probabilities_AKT_20, "AKT Pathway")
plotBinReg(probabilities_RAF_350, "RAF Pathway")
plotBinReg(probabilities_IGF1R_100, "IGF1R Pathway")
plotBinReg(probabilities_EGFR_50, "EGFR Pathway")
plotBinReg(probabilities_KRAS_100, "KRAS Pathway")
dev.off()

# plot and get the p.val from the ttest 

pdf("~/Documents/PhDProjects/PD1_Resistance_Signature/Results/Boxplots.pdf")


plotBinReg(probabilities_RAF_50) # 0.87
plotBinReg(probabilities_RAF_100) # 0.32
plotBinReg(probabilities_RAF_175) ##0.27

plotBinReg(probabilities_RAF_300) #  0.14
plotBinReg(probabilities_RAF_350) #0.13   #  go with this one 
plotBinReg(probabilities_RAF_500) #0.17

plotBinReg(probabilities_AKT_4)

plotBinReg(probabilities_AKT_10)  # 0.16
plotBinReg(probabilities_AKT_20) # 0.18 # go with this one 
plotBinReg(probabilities_AKT_50) # 0.7
plotBinReg(probabilities_AKT_100) # 0.99

plotBinReg(probabilities_IGF1R_50) # 0.38
plotBinReg(probabilities_IGF1R_100) # 0.009   #go with this one
plotBinReg(probabilities_IGF1R_300) #  0.01


plotBinReg(probabilities_EGFR_25 ) # #0.044
plotBinReg(probabilities_EGFR_50  ) #0.044  # go with this one
plotBinReg(probabilities_EGFR_100 ) # #0.06

plotBinReg(probabilities_KRAS_100) # 0.15   # go with this one
plotBinReg(probabilities_KRAS_250) # 0.34
plotBinReg(probabilities_KRAS_50) # 0.17 
plotBinReg(probabilities_HER2_10) #

plotBinReg(probabilities_HER2_10) #0.47
plotBinReg(probabilities_HER2_25) # 0.89

plotBinReg(probabilities_BAD_100) # #0.89   #go with this one
plotBinReg(probabilities_BAD_250) #
plotBinReg(probabilities_BAD_350) #
plotBinReg(probabilities_BAD_50) #  0.64

#clean the data for heatmaps 
rownames(probabilities_EGFR_50) = probabilities_EGFR_50[,1]
probabilities_EGFR_50=within(probabilities_EGFR_50, rm(Sample))
View(probabilities_EGFR_50)

rownames(probabilities_AKT_20) = probabilities_AKT_20[,1]
probabilities_AKT_20=within(probabilities_AKT_20, rm(Sample))
View(probabilities_AKT_20)

rownames(probabilities_IGF1R_100) = probabilities_IGF1R_100[,1]
probabilities_IGF1R_100=within(probabilities_IGF1R_100, rm(Sample))
View(probabilities_IGF1R_100)

rownames(probabilities_KRAS_100) = probabilities_KRAS_100[,1]
probabilities_KRAS_100=within(probabilities_KRAS_100, rm(Sample))
View(probabilities_KRAS_100)

rownames(probabilities_HER2_10) = probabilities_HER2_10[,1]
probabilities_HER2_10=within(probabilities_HER2_10, rm(Sample))
View(probabilities_HER2_10)

rownames(probabilities_BAD_100) =probabilities_BAD_100[,1]
probabilities_BAD_100=within(probabilities_BAD_100, rm(Sample))
View(probabilities_BAD_100)

rownames(probabilities_RAF_350) =probabilities_RAF_350[,1]
probabilities_RAF_350=within(probabilities_RAF_350, rm(Sample))
View(probabilities_RAF_350)



FilesToWrite=c("probabilities_EGFR_50", "probabilities_AKT_20", "probabilities_IGF1R_100", "probabilities_KRAS_100","probabilities_RAF_350", "probabilities_HER2_10", "probabilities_BAD_100")
FilesToWrite=probabilities_EGFR_50

for (i in 1:length(FilesToWrite)) {
variable= get(FilesToWrite[i])
variable
ResponseMerge=cbind(variable,response)
ResponseMerge
write.table(ResponseMerge, file = paste("Response_",FilesToWrite[i], ".txt",sep=""), sep='\t',quote=F, col.names = NA, row.names= T)
}

length(FilesToWrite)

all= cbind(probabilities_EGFR_50,probabilities_AKT_20,probabilities_IGF1R_100,probabilities_KRAS_100,probabilities_HER2_10, probabilities_BAD_100, probabilities_RAF_100   )
colnames(all) = c("EGFR", "AKT", "IGF1R", "KRAS", "HER2", "BAD")

View(all)
# function for ploting the responds vs. nonresponders





                                    ## make the heatmaps 
response
all_m=cbind(all, response)
all_rm=all_rm[,-5]
all_rm=all_rm[,-5]
all_rm= cbind(all_rm, response)
all_rm_or = all_rm[order(all_rm$Response, all_rm$EGFR), ]
all_rm_or_d=all_rm_or[,-5]

GeneList <- read.table("~/Documents/PhDProjects/PD1_Resistance_Signature/Data/PatientSortingPD1.txt", quote="\"", comment.char="")
class(GeneList)

row_EGFR <- HeatmapAnnotation(all_m[, "Response", drop =F], name="Response",
                              col = list(Response = c("R" =  "deeppink", "NR" = "darkolivegreen1")),
                              which = "row",width = unit(0.3, "cm"),show_legend = T)

heat_s <- Heatmap(scale(all_rm), cluster_rows = T,
                  cluster_columns = T, show_row_names = T, show_column_names = T,
                  row_title_gp = gpar(fontsize =10), combined_name_fun = NULL,
                  name="Pathway Predictions", column_title = "ALL",col=my_palette,
                  show_heatmap_legend = T, column_title_gp = gpar(fontsize = 12, fontface = "bold"),  width = unit(70, "mm"))


draw(heat_s + row_EGFR)


R= c("Pt28.baseline", "Pt19.baseline", "Pt35.baseline","Pt38.baseline","Pt27A.baseline", "Pt6.baseline", "Pt13.baseline","Pt9.baseline", "Pt8.baseline", "Pt2.baseline", "Pt15.baseline", "Pt27B.baseline", "Pt37.baseline", "Pt5.baseline", "Pt4.baseline"   )
R
NR= c("Pt10.baseline", "Pt12.baseline", "Pt16.OnTx", "Pt31.baseline", "Pt23.baseline", "Pt25.baseline", "Pt14.baseline", "Pt7.baseline", "Pt32.baseline", "Pt20.baseline", "Pt1.baseline", "Pt29.baseline")
NR
new_order=all_rm[ c(R,NR), ]
new_order


heat <- Heatmap(all_rm, cluster_rows = F,
                cluster_columns = T, show_row_names = F, show_column_names = T,
                row_title_gp = gpar(fontsize =10), combined_name_fun = NULL,
                name="Pathway Predictions", column_title = "ALL",col=my_palette,
                show_heatmap_legend = T, column_title_gp = gpar(fontsize = 12, fontface = "bold"),  width = unit(70, "mm"))

heat_s <- Heatmap(scale(all_rm), cluster_rows = F,
                cluster_columns = T, show_row_names = T, show_column_names = T,
                row_title_gp = gpar(fontsize =10), combined_name_fun = NULL,
                name="Pathway Predictions", column_title = "ALL",col=my_palette,
                show_heatmap_legend = T, column_title_gp = gpar(fontsize = 12, fontface = "bold"),  width = unit(70, "mm"))



##### Ones for Andrea 
pdf("~/Documents/PhDProjects/PD1_Resistance_Signature/Results/PD1_HeatMap_PathwayActivity_NotScaled.pdf")




heat_s_sort <- Heatmap(new_order, cluster_rows = F,
                       cluster_columns = T, show_row_names = F, show_column_names = T,
                       row_title_gp = gpar(fontsize =10), combined_name_fun = NULL,
                       name="Pathway Predictions", column_title = "Pathway Activation",,col=my_palette, row_names_gp = gpar(fontsize = 8),
                       show_heatmap_legend = T, column_title_gp = gpar(fontsize = 12, fontface = "bold"),  width = unit(90, "mm"))

draw(heat_s_sort + row_EGFR)
dev.off()

pdf("~/Documents/PhDProjects/PD1_Resistance_Signature/Results/PD1_HeatMap_PathwayActivity_Scaled.pdf")

##### Ones for Andrea 

pdf("~/Documents/PhDProjects/PD1_Resistance_Signature/Results/PD1_HeatMap_PathwayActivity_Scaled.pdf")

row_EGFR <- HeatmapAnnotation(all_m[, "Response", drop =F], name="Response",
                              col = list(Response = c("R" =  "deeppink", "NR" = "darkolivegreen1")),
                              which = "row",width = unit(0.3, "cm"),show_legend = T)


heat_s_sort <- Heatmap(scale(new_order), cluster_rows = F,
                       cluster_columns = T, show_row_names = F, show_column_names = T,show_column_dend = FALSE,
                       row_title_gp = gpar(fontsize =10), combined_name_fun = NULL,
                       name="Scaled\n Pathway\nPredictions", column_title = "  ",col=my_palette, row_names_gp = gpar(fontsize = 8),
                       show_heatmap_legend = T, column_title_gp = gpar(fontsize = 12, fontface = "bold"),  width = unit(90, "mm"))


draw(heat_s_sort + row_EGFR)
dev.off()



draw(heat_s_sort + row_EGFR)










GeneList


GeneList
responsers
nonresponers






heat_s_sort <- Heatmap(scale(all_rm_or_d), cluster_rows = T,
                       cluster_columns = T, show_row_names = T, show_column_names = T,
                       row_title_gp = gpar(fontsize =10), combined_name_fun = NULL,
                       name="Pathway Predictions", column_title = "ALL",col=my_palette, row_names_gp = gpar(fontsize =6),
                       show_heatmap_legend = T, column_title_gp = gpar(fontsize = 12, fontface = "bold"),  width = unit(50, "mm"))



#pdf("~/Documents/PhDProjects/PD1_Resistance_Signature/Results/PD1_HeatMap_PathwayActivity_Names.pdf")
draw(heat_s_sort + row_EGFR)
dev.off()

draw(heat + row_EGFR)
draw(heat_s + row_EGFR)

hclust(scale(all_rm_or_d), method = "complete", members = NULL)

hc <- hclust(dist(USArrests)^2, "cen")  
memb <- cutree(hc, k = 10)

make_heatmap(probabilities_EGFR_50 , "EGFR 50")
View(probabilities_EGFR_50)


pdf("~/Documents/PhDProjects/PD1_Resistance_Signature/Results/Top50DEGenes_PD1.pdf")
draw(heat+ row_EGFR)
#dev.off()




}



