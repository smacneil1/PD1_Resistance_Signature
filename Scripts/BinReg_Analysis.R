######
# Runs the t.tests on BinReg predicts for the PD1 data
# Shelley Macneil
# May-16-2016
#
######


library(ComplexHeatmap)
source('~/Documents/PhDProjects/GFRN_signatures/Key_ASSIGN_functions_balancedsig.R')
setwd("~/Documents/PhDProjects/PD1_Resistance_Signature/Results/BinReg/Probabilities/")


col_breaks = c(seq(0,0.2,length=100), seq(0.2,0.4,length=100), seq(0.4,1,length=100))

files=system("ls", intern=TRUE)
files
response <- read.table("~/Documents/PhDProjects/PD1_Resistance_Signature/Data/GSE78220_Response.txt",row.names=1, header = T)
dim(response)

# function to read in a list of filenames and store as a variable (add to list)
for(i in 1:length(files)){
  ###reading in the filess one at a time
  f<-read.table(files[i], header=1,sep='\t')
  var_name=gsub('\\..*$','',files[i])
  print(var_name)
  assign(var_name, f) 
}  
length(files)
dim(probabilities_KRAS_90)

dim(probabilities_EGFR_50)
dim(probabilities_EGFR_100)
dim(probabilities_EGFR_300)



# make the boxplot function

plotBinReg= function (probs, main){
  rownames(probs)=probs[,1]
  probs=merge_drop(probs, response, by=0)
  probs=probs[order(probs$Response),] 
  par(cex.lab=1)
  par(cex.axis=0.75)
  margins=c(8,8)
  boxplot(probs$Probability~probs$Response, main = main, ylab="Pathway Activity", col=c("darkolivegreen1","deeppink" ), names=c("Non-Responders", "Responders"))
  ttest=t.test(probs$Probability~probs$Response)
  title(sub=paste("P-value = ", format(ttest$p.value, digits=3)), adj=1, line=3, font=1)
  return(probs)
}

# ----------------------------- KRAS

# KRAS 300 

# stable from 100-50, go with 50? 
probabilities_KRAS_300=probabilities_KRAS_300[37:64,c(2,5) ]
probabilities_KRAS_100=probabilities_KRAS_100[37:64,c(2,5) ] #one we used
probabilities_KRAS_90=probabilities_KRAS_90[37:64,c(2,5) ] 
probabilities_KRAS_80=probabilities_KRAS_80[37:64,c(2,5) ] 
#probabilities_KRAS_250=probabilities_KRAS_250[37:64,c(2,5) ]
probabilities_KRAS_60=probabilities_KRAS_60[37:64,c(2,5) ]
probabilities_KRAS_50=probabilities_KRAS_50[37:64,c(2,5) ]
probabilities_KRAS_75=probabilities_KRAS_75[37:64,c(2,5) ]
probabilities_KRAS_25=probabilities_KRAS_25[37:64,c(2,5) ]

pdf("~/Documents/PhDProjects/PD1_Resistance_Signature/Results/Figures/PD1_Boxplots_KRAS.pdf")
par(mfrow = c(2, 2))
plotBinReg(probabilities_KRAS_25, "KRAS Pathway 25")  #0.47
plotBinReg(probabilities_KRAS_50, "KRAS Pathway 50") # 0.17 , #0.192
plotBinReg(probabilities_KRAS_60, "KRAS Pathway 60") #0.227
plotBinReg(probabilities_KRAS_75, "KRAS Pathway 75") #0.118
par(mfrow = c(2, 2))
plotBinReg(probabilities_KRAS_80, "KRAS Pathway 80") # 0.106
plotBinReg(probabilities_KRAS_90, "KRAS Pathway 90") # 0.0987
plotBinReg(probabilities_KRAS_100, "KRAS Pathway 100") # 0.15   # go with this one #0.11
#plotBinReg(probabilities_KRAS_250, "KRAS Pathway 250") # 0.34
plotBinReg(probabilities_KRAS_300, "KRAS Pathway 300") # 0.38
dev.off()

#------------------------------ RAF

probabilities_RAF_10=probabilities_RAF_10[37:64, c(2,5) ]
probabilities_RAF_20=probabilities_RAF_20[37:64, c(2,5) ]
probabilities_RAF_40=probabilities_RAF_40[37:64, c(2,5) ]
probabilities_RAF_50=probabilities_RAF_50[37:64, c(2,5) ]
probabilities_RAF_75=probabilities_RAF_75[37:64, c(2,5) ]
probabilities_RAF_100=probabilities_RAF_100[37:64,c(2,5) ]
#probabilities_RAF_175=probabilities_RAF_175[37:64,c(2,5) ]
probabilities_RAF_200=probabilities_RAF_200[37:64, c(2,5) ]
probabilities_RAF_250=probabilities_RAF_250[37:64, c(2,5) ]
probabilities_RAF_275=probabilities_RAF_275[37:64, c(2,5) ]
probabilities_RAF_300=probabilities_RAF_300[37:64,c(2,5) ]
probabilities_RAF_350=probabilities_RAF_350[37:64,c(2,5) ]
probabilities_RAF_500=probabilities_RAF_500[37:64, c(2,5)]
probabilities_RAF_1000=probabilities_RAF_1000[37:64, c(2,5)]

# I dont know if I trust it it just stays stable as you get higher and higher.
# but if we want to use prob use 200 ot 275
plotBinReg(probabilities_RAF_10, "RAF Pathway 10") # 0.236
plotBinReg(probabilities_RAF_20, "RAF Pathway 20") # 0.726
plotBinReg(probabilities_RAF_40, "RAF Pathway 40") # 0.633
plotBinReg(probabilities_RAF_50, "RAF Pathway 50")# 0.87, now #0.719
plotBinReg(probabilities_RAF_75, "RAF Pathway 75") #0.346
plotBinReg(probabilities_RAF_100, "RAF Pathway 100") # 0.32, now #0.265
#plotBinReg(probabilities_RAF_175, "RAF Pathway") ##0.27
plotBinReg(probabilities_RAF_200, "RAF Pathway 200") #0.193
plotBinReg(probabilities_RAF_250, "RAF Pathway 250") #0.12
plotBinReg(probabilities_RAF_275, "RAF Pathway 275") #0.113
plotBinReg(probabilities_RAF_300, "RAF Pathway 300") #  0.14, now 0.125
plotBinReg(probabilities_RAF_350, "RAF Pathway 350") #0.13   #  go with this one 
plotBinReg(probabilities_RAF_500, "RAF Pathway") #0.17, not 0.161
plotBinReg(probabilities_RAF_1000, "RAF Pathway") #0.161

#------------------------------ EGFR

#can go as low as 20 because 20,25,and 50 are all the same
probabilities_EGFR_10=probabilities_EGFR_10[25:52,c(2,5)] 
probabilities_EGFR_15=probabilities_EGFR_15[25:52,c(2,5)] 
probabilities_EGFR_20=probabilities_EGFR_20[25:52,c(2,5)] 
probabilities_EGFR_25=probabilities_EGFR_25[25:52,c(2,5)]  
probabilities_EGFR_50=probabilities_EGFR_50[25:52,c(2,5)]
probabilities_EGFR_100=probabilities_EGFR_100[25:52,c(2,5)]  
probabilities_EGFR_300=probabilities_EGFR_300[25:52,c(2,5)]  
dim(response)
pdf("~/Documents/PhDProjects/PD1_Resistance_Signature/Results/Figures/PD1_Boxplots_EGFR.pdf")

par(mfrow = c(2, 2))
plotBinReg(probabilities_EGFR_10, "EGFR Pathway 10" ) #0.006
plotBinReg(probabilities_EGFR_15, "EGFR Pathway 15" ) #0. 0244
plotBinReg(probabilities_EGFR_20, "EGFR Pathway 20" ) #0.0431
plotBinReg(probabilities_EGFR_25, "EGFR Pathway 25" ) # #0.044 before, 0.0321
par(mfrow = c(2, 2))
plotBinReg(probabilities_EGFR_50, "EGFR Pathway 50"  ) #0.044  # go with this one #0.057
plotBinReg(probabilities_EGFR_100, "EGFR Pathway 100" ) # 0.069 before
plotBinReg(probabilities_EGFR_300, "EGFR Pathway 300" ) #0.057

dev.off()
# ------------------------------ IGFR1


#probabilities_IGF1R_20=probabilities_IGF1R_15[37:64,c(2,5)]
probabilities_IGF1R_20=probabilities_IGF1R_20[37:64,c(2,5)]
probabilities_IGF1R_25=probabilities_IGF1R_25[37:64,c(2,5)]
probabilities_IGF1R_25
probabilities_IGF1R_40=probabilities_IGF1R_40[37:64,c(2,5)]
probabilities_IGF1R_50=probabilities_IGF1R_50[37:64,c(2,5)]
probabilities_IGF1R_50
probabilities_IGF1R_60=probabilities_IGF1R_60[37:64,c(2,5)]
probabilities_IGF1R_75=probabilities_IGF1R_75[37:64,c(2,5)]
probabilities_IGF1R_75
probabilities_IGF1R_80=probabilities_IGF1R_80[37:64,c(2,5)]
probabilities_IGF1R_80
probabilities_IGF1R_85=probabilities_IGF1R_85[37:64,c(2,5)]
probabilities_IGF1R_90=probabilities_IGF1R_90[37:64,c(2,5)]
probabilities_IGF1R_90
probabilities_IGF1R_100=probabilities_IGF1R_100[37:64,c(2,5)]
probabilities_IGF1R_300=probabilities_IGF1R_300[37:64,c(2,5)]
probabilities_IGF1R_500=probabilities_IGF1R_500[37:64,c(2,5)]


pdf("~/Documents/PhDProjects/PD1_Resistance_Signature/Results/Figures/PD1_Boxplots_IGF1R.pdf")
#plotBinReg(probabilities_IGF1R_15, "IGF1R Pathway: 15 genes") # dont trust it 
par(mfrow = c(2, 2))
plotBinReg(probabilities_IGF1R_20, "IGF1R Pathway: 20 genes") # 0.187 #not sure I trust this one
plotBinReg(probabilities_IGF1R_25, "IGF1R Pathway: 25 genes") # 0.38 before, 0.63 now
plotBinReg(probabilities_IGF1R_40, "IGF1R Pathway: 40 genes") # 0.54
plotBinReg(probabilities_IGF1R_50, "IGF1R Pathway: 50 genes") # 0.38 before, 0.389 now

par(mfrow = c(2, 2))
plotBinReg(probabilities_IGF1R_60, "IGF1R Pathway: 60 genes") # 0.38 before, 0.193 now
plotBinReg(probabilities_IGF1R_75, "IGF1R Pathway: 75 genes") #0.07
# maybe go with 80 or 85
plotBinReg(probabilities_IGF1R_80, "IGF1R Pathway: 80 genes") #0.0415
plotBinReg(probabilities_IGF1R_85, "IGF1R Pathway: 85 genes") #0.0285
par(mfrow = c(2, 2))
plotBinReg(probabilities_IGF1R_90, "IGF1R Pathway: 90 genes") #0.0152 before, 0.01
plotBinReg(probabilities_IGF1R_100, "IGF1R Pathway: 100") # 0.009   #go with this one (before), 0.009
plotBinReg(probabilities_IGF1R_300, "IGF1R Pathway: 300") #  0.01 before, 0.018
plotBinReg(probabilities_IGF1R_500, "IGF1R Pathway: 300") #  0.03

dev.off()

#----------------- AKT

probabilities_AKT_5=probabilities_AKT_5[37:64,c(2,5)]
probabilities_AKT_10=probabilities_AKT_10[37:64,c(2,5)]
probabilities_AKT_15=probabilities_AKT_15[37:64,c(2,5)]
probabilities_AKT_16=probabilities_AKT_16[37:64,c(2,5)]
probabilities_AKT_20=probabilities_AKT_20[37:64,c(2,5)]
probabilities_AKT_25=probabilities_AKT_25[37:64,c(2,5)]
probabilities_AKT_30=probabilities_AKT_30[37:64,c(2,5)]
probabilities_AKT_40=probabilities_AKT_40[37:64,c(2,5)]
probabilities_AKT_50=probabilities_AKT_50[37:64,c(2,5)]
probabilities_AKT_75=probabilities_AKT_75[37:64,c(2,5)]
probabilities_AKT_100=probabilities_AKT_100[37:64,c(2,5)]

# I think its more stable with 10-20 genes, say we go with 10. 
pdf("~/Documents/PhDProjects/PD1_Resistance_Signature/Results/Figures/PD1_Boxplots_AKT.pdf")
#plotBinReg(probabilities_AKT_5, "AKT Pathway") #0.91 before, # 0.99 now
par(mfrow = c(2, 2))
plotBinReg(probabilities_AKT_10, "AKT Pathway:10 genes")  # 0.16 before, # 0.164 now
plotBinReg(probabilities_AKT_15, "AKT Pathway:15 genes")  # 0.21
#plotBinReg(probabilities_AKT_16, "AKT Pathway")  #0.209
plotBinReg(probabilities_AKT_20, "AKT Pathway: 20 genes") # 0.18 # go with this one (before), #0.187
plotBinReg(probabilities_AKT_25, "AKT Pathway: 25 genes") #0.63

# not good 
par(mfrow = c(2, 2))
#plotBinReg(probabilities_AKT_30, "AKT Pathway ") #0.76
#plotBinReg(probabilities_AKT_40, "AKT Pathway")  #0.75
plotBinReg(probabilities_AKT_50, "AKT Pathway: 50") # 0.7 before
#plotBinReg(probabilities_AKT_75, "AKT Pathway")  #0.83
#par(mfrow = c(2, 2))
plotBinReg(probabilities_AKT_100, "AKT Pathway: 100") # 0.99 before, 0.95
dev.off()

#probabilities_AKT_4=probabilities_AKT_4[37:64,c(2,5)]

#-------------HER

# none really worked, dont use it. 
probabilities_HER2_10=probabilities_HER2_10[35:62,c(2,5) ]
probabilities_HER2_15=probabilities_HER2_15[35:62,c(2,5) ]
probabilities_HER2_20=probabilities_HER2_20[35:62,c(2,5) ]
probabilities_HER2_50=probabilities_HER2_50[35:62,c(2,5) ]
probabilities_HER2_100=probabilities_HER2_100[35:62,c(2,5) ]
probabilities_HER2_200=probabilities_HER2_200[35:62,c(2,5) ]


plotBinReg(probabilities_HER2_10, "HER2 Pathway 10") # 0.875
plotBinReg(probabilities_HER2_15, "HER2 Pathway 15") # 0.61
plotBinReg(probabilities_HER2_20, "HER2 Pathway 20") # 0.957
#plotBinReg(probabilities_HER2_25) # 0.89
plotBinReg(probabilities_HER2_50, "HER2 Pathway 50") # 0.458
plotBinReg(probabilities_HER2_100, "HER2 Pathway 100") #0.529
plotBinReg(probabilities_HER2_200, "HER2 Pathway 200") # 0.362

# ---------------BAD

# BAD didn't work dont use this one

probabilities_BAD_20=probabilities_BAD_20[37:64,c(2,5) ]
probabilities_BAD_30=probabilities_BAD_30[37:64,c(2,5) ]
probabilities_BAD_100=probabilities_BAD_100[37:64,c(2,5) ]
probabilities_BAD_200=probabilities_BAD_200[37:64,c(2,5) ]

#probabilities_BAD_250=probabilities_BAD_250[37:64,c(2,5) ]
#probabilities_BAD_350=probabilities_BAD_350[37:64,c(2,5) ]
#probabilities_BAD_50=probabilities_BAD_50[37:64,c(2,5) ]


plotBinReg(probabilities_BAD_20, "BAD Pathway 20") #0.342
plotBinReg(probabilities_BAD_30, "BAD Pathway 30") #0.502
plotBinReg(probabilities_BAD_100, "BAD Pathway 100") # #0.89, 0.889
plotBinReg(probabilities_BAD_200, "BAD Pathway 200") # 0.881

#plotBinReg(probabilities_BAD_250) #
#plotBinReg(probabilities_BAD_350) #
#plotBinReg(probabilities_BAD_50) #  0.64

pdf("~/Documents/PhDProjects/PD1_Resistance_Signature/Results/PD1_Boxplots_Andrea.pdf")

# final ones from paper
plotBinReg(probabilities_AKT_20, "AKT Pathway")
#plotBinReg(probabilities_RAF_350, "RAF Pathway")
plotBinReg(probabilities_IGF1R_100, "IGF1R Pathway")
plotBinReg(probabilities_EGFR_50, "EGFR Pathway")
#plotBinReg(probabilities_KRAS_100, "KRAS Pathway")

dev.off()

# plot and get the p.val from the ttest 

pdf("~/Documents/PhDProjects/PD1_Resistance_Signature/Results/Boxplots.pdf")



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

#all= cbind(probabilities_EGFR_50,probabilities_AKT_20,probabilities_IGF1R_100,probabilities_KRAS_100,probabilities_HER2_10, probabilities_BAD_100, probabilities_RAF_100   )
#colnames(all) = c("EGFR", "AKT", "IGF1R", "KRAS", "HER2", "BAD")

all= cbind(probabilities_EGFR_50,probabilities_AKT_20,probabilities_IGF1R_100  )
colnames(all) = c("EGFR", "AKT", "IGF1R")

View(all)
# function for ploting the responds vs. nonresponders



                                    ## make the heatmaps 
response
all_m=cbind(all, response)
View(all_m)
all_rm=all_m[,-4]
all_rm
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
my_palette <- colorRampPalette(c("darkblue","aliceblue","brown4"))(n = 299)
pdf("~/Documents/PhDProjects/PD1_Resistance_Signature/Results/PD1_HeatMap_PathwayActivity_Scaled_Gaju.pdf")

View(all_m)
row_EGFR <- HeatmapAnnotation(all_m[, "Response", drop =F], name="Response",
                              col = list(Response = c("R" =  "deeppink", "NR" = "darkolivegreen1")),
                              which = "row",width = unit(0.3, "cm"),show_legend = T)


View(new_order)
View(response)
heat_s_sort <- Heatmap(new_order, cluster_rows = F,
                       cluster_columns = T, show_row_names = F, show_column_names = T,show_column_dend = FALSE,
                       row_title_gp = gpar(fontsize =10), combined_name_fun = NULL,
                       name="Scaled\n Pathway\nPredictions", column_title = "  ",col=my_palette, row_names_gp = gpar(fontsize = 8),
                       show_heatmap_legend = T, column_title_gp = gpar(fontsize = 12, fontface = "bold"),  width = unit(90, "mm"))


draw(heat_s_sort + row_EGFR)
dev.off()



draw(heat_s_sort + row_EGFR)




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



