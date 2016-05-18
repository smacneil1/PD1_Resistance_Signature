######
# Runs the t.tests on BinReg predicts for the PD1 data
# Shelley Macneil
# May-16-2016
#
######

setwd("~/Documents/PhDProjects/PD1_Resistance_Signature/Results/BinReg/Probabilities/")


probs_raf_100 = read.table("probabilities_RAF_100.txt",sep='\t',header=1)
probs_raf_100[,2]
probs_raf_100=probs_raf_100[37:64, ]
row.names(probs_raf_100)=probs_raf_100[,2]

View(probs_raf_100)
row.names(probs_raf_100)=probs_raf_100[,2]

probs_akt_100 = read.table("probabilities_AKT_20.txt",sep='\t',row.names=F,header=1)
probs_egfr_50 = read.table("probabilities_EGFR_50.txt",sep='\t',row.names=F,header=1)

boxplot(probs_raf_100$Probability~probs_raf_100$)