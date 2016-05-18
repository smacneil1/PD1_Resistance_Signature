
library(devtools)
source("https://bioconductor.org/biocLite.R")
biocLite("BiocInstaller")
install_github("wevanjohnson/sva-devel", ref="2e87233c3797c3a630d3cf662e2ee08dfe42c7b5")
library(sva)
library(data.table)

# Name:    ASSIGN_merge_and_combat.R
#
# Purpose: Merge together the signature data, and a test dataset and perform
#          the reference version of ComBat and save a session that can be used
#          with ASSIGN_run_predictions_single.R to run ASSIGN
#
# Usage:   Rscript ASSIGN_merge_and_combat.R
#
# Author:  David Jenkins (modified from ASSIGN scripts from Mumtahena Rahman)
# Date:    2015-09-29
#
################################################################################

#----------------------------------------------------#
#Input Files (modify these locations for your system)#
#----------------------------------------------------#
signatures_dir      <- "~/Documents/PhDProjects/SignatureData/"
expr_file           <- paste(signatures_dir,"GFP18_AKT_BAD_HER2_IGF1R_RAF_ERK.tpmlog",sep="/")
control_egfr_l_file <- paste(signatures_dir,"18_GFP_EGFR_TPMlog2.txt",sep="/")
gfp_kras_file       <- paste(signatures_dir,"GFP30_KRAS-GV_KRAS-QH_KRAS-WT_tpmlog.txt",sep="/")
key_assign_file     <- "~/Documents/PhDProjects/GFRN_signatures-master/Key_ASSIGN_functions_balancedsig.R"
testFile            <- "~/Documents/PhDProjects/PD1_Resistance_Signature/Data/GSE78220_PatientFPKM_withResponse_Sorted.txt"

#--------------------------------------#
#Output Files (modify these every time)#
#--------------------------------------#
working_dir         <- "~/Documents/PhDProjects/PD1_Resistance_Signature/"
output_rda          <- "~/Documents/PhDProjects/PD1_Resistance_Signature/Results/PD1_results_onestep.rda"

#---------#
#Load Data#
#---------#
source(key_assign_file)
setwd(working_dir)


#load in test data 
test_unlogged<-read.table(testFile,row.names=1,header=TRUE)
test_unlogged_1= (test_unlogged +1)
test=log2(test_unlogged_1)
#View(test)
dim(test)
expr_test <-test[apply(test[,1:28]==0,1,mean) < 0.85,]
dim(expr_test)
View(expr_test)
pcaplot(expr_test, sub = rep(1,28))
write.table(expr_test, "~/Documents/PhDProjects/PD1_Resistance_Signature/Data/GSE78220_FPKMlog2_filtered.txt", sep='\t',quote=F, col.names = NA, row.names= T)
## EGFR
control_egfr_l<-read.table(control_egfr_l_file, sep='\t', header=1, row.names=1)
colnames(control_egfr_l)
sub<-c(6,6)
dim(control_egfr_l)
control_egfr_l <-control_egfr_l[apply(control_egfr_l[,1:12]==0,1,mean) < 0.85,]
dim(control_egfr_l)
egfr_gfp <- control_egfr_l[,1:6]
head(egfr )
egfr <- control_egfr_l[,7:12]
head(egfr_gfp)
write.table(egfr, "~/Documents/PhDProjects/PD1_Resistance_Signature/Data/Signature_egfr.txt", sep='\t',quote=F, col.names = NA, row.names= T)
write.table(egfr_gfp, "~/Documents/PhDProjects/PD1_Resistance_Signature/Data/Signature_gfp_egfr.txt", sep='\t',quote=F, col.names = NA, row.names= T)
colnames(control_egfr_l)
pcaplot(control_egfr_l, sub)
# merge EGFR with test data
control_egfr_l_test=merge_drop(control_egfr_l, expr_test)
View(control_egfr_l_test)
#gfp_egfr_multi_f <- merge_drop(control_egfr_l,expr_all_f)
sub<-c(6,6,ncol(test))
sub
pcaplot(control_egfr_l_test,sub)
bat<-c(rep(1,12),rep(2,28))
bat
#colnames(bat)<-c("Batch","Model")
rownames(bat)<-colnames(control_egfr_l_test)
#mod <- model.matrix(~as.factor(bat$Model))
#  mod
combat_control_egfr_l_test<-ComBat(dat=control_egfr_l_test, batch=bat, mod=NULL, ref.batch=1) 
colnames(combat_control_egfr_l_test)
dim(combat_control_egfr_l_test)
train_egfr <- combat_control_egfr_l_test[,1:12]
View(train_egfr)
c_egfr <- combat_control_egfr_l_test[,7:12]
c_test <- combat_control_egfr_l_test[,13:40]


save.image(output_rda)


#KRAS alone

gfp_kras<-read.table(gfp_kras_file, sep='\t', header=1, row.names=1)
dim(gfp_kras)
#gfp_egfr_kras_multi_f<-merge_drop(gfp_egfr_multi_f,gfp_kras)
#sapply(gfp_egfr_kras_multi_f, class)
sum(sub)
colnames(gfp_kras)
dim(gfp_kras)
dim(gfp_kras)

expr_kras  <-gfp_kras[apply(gfp_kras[,1:35]==0,1,mean) < 0.85,]
gfp_kras=expr_kras[,1:9]
head(gfp_kras)
krasgv=expr_kras[,10:18]
head(krasgv)
krasqh=expr_kras[,19:27]
head(krasqh)

write.table(gfp_kras, "~/Documents/PhDProjects/PD1_Resistance_Signature/Data/Signature_gfp_kras.txt", sep='\t',quote=F, col.names = NA, row.names= T)
write.table(krasgv, "~/Documents/PhDProjects/PD1_Resistance_Signature/Data/Signature_krasgv.txt", sep='\t',quote=F, col.names = NA, row.names= T)
write.table(krasqh, "~/Documents/PhDProjects/PD1_Resistance_Signature/Data/Signature_krasqh.txt", sep='\t',quote=F, col.names = NA, row.names= T)




sub<-c(9,9,9)
pcaplot(expr_kras, sub)
View(gfp_kras)
expr_kras_test=merge_drop(expr_kras, expr_test)
View(expr_kras_test)
sub<-c(9,9,9, ncol(test))
sub
pcaplot(expr_kras_test,sub)
dim(expr_kras_test)
bat<-as.data.frame(cbind(c(rep(1,27),rep(2,28)),c(rep(1,9),rep(2,9),rep(2,9), rep(1,28))))
colnames(bat)<-c("Batch","Model")
bat
rownames(bat)<-colnames(expr_kras_test)
mod <- model.matrix(~as.factor(bat$Model))
mod
combat_control_expr_kras_test<-ComBat(dat=expr_kras_test batch=bat[,1], mod=mod, ref.batch=1) 
colnames(combat_control_expr_kras_test)
dim(combat_control_expr_kras_test)

c_krasgv <- combat_control_expr_kras_test[,10:18]
c_kras_gfp <- combat_control_expr_kras_test[,1:9]
c_test <- combat_control_expr_kras_test[,28:55]


# all others
expr<-as.matrix(read.table(expr_file,sep='\t',row.names=1,header=1))
dim(expr)
expr_all_f <-expr_all[apply(expr_all[,1:41]==0,1,mean) < 0.85,]
dim(expr_all_f)

control_raf_her2_akt_bad_igf1r<-subset(expr_all_f , select=GFP.1:GFP.12)
raf<-subset(expr_all_f ,select=RAF.1:RAF.6)
her2<-subset(expr_all_f , select=HER2.1:HER2.6)
akt<-subset(expr_all_f ,select=AKT.1:AKT.6)
bad<-subset(expr_all_f ,select=BAD.1:BAD.6)

igf1r<-subset(expr_all_f ,select=IGF1R.1:IGF1R.6)

write.table(control_raf_her2_akt_bad_igf1r, "~/Documents/PhDProjects/PD1_Resistance_Signature/Data/Signature_control_batch1.txt", sep='\t',quote=F, col.names = NA, row.names= T)
write.table(raf, "~/Documents/PhDProjects/PD1_Resistance_Signature/Data/Signature_raf.txt", sep='\t',quote=F, col.names = NA, row.names= T)
write.table(her2, "~/Documents/PhDProjects/PD1_Resistance_Signature/Data/Signature_her2.txt", sep='\t',quote=F, col.names = NA, row.names= T)
write.table(akt, "~/Documents/PhDProjects/PD1_Resistance_Signature/Data/Signature_akt.txt", sep='\t',quote=F, col.names = NA, row.names= T)
write.table(bad, "~/Documents/PhDProjects/PD1_Resistance_Signature/Data/Signature_bad.txt", sep='\t',quote=F, col.names = NA, row.names= T)
write.table(igf1r, "~/Documents/PhDProjects/PD1_Resistance_Signature/Data/Signature_igf1r.txt", sep='\t',quote=F, col.names = NA, row.names= T)




expr_all<-cbind(control,akt,bad,her2,igf1r,raf)














#------#
#ComBat#
#------#






#### Before, doing them all together. 

pdf("~/Documents/PhDProjects/PD1_Resistance_Signature/Results/PD1_batch_twostep_postlog.pdf")
pcaplot(gfp_egfr_kras_multi_f,sub)
colnames(gfp_egfr_kras_multi_f)
bat<-as.data.frame(cbind(c(rep(1,12),rep(2,41),rep(3,36)),c(rep(1,6),rep(2,6),rep(1,12),rep(3,6),rep(4,6),rep(5,5),rep(6,6),rep(7,6),rep(1,9),rep(8,9),rep(9,9),rep(10,9))))
bat
colnames(bat)<-c("Batch","Model")
rownames(bat)<-colnames(gfp_egfr_kras_multi_f)
mod <- model.matrix(~as.factor(bat$Model))
mod
bat

combat_expr<-ComBat(dat=gfp_egfr_kras_multi_f, batch=bat[,1], mod=mod, ref.batch=2)
pcaplot(combat_expr,sub)
# combat with test data
dat<-merge_drop(combat_expr,test)
dim(dat)
ncol(test)

sub<-c(6,6,12,6,6,5,6,6,9,9,9,9,ncol(test))
dim(dat)
colnames(dat)
pcaplot(dat,sub)
bat<-as.matrix(cbind(colnames(dat),c(rep(1,ncol(gfp_egfr_kras_multi_f)),rep(2,ncol(test)))))
bat
combat_expr1<-ComBat(dat=dat,batch=bat[,2], mod=NULL, ref.batch=1)
pcaplot(combat_expr1,sub)
dev.off()
c_gfp<-subset(combat_expr1, select=GFP.1:GFP.12)
c_gfp
c_akt<-subset(combat_expr1, select=AKT.1:AKT.6)
c_bad<-subset(combat_expr1, select=BAD.1:BAD.6)
c_her2<-subset(combat_expr1, select=HER2.1:HER2.6)
c_igf1r<-subset(combat_expr1, select=IGF1R.1:IGF1R.6)
c_raf<-subset(combat_expr1, select=RAF.1:RAF.6)
train_egfr<-combat_expr1[,1:12]
c_egfr_gfp <- train_egfr[,1:6]
c_egfr <- train_egfr[,7:12]
c_kras_gfp<-subset(combat_expr1,select=GFP30.1:GFP30.9)
c_kraswt<-subset(combat_expr1,select=KRAS_WT.1:KRAS_WT.9)
c_krasqh<-subset(combat_expr1,select=KRAS_QH.1:KRAS_QH.9)
c_krasgv<-subset(combat_expr1,select=KRAS_GV.1:KRAS_GV.9)
c_test<-combat_expr1[,(ncol(gfp_egfr_kras_multi_f)+1):ncol(combat_expr1)]
View(c_test)

# test files
GSE78220_Response <- read.table("~/Documents/PhDProjects/PD1_Resistance_Signature/Data/GSE78220_Response.txt",row.names=1, header = T)
GSE78220_Response
test_t=t(test)
View(test)
test_resp=cbind(GSE78220_Response, test_t)
test_resp_t=t(test_resp)

responders <-test_resp_t[,1:15]
responders=responders[-1,]
nonresponders<-test_resp_t[,16:28]
nonresponders=nonresponders[-1,]
write.table(responders, "~/Documents/PhDProjects/PD1_Resistance_Signature/Data/GSE78220_PatientFPKMlog2_Responders.txt", sep='\t',quote=F, col.names = NA, row.names= T)
write.table(nonresponders, "~/Documents/PhDProjects/PD1_Resistance_Signature/Data/GSE78220_PatientFPKMlog2_NonResponders.txt", sep='\t',quote=F, col.names = NA, row.names= T)




save.image(file=output_rda)
