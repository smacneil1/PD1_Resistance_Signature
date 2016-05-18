
# Name:    ASSIGN_run_predictions.R
#
# Purpose: Load a session created with ASSIGN_merge_and_combat.R and run a 
#          pathway prediction. This script runs one pathway at a time depending
#          on a number provided when running this script.
#
# Usage:   Rscript ASSIGN_run_predictions.R <pathway_number>
#
# Author:  David Jenkins (modified from ASSIGN scripts from Mumtahena Rahman)
# Date:    2015-09-30
#
################################################################################

# if (length(commandArgs(trailingOnly=T)) < 1){
#   print("ERROR: Submit the correct options!")
#   print("Rscript ASSIGN_run_predictions.R <pathway_number>")
#   quit(save = "no", status = 1, runLast = FALSE)
# }

library(devtools)
install_github("wevanjohnson/ASSIGN", ref="fa6621a31629d9a122b9a2bcaa724f7d8f2c02ce")
library(ASSIGN)
source('~/Documents/PhDProjects/GFRN_signatures-master/Key_ASSIGN_functions_balancedsig.R')

#----------------------------------------------------#
#Input Files (modify these locations for your system)#
#----------------------------------------------------#
working_dir      <- "~/Documents/PhDProjects/PD1_Resistance_Signature"
basedir          <- "Results/ASSIGN"
input_rda        <- "~/Documents/PhDProjects/PD1_Resistance_Signature/Results/PD1_results.rda"
#run_pathway      <- as.numeric(commandArgs(trailingOnly = T)[1])
genelists        <- "~/Documents/PhDProjects/PD1_Resistance_Signature/Data/pathway_gene_lists.rda"


#----------------------------------------------------#
#Parameters (modify these to change ASSIGN functions)#
#----------------------------------------------------#
num_genes_akt    <- 20
num_genes_bad    <- 250
num_genes_egfr   <- 50
num_genes_her2   <- 10
num_genes_igf1r  <- 100
num_genes_krasgv <- 175
num_genes_krasqh <- 300
num_genes_raf    <- 350
sigma_sZero      <- 0.05
sigma_sNonZero   <- 0.5
S_zeroPrior      <- FALSE

#---------#
#Load Data#
#---------#
setwd(working_dir)
load(input_rda)
c_egfr_gfp <- train_egfr[,1:6]
c_egfr <- train_egfr[,7:12]
load(genelists)


#------------------------------#
# Single-Pathway, Optimized    #
#------------------------------#

# BAD 250
trainingLabela <- list(control=list(bad=1:12),bad=13:18)
sub_dir <- paste(basedir,paste("bad2_",num_genes_bad,"_gene_list", sep=""),sep='/')
dir.create(sub_dir)

assign_easy_multi_anchor_exclude(trainingData=cbind(c_gfp,c_bad),
                                 test=c_test,
                                 trainingLabel1=trainingLabela,
                                 anchorGenes=list(bad=c("BAD")),
                                 g=5,
                                 geneList=list(bad=c(bad_gene_list_DOWN, bad_gene_list_UP)),
                                 out_dir_base=sub_dir,
                                 single=1,
                                 sigma_sZero = sigma_sZero,
                                 sigma_sNonZero = sigma_sNonZero,
                                 S_zeroPrior=S_zeroPrior)


# AKT
  trainingLabela <- list(control=list(akt=1:12),akt=13:18)
  sub_dir <- paste(basedir,paste("akt_",num_genes_akt,"_gene_list", sep=""),sep='/')
  dir.create(sub_dir)
  assign_easy_multi_anchor_exclude(trainingData=cbind(c_gfp,c_akt),
                                   test=c_test,
                                   trainingLabel1=trainingLabela,
                                   anchorGenes=list(akt=c("AKT1")),
                                   g=NULL,
                                   geneList=list(akt=c(akt_gene_list_DOWN, akt_gene_list_UP)),
                                   out_dir_base=sub_dir,
                                   single=1,
                                   sigma_sZero = sigma_sZero,
                                   sigma_sNonZero = sigma_sNonZero,
                                   S_zeroPrior=S_zeroPrior)
  
  

# HER2
  trainingLabela <- list(control=list(her2=1:12),her2=13:17)
  sub_dir <- paste(basedir,paste("her2_",num_genes_her2,"_gene_list", sep=""),sep='/')
  dir.create(sub_dir)
  assign_easy_multi_anchor_exclude(trainingData=cbind(c_gfp,c_her2),
                                   test=c_test,
                                   trainingLabel1=trainingLabela,
                                   anchorGenes=list(her2=c("ERBB2")),
                                   g=NULL,
                                   geneList=list(her2=c(her2_gene_list_DOWN, her2_gene_list_UP)),
                                   out_dir_base=sub_dir,
                                   single=1,
                                   sigma_sZero = sigma_sZero,
                                   sigma_sNonZero = sigma_sNonZero,
                                   S_zeroPrior=S_zeroPrior)
  
# IGF1R
  trainingLabela <- list(control=list(igf1r=1:12),igf1r=13:18)
  sub_dir <- paste(basedir,paste("igf1r_",num_genes_igf1r,"_gene_list", sep=""),sep='/')
  dir.create(sub_dir)
  assign_easy_multi_anchor_exclude(trainingData=cbind(c_gfp,c_igf1r),
                                   test=c_test,
                                   trainingLabel1=trainingLabela,
                                   anchorGenes=list(igf1r=c("IGF1R")),
                                   g=NULL,
                                   geneList=list(igf1r=c(igf1r_gene_list_DOWN, igf1r_gene_list_UP)),
                                   out_dir_base=sub_dir,
                                   single=1,
                                   sigma_sZero = sigma_sZero,
                                   sigma_sNonZero = sigma_sNonZero,
                                   S_zeroPrior=S_zeroPrior)

  # KRAS GV
  trainingLabela <- list(control=list(krasgv=1:9),krasgv=10:18)
  sub_dir <- paste(basedir,paste("krasgv_",num_genes_krasgv,"_gene_list", sep=""),sep='/')
  dir.create(sub_dir)
  assign_easy_multi_anchor_exclude(trainingData=cbind(c_gfp,c_krasgv),
                                   test=c_test,
                                   trainingLabel1=trainingLabela,
                                   anchorGenes=list(krasgv=c("KRAS")),
                                   g=NULL,
                                   geneList=list(krasgv=c(krasgv_gene_list_DOWN, krasgv_gene_list_UP)),
                                   out_dir_base=sub_dir,
                                   single=1,
                                   sigma_sZero = sigma_sZero,
                                   sigma_sNonZero = sigma_sNonZero,
                                   S_zeroPrior=S_zeroPrior)
  
# RAF
  trainingLabela <- list(control=list(raf=1:12),raf=13:18)
  sub_dir <- paste(basedir,paste("raf_",num_genes_raf,"_gene_list", sep=""),sep='/')
  dir.create(sub_dir)
  
  assign_easy_multi_anchor_exclude(trainingData=cbind(c_gfp,c_raf),
                                   test=c_test,
                                   trainingLabel1=trainingLabela,
                                   anchorGenes=list(raf=c("RAF1")),
                                   g=NULL,
                                   geneList=list(raf=c(raf_gene_list_DOWN, raf_gene_list_UP)),
                                   out_dir_base=sub_dir,
                                   single=1,
                                   sigma_sZero = sigma_sZero,
                                   sigma_sNonZero = sigma_sNonZero,
                                   S_zeroPrior=S_zeroPrior)
  
  

#EGFR
  trainingLabela <- list(control=list(egfr=1:6),egfr=7:12)
  sub_dir <- paste(basedir,paste("egfr_",num_genes_egfr,"_gene_list", sep=""),sep='/')
  dir.create(sub_dir) 
  
  assign_easy_multi_anchor_exclude(trainingData=cbind(c_gfp,c_egfr),
                                   test=c_test,
                                   trainingLabel1=trainingLabela,
                                   anchorGenes=list(egfr=c("EGFR")),
                                   g=NULL,
                                   geneList=list(egfr=c(egfr_gene_list_DOWN, egfr_gene_list_UP)),
                                   out_dir_base=sub_dir,
                                   single=1,
                                   sigma_sZero = sigma_sZero,
                                   sigma_sNonZero = sigma_sNonZero,
                                   S_zeroPrior=S_zeroPrior)
  
  

  #KRAS QH
  trainingLabela <- list(control=list(krasqh=1:9),krasqh=10:18)
  sub_dir <- paste(basedir,paste("krasqh_",num_genes_krasqh,"_gene_list", sep=""),sep='/')
  dir.create(sub_dir)
  
  
  assign_easy_multi_anchor_exclude(trainingData=cbind(c_gfp,c_krasqh),
                                   test=c_test,
                                   trainingLabel1=trainingLabela,
                                   anchorGenes=list(krasqh=c("KRAS")),
                                   g=NULL,
                                   geneList=list(krasqh=c(krasqh_gene_list_DOWN, krasqh_gene_list_UP)),
                                   out_dir_base=sub_dir,
                                   single=1,
                                   sigma_sZero = sigma_sZero,
                                   sigma_sNonZero = sigma_sNonZero,
                                   S_zeroPrior=S_zeroPrior)
  
  

