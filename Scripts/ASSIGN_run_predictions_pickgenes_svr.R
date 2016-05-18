
# Name:    ASSIGN_run_predictions_single.R
#
# Purpose: Load a session created with ASSIGN_merge_and_combat.R and run a 
#          pathway prediction. This script runs one pathway at a time depending
#          on a number provided when running this script.
#
# Usage:   Rscript ASSIGN_run_predictions.R <pathway_number> <num_genes>
#
# Author:  David Jenkins (modified from ASSIGN scripts from Mumtahena Rahman)
# Date:    2015-09-30
#
################################################################################

if (length(commandArgs(trailingOnly=T)) < 1){
  print("ERROR: Submit the correct options!")
  print("Rscript ASSIGN_run_predictions.R <pathway_number> <num_genes>")
  quit(save = "no", status = 1, runLast = FALSE)
}

library(devtools)
install_github("wevanjohnson/ASSIGN", ref="fa6621a31629d9a122b9a2bcaa724f7d8f2c02ce")
library(ASSIGN)
source('~/Documents/PhDProjects/GFRN_signatures-master/Key_ASSIGN_functions_balancedsig.R')


#----------------------------------------------------#
#Input Files (modify these locations for your system)#
#----------------------------------------------------#

#base_dir="~/data3/Shelley/PD1Resistance"

working_dir      <- "~/Documents/PhDProjects/PD1_Resistance_Signature"
#working_dir      <- "/data3/Shelley/PD1Resistance"
basedir          <- "Results/ASSIGN"
#basedir          <- "Results/ASSIGN"

input_rda        <- "~/Documents/PhDProjects/PD1_Resistance_Signature/Results/PD1_results.rda"
input_rda          <- "~/Documents/PhDProjects/PD1_Resistance_Signature/Results/PD1_results_onestep.rda"
#input_rda        <- "~/data3/Shelley/PD1Resistance/Data/PD1_results.rda"

run_pathway         <- as.numeric(commandArgs(trailingOnly = TRUE)[1]) 
run_pathway
pathwaynum= run_pathway

num_genes           <- as.numeric(commandArgs(trailingOnly = TRUE)[2])
#print(num_genes)
gene_number = num_genes
gene_number
#----------------------------------------------------#
#Parameters (modify these to change ASSIGN functions)#
#----------------------------------------------------#
sigma_sZero    <- 0.05
sigma_sNonZero <- 0.5
S_zeroPrior    <- FALSE

#---------#
#Load Data#
#---------#
setwd(working_dir)
load(input_rda)

#-------------------------------------------#
#running ASSIGN with helper function testSig#
#-------------------------------------------#

#gene_number=350
if(pathwaynum == 1){
  trainingLabela <- list(control=list(akt=1:12),akt=13:18)
  sub_dir <- paste(basedir,paste("akt2_",gene_number,"_assign_picks", sep=""),sep='/')
  dir.create(sub_dir)
  assign_easy_multi_anchor_exclude(trainingData=cbind(c_gfp,c_akt),
                    test=c_test,
                    trainingLabel1=trainingLabela,
                    anchorGenes=list(akt=c("AKT1")),
                    g=gene_number,
                    out_dir_base=sub_dir,
                    single=1,
                    sigma_sZero = sigma_sZero,
                    sigma_sNonZero = sigma_sNonZero,
                    S_zeroPrior=S_zeroPrior)

}


if(pathwaynum == 2){
  trainingLabelb <- list(control=list(bad=1:12),bad=13:18)
  sub_dir <- paste(basedir,paste("bad_",gene_number,"_assign_picks", sep=""),sep='/')
  sub_dir
  dir.create(sub_dir)
  assign_easy_multi_anchor_exclude(trainingData=cbind(c_gfp,c_bad),
                    test=c_test,
                    trainingLabel1=trainingLabelb,
                    anchorGenes=list(bad=c("BAD")),
                    g=gene_number,
                    out_dir_base=sub_dir,
                    single=1,
                    sigma_sZero = sigma_sZero,
                    sigma_sNonZero = sigma_sNonZero,
                    S_zeroPrior=S_zeroPrior)
}

gene_number=5
if(pathwaynum == 3){
  trainingLabel <- list(control=list(egfr=1:6),egfr=7:12)
  trainingLabel
  train_egfr
  sub_dir <- paste(basedir,paste("egfr_",gene_number,"_assign_picks_onestep", sep=""),sep='/')
  sub_dir
  dir.create(sub_dir)
  assign_easy_multi_anchor_exclude(trainingData=train_egfr,
                    test=c_test,
                    trainingLabel1=trainingLabel,
                    anchorGenes=list(egfr=c("EGFR")),
                    g=gene_number,
                    out_dir_base=sub_dir,
                    single=1,
                    sigma_sZero = sigma_sZero,
                    sigma_sNonZero = sigma_sNonZero,
                    S_zeroPrior=S_zeroPrior)
}

if(pathwaynum == 4){
  trainingLabelh <- list(control=list(her2=1:12),her2=13:17)
  sub_dir <- paste(basedir,paste("her2_",gene_number,"_assign_picks", sep=""),sep='/')
  sub_dir
  dir.create(sub_dir)
  assign_easy_multi_anchor_exclude(trainingData=cbind(c_gfp,c_her2),
                    test=c_test,
                    trainingLabel1=trainingLabelh,
                    anchorGenes=list(her2=c("ERBB2")),
                    g=gene_number,
                    out_dir_base=sub_dir,
                    single=1,
                    sigma_sZero = sigma_sZero,
                    sigma_sNonZero = sigma_sNonZero,
                    S_zeroPrior=S_zeroPrior)
}

if(pathwaynum == 5){
  trainingLabeli <- list(control=list(igf1r=1:12),igf1r=13:18)
  sub_dir <- paste(basedir,paste("igf1r_",gene_number,"_assign_picks", sep=""),sep='/')
  sub_dir
   dir.create(sub_dir)
  assign_easy_multi_anchor_exclude(trainingData=cbind(c_gfp,c_igf1r),
                    test=c_test,
                    trainingLabel1=trainingLabeli,
                    anchorGenes=list(igf1r=c("IGF1R")),
                    g=gene_number,
                    out_dir_base=sub_dir,
                    single=1,
                    sigma_sZero = sigma_sZero,
                    sigma_sNonZero = sigma_sNonZero,
                    S_zeroPrior=S_zeroPrior)
}

if(pathwaynum == 6){
  trainingLabel <- list(control=list(krasgv=1:9),krasgv=10:18)
  sub_dir <- paste(basedir,paste("krasgv_",gene_number,"_assign_picks", sep=""),sep='/')
  sub_dir
  dir.create(sub_dir)
  assign_easy_multi_anchor_exclude(trainingData=cbind(c_kras_gfp,c_krasgv),
                    test=c_test,
                    trainingLabel1=trainingLabel,
                    anchorGenes=list(krasgv=c("KRAS")),
                    g=gene_number,
                    out_dir_base=sub_dir,
                    single=1,
                    sigma_sZero = sigma_sZero,
                    sigma_sNonZero = sigma_sNonZero,
                    S_zeroPrior=S_zeroPrior)
}


if(pathwaynum == 7){
  trainingLabelr <- list(control=list(raf=1:12),raf=13:18)
  sub_dir <- paste(basedir,paste("raf_", gene_number,"_assign_picks", sep=""),sep='/')
  print(sub_dir)
  dir.create(sub_dir)
  assign_easy_multi_anchor_exclude(trainingData=cbind(c_gfp,c_raf),
                    test=c_test,
                    trainingLabel1=trainingLabelr,
                    anchorGenes=list(raf=c("RAF1")),
                    g=gene_number,
                    out_dir_base=sub_dir,
                    single=1,
                    sigma_sZero = sigma_sZero,
                    sigma_sNonZero = sigma_sNonZero,
                    S_zeroPrior=S_zeroPrior)
}

# if(pathwaynum == 7){
#   trainingLabel <- list(control=list(krasqh=1:9),krasqh=10:18)
#   sub_dir <- paste(basedir,paste("krasqh_",gene_number,"_assign_picks", sep=""),sep='/')
#   sub_dir
#   dir.create(sub_dir)
#   assign_easy_multi_anchor_exclude(trainingData=cbind(c_kras_gfp,c_krasqh),
#                                    test=c_test,
#                                    trainingLabel1=trainingLabel,
#                                    anchorGenes=list(krasqh=c("KRAS")),
#                                    g=gene_number,
#                                    out_dir_base=sub_dir,
#                                    single=1,
#                                    sigma_sZero = sigma_sZero,
#                                    sigma_sNonZero = sigma_sNonZero,
#                                    S_zeroPrior=S_zeroPrior)
# }
