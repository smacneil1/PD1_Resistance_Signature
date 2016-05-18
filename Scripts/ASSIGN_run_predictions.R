
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

if (length(commandArgs(trailingOnly=T)) < 1){
  print("ERROR: Submit the correct options!")
  print("Rscript ASSIGN_run_predictions.R <pathway_number>")
  quit(save = "no", status = 1, runLast = FALSE)
}

library(ASSIGN)

#----------------------------------------------------#
#Input Files (modify these locations for your system)#
#----------------------------------------------------#
working_dir      <- "~/Documents/PhDProjects/PD1_Resistance_Signature"
basedir          <- "Results"
input_rda        <- "~/Documents/PhDProjects/PD1_Resistance_Signature/Results/PD1_results.rda"
run_pathway      <- as.numeric(commandArgs(trailingOnly = T)[1])
genelists        <- "~/Documents/PhDProjects/PD1_Resistance_Signature/Data/pathway_gene_lists.rda"


#----------------------------------------------------#
#Parameters (modify these to change ASSIGN functions)#
#----------------------------------------------------#
num_genes_akt    <- 20
num_genes_bad    <- 250
num_genes_egfr   <- 15
num_genes_her2   <- 10
num_genes_igf1r  <- 125
num_genes_krasgv <- 125
num_genes_krasqh <- 150
num_genes_raf    <- 175
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

trainingLabela <- list(control=list(bad=1:12),bad=13:18)
sub_dir <- paste(basedir,paste("bad_",num_genes_bad,"_gene_list", sep=""),sep='/')
dir.create(sub_dir)
assign_easy_multi_anchor_exclude(trainingData=cbind(c_gfp,c_bad),
                                 test=c_test,
                                 trainingLabel1=trainingLabela,
                                 anchorGenes=list(bad=c("BAD")),
                                 g=NULL,
                                 geneList=list(bad=c(bad_gene_list_DOWN, bad_gene_list_UP)),
                                 out_dir_base=sub_dir,
                                 single=1,
                                 sigma_sZero = sigma_sZero,
                                 sigma_sNonZero = sigma_sNonZero,
                                 S_zeroPrior=S_zeroPrior)



if(run_pathway == 1){
  trainingLabelb <- list(control=list(bad=1:12),bad=13:18)
  sub_dir <- paste(basedir,paste("bad_",num_genes_bad,"_gene_list", sep=""),sep='/')
  dir.create(sub_dir)
  assign_easy_multi(trainingData=NULL,
                    test=c_test,
                    trainingLabel1=NULL,
                    g=NULL,
                    geneList=list(bad=c(bad_gene_list_DOWN, bad_gene_list_UP)),
                    out_dir_base=sub_dir,
                    single=1,
                    sigma_sZero = sigma_sZero,
                    sigma_sNonZero = sigma_sNonZero,
                    S_zeroPrior=S_zeroPrior)
}

if(run_pathway == 2){
  trainingLabel <- list(control=list(egfr=1:6),egfr=7:12)
  sub_dir <- paste(basedir,paste("egfr_",num_genes_egfr,"_gene_list", sep=""),sep='/')
  dir.create(sub_dir)
  assign_easy_multi(trainingData=NULL,
                    test=c_test,
                    trainingLabel1=NULL,
                    g=NULL,
                    geneList=list(egfr=egfr_gene_list[[1]]),
                    out_dir_base=sub_dir,
                    single=1,
                    sigma_sZero = sigma_sZero,
                    sigma_sNonZero = sigma_sNonZero,
                    S_zeroPrior=S_zeroPrior)
}

if(run_pathway == 3){
  trainingLabeli <- list(control=list(igf1r=1:12),igf1r=13:18)
  sub_dir <- paste(basedir,paste("igf1r_",num_genes_igf1r,"_gene_list", sep=""),sep='/')
  dir.create(sub_dir)
  assign_easy_multi(trainingData=NULL,
                    test=c_test,
                    trainingLabel1=NULL,
                    g=NULL,
                    geneList=list(igf1r=igf1r_gene_list[[1]]),
                    out_dir_base=sub_dir,
                    single=1,
                    sigma_sZero = sigma_sZero,
                    sigma_sNonZero = sigma_sNonZero,
                    S_zeroPrior=S_zeroPrior)
}

if(run_pathway == 4){
  trainingLabel <- list(control=list(krasgv=1:9),krasgv=10:18)
  sub_dir <- paste(basedir,paste("krasgv_",num_genes_krasgv,"_gene_list", sep=""),sep='/')
  dir.create(sub_dir)
  assign_easy_multi(trainingData=NULL,
                    test=c_test,
                    trainingLabel1=NULL,
                    g=NULL,
                    geneList=list(krasgv=krasgv_gene_list[[1]]),
                    out_dir_base=sub_dir,
                    single=1,
                    sigma_sZero = sigma_sZero,
                    sigma_sNonZero = sigma_sNonZero,
                    S_zeroPrior=S_zeroPrior)
}

if(run_pathway == 5){
  trainingLabel <- list(control=list(krasqh=1:9),krasqh=10:18)
  sub_dir <- paste(basedir,paste("krasqh_",num_genes_krasqh,"_gene_list", sep=""),sep='/')
  dir.create(sub_dir)
  assign_easy_multi(trainingData=NULL,
                    test=c_test,
                    trainingLabel1=NULL,
                    g=NULL,
                    geneList=list(krasqh=krasqh_gene_list[[1]]),
                    out_dir_base=sub_dir,
                    single=1,
                    sigma_sZero = sigma_sZero,
                    sigma_sNonZero = sigma_sNonZero,
                    S_zeroPrior=S_zeroPrior)
}

if(run_pathway == 6){
  trainingLabelr <- list(control=list(raf=1:12),raf=13:18)
  sub_dir <- paste(basedir,paste("raf_",num_genes_raf,"_gene_list", sep=""),sep='/')
  dir.create(sub_dir)
  assign_easy_multi(trainingData=NULL,
                    test=c_test,
                    trainingLabel1=NULL,
                    g=NULL,
                    geneList=list(raf=raf_gene_list[[1]]),
                    out_dir_base=sub_dir,
                    single=1,
                    sigma_sZero = sigma_sZero,
                    sigma_sNonZero = sigma_sNonZero,
                    S_zeroPrior=S_zeroPrior)
}




