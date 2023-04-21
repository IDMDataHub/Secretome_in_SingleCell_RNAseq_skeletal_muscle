#!/usr/bin/env R

###>Author: Maria Kondili
###>Date : May 2022
###>Subject: Analyse Outcyte Results of Oprescu-FAPs-markers,
###>         Add GeneName in results with EntrezID

library(tidyverse)
library(biomaRt)
library(stringr)
library(httr)
library(seqinr)


work_dir <- "Secretome_in_Single_Cell_RNAseq_skeletal_muscle/"
secretome_dir <- paste0(work_dir,"Oprescu_mouse/Secretome/")
outcyte_dir <- paste0(secretome_dir ,"Outcyte/")

source(paste0(work_dir,"functions/convert_GeneIds.R")


###> All FAPs:


outcyte_faps <- readr::read_delim(paste0(outcyte_dir,"Results/PepSeq_FAPsMarkers_Oprescu_pipeline.txt"),col_names=F)
colnames(outcyte_faps) <- c("EntrezID","Prediction","Score")

outcyte_faps$GeneSymbol  <- convert_gene_ids.mm(id_list=as.character(outcyte_faps$EntrezID),from="ENTREZID",to="SYMBOL")

write_delim(outcyte_faps,
            paste0(outcyte_dir,"Results/Outcyte_predictions_FAPsMarkers_with_GeneName.tsv"),
            delim="\t",col_names = T)

###> Fibroblasts :

outcyte_fibro <- readr::read_delim(paste0(outcyte_dir,"Results/PepSeq_FibroMarkers_Oprescu_pipeline.txt"),col_names=F)
colnames(outcyte_fibro) <- c("EntrezID","Prediction","Score")

outcyte_fibro$GeneSymbol  <- convert_gene_ids.mm(as.character(outcyte_fibro$EntrezID),from="ENTREZID",to="SYMBOL")

write_delim(outcyte_fibro,
            paste0(outcyte_dir,"Results/Outcyte_predictions_FibroMarkers_with_GeneName.tsv"),
            delim="\t",col_names = T)


##
##> Check #GENES -> #Peptides predicted
##

outcyte_faps %>% subset(Prediction =="UPS") %>% nrow #> 95
outcyte_faps %>% subset(Prediction =="Signal-peptide") %>% nrow #> 326

outcyte_faps$GeneSymbol %>% unique %>% length()
## 203

outcyte_fibro %>% subset(Prediction =="UPS") %>% nrow #> 15
outcyte_fibro %>% subset(Prediction =="Signal-peptide") %>% nrow #>84

outcyte_fibro$GeneSymbol %>% unique %>% length
##  44
