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


convert_Entrez2Symbol.mm <- function(EntrezID) {

    library(org.Mm.eg.db) # for Human IDs
    library(AnnotationDbi)
    converted_ids <- mapIds(x=org.Mm.eg.db,
                            keys=EntrezID, #to give
                            column="SYMBOL",  # to receive
                            keytype="ENTREZID",
                            multiVals = "first")
    return(converted_ids)
}

work_dir <- "/projects/single_cell_skeletal_muscle/Oprescu_mouse/Secretome/"

###> All FAPs:

outcyte_dir <- paste0(work_dir,"Outcyte_Results/")

outcyte_faps <- readr::read_delim(paste0(outcyte_dir,"PepSeq_FAPsMarkers_Oprescu_pipeline.txt"),col_names=F)
colnames(outcyte_faps) <- c("EntrezID","Prediction","Score")

outcyte_faps$GeneSymbol  <- convert_Entrez2Symbol.mm(as.character(outcyte_faps$EntrezID))

write_delim(outcyte_faps,
            paste0(outcyte_dir,"Outcyte_predictions_FAPsMarkers_with_GeneName.tsv"),
            delim="\t",col_names = T)

###> Fibroblasts :

outcyte_fibro <- readr::read_delim(paste0(outcyte_dir,"PepSeq_FibroMarkers_Oprescu_pipeline.txt"),col_names=F)
colnames(outcyte_fibro) <- c("EntrezID","Prediction","Score")

outcyte_fibro$GeneSymbol  <- convert_Entrez2Symbol.mm(as.character(outcyte_fibro$EntrezID))

write_delim(outcyte_fibro,
            paste0(outcyte_dir,"Outcyte_predictions_FibroMarkers_with_GeneName.tsv"),
            delim="\t",col_names = T)


##
##> Check #GENES -> #Peptides predicted
##

outcyte_faps %>% subset(Prediction =="UPS") %>% nrow #> 95
outcyte_faps %>% subset(Prediction =="Signal-peptide") %>% nrow #> 326
## total = 421
outcyte_faps$GeneSymbol %>% unique %>% length()
##unique = 203

outcyte_fibro %>% subset(Prediction =="UPS") %>% nrow #> 15
outcyte_fibro %>% subset(Prediction =="Signal-peptide") %>% nrow #>84
## total = 99
outcyte_fibro$GeneSymbol %>% unique %>% length
## unique = 44
