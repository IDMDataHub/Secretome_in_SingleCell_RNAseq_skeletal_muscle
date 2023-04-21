#!/usr/bin/env R

###>Author: Maria Kondili
###>Date : April 2022
###>Subject : Read Outcyte results(DeMicheli-human Fibroblasts) of each marker-gene-set and add GeneName column

library(tidyverse)
library(biomaRt)
library(stringr)
library(httr)
library(seqinr)

## original_fibro_pred <- read_delim("Outcyte/Results/PepSeq_Fibroblasts_Markers_human_pipeline.txt",col_names=F,delim="\t")

convert_Entrez2Symbol.hs <- function(entrezid) {

  library(org.Hs.eg.db)
  library(AnnotationDbi)
  converted_ids <- mapIds(x=org.Hs.eg.db,
                          keys=entrezid, #to give
                          column="SYMBOL",  # to receive
                          keytype="ENTREZID",
                          multiVals="first")
  return(converted_ids)
}

annotate_outcyte_genes <- function(outcyte_dir,filename,outname) {

    suppressPackageStartupMessages(library(tidyverse))

    pred <- readr::read_delim(paste0(outcyte_dir,filename),col_names=F)

    colnames(pred) <- c("EntrezID","Prediction","Score")

    pred$GeneSymbol  <- convert_Entrez2Symbol.hs(as.character(pred$EntrezID))

    readr::write_delim(pred,paste0(outcyte_dir,outname), delim="\t",col_names = T)

    return(pred)
}


###> All Fibroblasts Predictions
work_dir    <- "/projects/single_cell_skeletal_muscle/DeMicheli_human/Secretome/"
outcyte_dir <- paste0(work_dir,"Outcyte/")
cat("\nOutcyte_dir=", outcyte_dir,"\n")


# all_fibro_pred <- annotate_outcyte_genes(outcyte_dir,"PepSeq_Fibroblasts_Markers_human_pipeline.txt",
#                   "OutcytePred_FibroblastsMarkers_DeMicheli_human_with_GeneName.tsv")
#
# all_fibro_pred %>% dim

###>
###> Fibro-Subgroups-Predictions
###>

col1a1_markers_pred <- annotate_outcyte_genes(outcyte_dir,"PepSeq_COL1A1_Fibro_Genes_human_pipeline.txt",
                                              "OutcytePred_COL1A1_FibroMarkers_DeMicheli_human_with_GeneName.tsv")

fbn1_markers_pred   <- annotate_outcyte_genes(outcyte_dir,"PepSeq_FBN1_Fibro_Genes_human_pipeline.txt",
                                              "OutcytePred_FBN1_Fibro_DeMicheli_human_with_GeneName.tsv")

myoc_markers_pred   <- annotate_outcyte_genes(outcyte_dir,
                                              "PepSeq_MYOC_Fibro_Genes_human_pipeline.txt",
                                              "OutcytePred_MYOC_Fibro_DeMicheli_human_with_GeneName.tsv")



col1a1_markers_pred %>% dim
#[1] 4009    4

fbn1_markers_pred %>% dim
#[1] 10079    4

myoc_markers_pred %>% dim
#[1] 6155    4



#### VENN of Secreted Genes :
library(ggvenn)

secreted_genes_col1a1 <- subset(col1a1_markers_pred, Prediction %in% c("Signal-peptide","UPS"))
secreted_genes_fbn1   <- subset(fbn1_markers_pred , Prediction %in% c("Signal-peptide","UPS"))
secreted_genes_myoc   <- subset(myoc_markers_pred ,Prediction %in% c("Signal-peptide", "UPS"))


ggvenn( list("Fibroblasts-1"=secreted_genes_col1a1$GeneSymbol,
             "Fibroblasts-2"=secreted_genes_fbn1$GeneSymbol ,
             "Fibroblasts-3"=secreted_genes_myoc$GeneSymbol),
        fill_color = c("#CD534CFF","#0073C2FF", "#EFC000FF"),
        text_color = "black",stroke_size = 0.5,show_percentage = TRUE,
        set_name_size = 4,text_size = 4 )
