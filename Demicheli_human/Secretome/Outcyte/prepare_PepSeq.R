#!/usr/bin/env R

###> Author: Maria Kondili
###> Date: 26/04/2022
###> Subject : Read gene-Markers of Fibroblasts clusters (analyse_with_Seurat.Rmd) and extract peptide-sequences ,to be predicted
###> in next step by Outcyte for Secretome.

suppressPackageStartupMessages(library(tidyverse))

work_dir <- "Secretome_in_SingleCell_RNAseq_skeletal_muscle/"
secretome_dir <- paste0(work_dir,"DeMicheli_human/Secretome/")
#setwd(secretome_dir)
out_dir <- paste0(secretome_dir,"PepSeq_For_Outcyte/")



####
#### Search Peptide-sequences via BioMart  for each gene-set ###
####

source(paste0(work_dir,"functions/prepare_PepSeq_for_Outcyte_function.R")


#> 1/
fibro_col1a1_markers <- readr::read_delim("markers_of_COL1A1_fibroblasts_hs.tsv",col_names=T,delim="\t")

fibro_col1a1_PepSeqs <- prepare_seq_Outcyte(fibro_col1a1_markers$GeneName,
                                           org_dataset="hsapiens_gene_ensembl",
                                           gene_format="external_gene_name",
                                           out_dir= out_dir,
                                           suffix_outfile="PepSeq_COL1A1_Fibroblasts_Genes_human")


#> 2/
fibro_myoc_markers <- readr::read_delim("markers_of_MYOC_fibroblasts_hs.tsv",col_names=T,delim="\t")

fibro_myoc_PepSeqs <- prepare_seq_Outcyte(fibro_myoc_markers$GeneName,
                                           org_dataset="hsapiens_gene_ensembl",
                                           gene_format="external_gene_name",
                                           out_dir=out_dir ,
                                           suffix_outfile="PepSeq_MYOC_Fibro_Genes_human")

#>3 /
fibro_fbn1_markers <- readr::read_delim("markers_of_FBN1_fibroblasts_hs.tsv",col_names=T,delim="\t")

fibro_fbn1_Pepseqs <- prepare_seq_Outcyte(fibro_fbn1_markers$GeneName,
                                          org_dataset="hsapiens_gene_ensembl",
                                          gene_format="external_gene_name",
                                          out_dir= out_dir,
                                          suffix_outfile="PepSeq_FBN1_Fibro_Genes_human")
