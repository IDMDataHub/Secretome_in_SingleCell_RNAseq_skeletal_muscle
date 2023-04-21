#!/usr/bin/env R

###>Author: Maria Kondili
###>Date : May 2022
###>Subject: Extract peptide-seq from Biomart, of FAPs-markers genes
###>         Launch Outcyte with the resulting .fasta files.


suppressPackageStartupMessages(library(tidyverse))

## Call function that converts geneNames->EntrezID & extracts Peptide-seq of each gene in Fasta from BiomaRt
source("../functions/prepare_PepSeq_for_Outcyte_function.R")

proj_dir  <- "/projects/single_cell_skeletal_muscle/Oprescu_mouse/"
work_dir  <- paste0(proj_dir,"Secretome/")
out_dir   <- paste0(work_dir,"PepSeq_For_Outcyte/")
dir.create(out_dir,showWarnings = F)

#>> FAPs-Genes only:

faps_markers_uniq <- read_delim(paste0(proj_dir, "Seurat_Analysis/FAPs_Markers_6subgroups_Oprescu.txt"),delim="\n")

PepSeq_faps <- prepare_seq_Outcyte(faps_markers_uniq$GeneName,
                                      org_dataset = "mmusculus_gene_ensembl",
                                      gene_format = "external_gene_name",
                                      out_dir     = out_dir,
                                      suffix_outfile = "PepSeq_FAPsMarkers_Oprescu" )

## PepSeq_faps %>% dim
#> 647

### Fibroblasts

fibro_markers <- read_delim( paste0(proj_dir,"Table_S2_FAP_gene_signature_list_per_subtype.tsv"),
                                delim="\t",col_names = T)

fibro_markers <- fibro_markers$Fibroblasts
#>50

PepSeq_fibro_mm <- prepare_seq_Outcyte(fibro_markers,
                                      org_dataset="mmusculus_gene_ensembl",
                                      gene_format = "external_gene_name",
                                      out_dir = out_dir ,
                                      suffix_outfile = "PepSeq_FibroMarkers_Oprescu")
PepSeq_fibro_mm %>% nrow
##> 154
