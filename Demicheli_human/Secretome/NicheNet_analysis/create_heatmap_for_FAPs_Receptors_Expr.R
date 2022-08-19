#!/bin/R


suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))


###---
###--- LOCATE YOURSELF :
###---

work_dir <- "/projects/single_cell_skeletal_muscle/DeMicheli_human/Secretome/NicheNet_analysis/"
input_dir <- paste0(work_dir, "Input/")
output_dir <- paste0(work_dir,"FN1_receptors_Graphs/")
dir.create(output_dir,showWarnings = F)

###---
###--- GIVE INPUT VARIABLES
###---

gene_list <- paste0(work_dir,"FN1_receptors.tsv")

cell_types <- c( "FBN1+ MFAP5+ CD55+ Fibroblasts",
                     "DCN+ GSN+ MYOC+ Fibroblasts",
                     "COL1A1+ Fibroblasts",
                     "ACTA1+ Mature skeletal muscle",
                     "CLDN5+ PECAM1+ Endothelial",
                     "ICAM1+ SELE+ VCAM1+ Endothelial",
                     "PAX7+ DLK1+ MuSCs and progenitors",
                     "PAX7low MYF5+ MuSCs and progenitors",
                     "S100A9+ LYZ+ Inflammatory macrophages",
                     "C1QA+ CD74+ Macrophages")

## Table of merged expr is very large ,so not given here ! (created with another script for all genes )
merged_expr <-  readr::read_delim(paste0(work_dir,"Merged_MeanExpr+sem_All_CellTypes_SC_DeMicheli_all_Genes_human.tsv"),
                           delim="\t")

merged_expr %>% dim

#Turn to 0 the negative values :
merged_expr$MeanExpr <- ifelse(merged_expr$MeanExpr< 0 ,merged_expr$MeanExpr*(-1), merged_expr$MeanExpr  )

# remove 1st column=rownames
merged_expr <- merged_expr %>% dplyr::select(- c("...1"))


### READ INPUT
gene_table <- read_delim(gene_list, delim="\n",col_names=T)
gene_name  <- gene_table$hsap_gene_name

###
### Find which patients are not measured for all cell-types
###

pt_list <- list()
for (ct in cell_types) {
 pt <- subset(merged_expr,CellType== ct ) %>% pluck("PatientID") %>% as.factor() %>%  levels
 pt_list[[ct]] <- pt
}

#> Which is only in few cell-types :
all_patients <- merged_expr$PatientID %>% as.factor %>% levels()

#all_patients
#[1] "061419" "062419" "062719" "080619" "092618" "092712" "101018" "111418" "120518" "B"

print( map( seq(1,10), ~(which(! all_patients %in% pt_list[[.x]] ) ) %>% unlist) )

# [1]  2 10
# [1]  2 10
# [1]  2 10
# [1] 2 6
# [1] 1 9
# [1]  2 10
# [1] 2 6
# [1] 2 6
# [1]  2 10
# [1]  2 10

all_patients[c(2,6,10)]
# [1] "062419" "092712" "B"

##>> patient[2], patient[6] patient[10] shouldn't be included in histograms
lines2rmv <- which(grepl( paste(all_patients[c(2,6,10)], collapse="|"), merged_expr$PatientID))

merged_expr <- merged_expr[-lines2rmv, ]

patients_used <- all_patients[-c(2,6,10)]


####
#### MERGE Cell-Clusters of same Type & Mean over all Patients to have one value per Cell-type:
####


#> MEAN per Gene for all patients of each cell-type group:

reduced_merged_with_mean_expr <- function(gene_name , celltypes_array, new_ct_name, merged_expr){

  ct_lines <- which(merged_expr$CellType %in% celltypes_array)

  merged_mean_ct <- map(gene_name , ~(merged_expr[ct_lines,] %>%
                                          subset(GeneName == .x) %>%
                                          pluck("MeanExpr") %>% mean(na.rm=TRUE))) %>% unlist()

  merged_SEM_ct <- map(gene_name , ~(merged_expr[ct_lines,] %>%
                                        subset(GeneName == .x) %>%
                                        pluck("SEM") %>% mean(na.rm=TRUE))) %>% unlist()


  merged_mean_ct_df <- data.frame("GeneName"=gene_name ,
                                  "MeanExpr"=merged_mean_ct,
                                  "SEM" = merged_SEM_ct,
                                  "CellType"=new_ct_name )

  ##data.frame("GeneName"=gene_receptors, "MeanExpr"=merged_mean_ct,  "CellType"="fibroblasts")
  return(merged_mean_ct_df)
}


mrg_cell_types <- c("Fibroblasts", "Endothelials","Mature_sk_muscle", "MuSC_Progenitors","Macrophages" )


# 1/ "FBN1+ MFAP5+ CD55+ Fibroblasts","DCN+ GSN+ MYOC+ Fibroblasts","COL1A1+ Fibroblasts"
#fibro_lines <- which(merged_expr$CellType %in% c("FBN1+ MFAP5+ CD55+ Fibroblasts","DCN+ GSN+ MYOC+ Fibroblasts","COL1A1+ Fibroblasts"))

merged_mean_fibro_expr <- reduced_merged_with_mean_expr(gene_name ,
                          celltypes_array=c("FBN1+ MFAP5+ CD55+ Fibroblasts","DCN+ GSN+ MYOC+ Fibroblasts","COL1A1+ Fibroblasts"),
                          new_ct_name=mrg_cell_types[1], merged_expr)

# 2/ "CLDN5+ PECAM1+ Endothelial","ICAM1+ SELE+ VCAM1+ Endothelial"
# endoth_lines <- which(merged_expr$CellType %in% c("CLDN5+ PECAM1+ Endothelial","ICAM1+ SELE+ VCAM1+ Endothelial"))
merged_mean_endoth_expr <- reduced_merged_with_mean_expr(gene_name ,
                           celltypes_array=c("CLDN5+ PECAM1+ Endothelial","ICAM1+ SELE+ VCAM1+ Endothelial"),
                           new_ct_name=mrg_cell_types[2], merged_expr)

#3/ "ACTA1+ Mature skeletal muscle" (no merging for ct, it's only one cluster)
merged_mean_matureMu_expr <- reduced_merged_with_mean_expr(gene_name ,
                           celltypes_array=c("ACTA1+ Mature skeletal muscle"),
                           new_ct_name=mrg_cell_types[3], merged_expr)


# 4/  "PAX7+ DLK1+ MuSCs and progenitors", "PAX7low MYF5+ MuSCs and progenitors"
# musc_lines <- which(merged_expr$CellType %in% c("PAX7+ DLK1+ MuSCs and progenitors", "PAX7low MYF5+ MuSCs and progenitors"))
merged_mean_MuSC_expr <- reduced_merged_with_mean_expr(gene_name ,
                         celltypes_array = c("PAX7+ DLK1+ MuSCs and progenitors", "PAX7low MYF5+ MuSCs and progenitors"),
                         new_ct_name=mrg_cell_types[4], merged_expr)


# 5/ "S100A9+ LYZ+ Inflammatory macrophages", "C1QA+ CD74+ Macrophages"
# macroph_lines <- which(merged_expr$CellType %in% c("S100A9+ LYZ+ Inflammatory macrophages", "C1QA+ CD74+ Macrophages"))

merged_mean_macroph_expr <- reduced_merged_with_mean_expr(gene_name ,
                            celltypes_array = c("S100A9+ LYZ+ Inflammatory macrophages","C1QA+ CD74+ Macrophages"),
                            new_ct_name=mrg_cell_types[5], merged_expr)


###  Concatenate all smaller merged-data.frames in One !
merged_mean_expr_gene_receptors <- rbind(merged_mean_fibro_expr,
                                         merged_mean_endoth_expr,
                                         merged_mean_matureMu_expr,
                                         merged_mean_MuSC_expr,
                                         merged_mean_macroph_expr )

###
### PLOT
###

##> with Facet-Grid :  All cell-types at once in multi-window graph

png(paste0(output_dir,"Histogram_Expr_per_CellType_facet_grid_DeMicheli_hsap.png"))
  p <- ggplot(merged_mean_expr_gene_receptors,
         aes(x=GeneName, y=MeanExpr,fill=GeneName)) +
        geom_bar(stat="identity", width=0.9, position="dodge") +
        geom_errorbar(aes(ymin=MeanExpr-SEM, ymax=MeanExpr+SEM, color=GeneName),
                      position="dodge",size=0.3) +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 90, size=6.5)) +
        labs(x="FN1 Receptors Genes", y="mean expression(from 7 patients)") +
        facet_grid(~CellType)
  print(p)
dev.off()


### HEATMAP
#> Transform the merged_mean_expr table in a zscore-column-wise counts
merged_mean_expr_gene_receptors$CellType <- factor(merged_mean_expr_gene_receptors$CellType,
                                                   levels=c("Fibroblasts","Endothelials",
                                                            "Mature_sk_muscle","MuSC_Progenitors",
                                                            "Macrophages") )

merged_mean_expr_receptors_spread <- merged_mean_expr_gene_receptors %>%
                                      dplyr::select(-"SEM") %>%
                                      spread(key = "CellType",value = "MeanExpr")


rownames(merged_mean_expr_receptors_spread)<- merged_mean_expr_receptors_spread$GeneName
merged_mean_expr_receptors_spread <- merged_mean_expr_receptors_spread %>% dplyr::select(-"GeneName")

colnames(merged_mean_expr_receptors_spread) <- c("Fibroblasts","Endothelials","Mature.sk.mu",
                                                 "MuSC/Progenitors","Macrophages")


write_delim(merged_mean_expr_receptors_spread,
            paste0(output_dir,"Merged_Mean_Expr_FN1_receptors_for_all_celltypes.txt"),
            delim="\t",col_names=T)

library(RColorBrewer) # display.brewer.all() ->show names and colours
library(pheatmap)
png("Heatmap_FN1_receptors_for_all_cell_types.png")
  p <- pheatmap(merged_mean_expr_receptors_spread,
               color = colorRampPalette(brewer.pal(n=9,name ="PuRd"))(50),
               clustering_distance_rows="euclidean",
               cluster_cols=FALSE,
               annotation_legend = FALSE,
               annotation_names_col = TRUE,
               show_rownames=TRUE,
               labels_row = rownames(merged_mean_expr_receptors_spread),
               angle_col = 45)
  print(p)
dev.off()
