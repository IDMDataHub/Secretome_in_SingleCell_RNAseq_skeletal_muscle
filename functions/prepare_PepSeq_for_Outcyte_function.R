
prepare_seq_Outcyte <- function( my_input_genes,
                                org_dataset="mmusculus_gene_ensembl",
                                gene_format="external_gene_name",out_dir, suffix_outfile){

    library(tidyverse)
    library(biomaRt)
    library(dplyr)
    library(stringr)
    library(httr)
    library(seqinr)

    #listMarts()
    ensembl = useMart("ensembl")
    ensembl = useDataset(org_dataset,mart=ensembl)

    filters = listFilters(ensembl)
    attributes_bm <- c("peptide",
                       "entrezgene_id",
                       "external_gene_name",
                       "uniprotswissprot",
                       "uniprotsptrembl",
                       "ensembl_gene_id",
                       "ensembl_transcript_id",
                       "ensembl_peptide_id")

    ### Transform gene-name --> EntrezID
    EntrezID <- getBM(attributes = c("entrezgene_id" ,gene_format), #what I want
                      filters = gene_format, # or "hgnc_symbol",  # what I give
                      values  = my_input_genes,
                      mart    = ensembl )

    ##> ? genes not converted:
    idx_not_convrt <-  which(! my_input_genes %in% EntrezID[,gene_format] )

    if (length(idx_not_convrt) > 0) {
      my_input_genes <- my_input_genes[-idx_not_convrt]
    }

    # remove duplicate gene
    if ( any(duplicated(EntrezID[,gene_format])))  {
      EntrezID <- EntrezID[-which(duplicated(EntrezID[,gene_format])),]
    }


    # listFilters(ensembl.hs )$name --> to find filter of gene-name-format you give ,e.g : external_gene_name, ensembl_gene_id

    timeout(100000) ##for giving time to getBM

    protein_bm <- getBM(attributes = attributes_bm, # what you want to get
                        filters = "entrezgene_id",  # db of ID  you give
                        values  = EntrezID$entrezgene_id, # ids to convert
                        mart    = ensembl  )


    ##Nettoyage de protein_bm
    protein_bm_curated <- subset(protein_bm, peptide != "Sequence unavailable")

    any(protein_bm_curated$peptide == "Sequence unavailable")


    protein_bm_curated$peptide <- as.vector(as.character(protein_bm_curated$peptide))
    protein_bm_curated$peptide <- str_remove(protein_bm_curated$peptide, pattern="\\*$")
    protein_bm_curated <- subset( protein_bm_curated, nchar(peptide)> 20 )

    ##Export FASTA  with ENTREZ ID & Gene Symbol
    ##Prepare dataframe for generation of multifasta : fuse information columns

    PepSeqs <- data.frame("Peptide"= protein_bm_curated$peptide,
                          "Info"= paste(protein_bm_curated$entrezgene_id,
                                               "|", protein_bm_curated$external_gene_name,
                                               "|", protein_bm_curated$uniprotswissprot,
                                               "|",protein_bm_curated$uniprotsptrembl,
                                               "|",protein_bm_curated$ensembl_gene_id,
                                               "|", protein_bm_curated$ensembl_transcript_id,
                                               "|",protein_bm_curated$ensembl_peptide_id,
                                               "|","stop"))

    exportFASTA(PepSeqs, paste0(out_dir , suffix_outfile,".fasta"))

    PepSeq_df <- data.frame("Peptide"= protein_bm_curated$peptide,
                            "EntrezID"= protein_bm_curated$entrezgene_id,
                            "GeneName"=protein_bm_curated$external_gene_name,
                            "TranscriptID"=protein_bm_curated$ensembl_transcript_id,
                            "PeptideID"= protein_bm_curated$ensembl_peptide_id)
    return(PepSeq_df)
}


##> if we want to define number of elements of array=number of fasta created :
##> length.out = total number of fasta created
