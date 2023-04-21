#!/usr/bin/env R

#library(clusterProfiler)
# convertedIDs_df <- clusterProfiler::bitr(gene_list,
#                                   fromType = from,
#                                   toType = to,
#                                   OrgDb  = org.Hs.eg.db)


convert_gene_ids.hs <- function(id_list,from="SYMBOL",to="ENTREZID") {

  library(org.Hs.eg.db) # installed via biocLite. Cannot call library from an R-variable/object
  require(AnnotationDbi)
  converted_id <- mapIds(x=org.Hs.eg.db,
                         keys=id_list,
                         keytype=from, # ENSEMBL, ENTREZID, SYMBOL
                         column=to,
                         multiVals = "first")
  return(converted_id)

}

convert_gene_ids.mm <- function(id_list,from="ENSEMBL", to="SYMBOL" ) {

  library(org.Mm.eg.db) # installed via biocLite. Cannot call library from an R-variable/object
  require(AnnotationDbi)
  converted_id <- mapIds(x=org.Mm.eg.db,
                         keys=id_list,
                         keytype=from,  # ENSEMBL, ENTREZID, SYMBOL
                         column=to,
                         multiVals = "first")
  return(converted_id)

}


convertID_EnsemblDB_mm <- function(id_list){
    gene_symbols <- ensembldb::select(EnsDb.Mmusculus.v79,
                                  keys= id_list,
                                  keytype = "GENEID",
                                  columns = "SYMBOL")

    return(gene_symbols)
}
