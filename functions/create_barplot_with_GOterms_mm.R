create_barplot_GOterms.mm <- function(genenames,gene_format = "SYMBOL",my_ontol="ALL",main="",in_png=TRUE){

  suppressPackageStartupMessages(library(clusterProfiler))
  suppressPackageStartupMessages(library(org.Mm.eg.db))
  suppressPackageStartupMessages(library(AnnotationDbi))
  suppressPackageStartupMessages(library(DOSE))
  suppressPackageStartupMessages(library(enrichplot))

  if (gene_format !="ENTREZID"){
    EntrezID <- mapIds(x=org.Mm.eg.db, keys=genenames,
                      column="ENTREZID",
                      keytype=gene_format,
                      multiVals = "first")
      ##to see accepted keytypes :  keytypes(org.Mm.eg.db)
      ## column = what I want to obtain
      ## keytype = what I give in keys
    } else {
      EntrezID <- genenames
    }

    go.enrichm <- enrichGO(EntrezID, ont=my_ontol,
                            OrgDb=org.Mm.eg.db,
                            keyType = "ENTREZID",
                            pvalueCutoff=0.1,
                            qvalueCutoff=0.1,
                            pAdjustMethod="BH",
                            readable=TRUE)

    #sort GeneRatio from higher to lower :
    go.enrichm@result$GeneRatio <- go.enrichm@result$GeneRatio[order(go@result$Count,decreasing = T)]


  if (my_ontol=="BP") {msg = "Biol_Processes" }
  if (my_ontol=="MF") {msg = "Molec_Functions"}
  if (my_ontol=="CC") {msg = "Cellular_Components" }
  if (my_ontol== "ALL") {msg="All_GO_terms"}

  dir.create("GO_terms_Plots/",showWarnings = F)

  p <- barplot(go.enrichm, showCategory = 15,
              x="GeneRatio", color="p.adjust", font.size=8,
              title=paste(msg," of ", main))

  if ( in_png==TRUE) {
      png(paste0("GO_terms_Plots/Barplot_",msg,"_",main,".png"))
          print(p)
      dev.off()
  } else {
      print(p)
    }

}
