ensemblToGeneId <- function(ensemblId,dataSet='hsapiens_gene_ensembl'){
  ensemblToGeneIdInternal <- function(x,ensembl){
    gene1 <- NA;
    try(gene1 <- getBM(attributes='external_gene_name',filters='ensembl_gene_id',values=x,mart=ensembl),silent=T)
    if(is.null(gene1)){
      gene1<-NA
    }
    if(length(unlist(gene1))==0){
      gene1=NA
    }
    return(gene1)
  }
  require(dplyr)
  require(biomaRt)
  ensembl=useMart("ensembl")
  ensembl = useDataset(dataSet,mart=ensembl)
  geneId <- ensemblId %>%
    sapply(ensemblToGeneIdInternal,ensembl) %>%
    unlist() %>%
    as.character()
  return(geneId)
}