#################################
## checkInputs(X,time,geneset) ##
#################################

checkInputs <- function(X, genesets){
  ## 'genesets' should be a list
  if (!is.list(genesets)){
    stop("'genesets' is not in required format.")
  }
  
  ## rownames(X) and genesets must use the same gene annotation
  if (sum(rownames(X) %in% unlist(genesets))<1){
    stop("'X' is not in required format. Genes in rows and time points in columns. Row names should be the gene IDs using the same annotation as in 'genesets'.")
  }
  
  return(list(X=X, genesets=genesets))
}
