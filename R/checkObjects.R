#################################
## checkInputs(X,time,geneset) ##
#################################

checkInputs <- function(X, tt, genesets){
  ## 'genesets' should be a list
  if (!is.list(genesets)){
    stop("'genesets' is not in required format.")
  }
  ## gene set names much be stored in the names of the list
  if (is.null(names(genesets))){
    stop("'genesets' must have names.")
  }
  
  ## tt should be numerical vector
  if (!is.numeric(tt)){
    stop("'tt' must numeric.")
  }
  
  ## rownames(X) and genesets must use the same gene annotation
  if (sum(rownames(X) %in% unlist(genesets))<1){
    stop("'X' is not in required format. Genes in rows and time points in columns. Row names should be the gene IDs using the same annotation as in 'genesets'.")
  }
  
  return(list(X=X, tt=tt, genesets=genesets))
}

