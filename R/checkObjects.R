#################################
## checkInputs(X,time,geneset) ##
#################################

checkInputs <- function(X,time,geneset){
  ## check ncol(X)==length(time)
  if (ncol(X)!=length(time)){
    stop("Number of columns in 'X' should be the same as length of 'time'.")
  }
  
  ## check 'geneset' should be a list
  if (!is.list(geneset)){
    stop("The 'geneset' should be a list.")
  }
  
  ## check rownames(X) and geneset use the same gene annotation
  if (sum(rownames(X) %in% unlist(geneset))<100){
    stop("Row names of 'X' should be denoted as gene names, and they should be using the same annotation as the genes in 'geneset'.")
  }
  
  return(list(X=X,time=time,geneset=geneset))
}
