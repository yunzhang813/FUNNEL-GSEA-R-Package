## Required packages: fda, quadrupen, pheatmap, RColorBrewer


#############################
## FUNCTION: FPCA.Fstats() ##
#############################

FPCA.Fstats <- function(X,tt,rr=rep(1,length(tt)),selection_k="FVE",FVE_threshold=0.9)
### INPUT ###
# X: n*m data matrix, with missing values denoted as NA.
# tt: length mm vector, unique time points.
# rr: number of repetitions at each unique time point.
# selection_k: the method of choosing the number of principal components;
#              "FVE" (fraction of variance explained) : use scree plot
#                           approach to select number of principal
#                           components), see "FVE_threshold" below;
#              positive integer K: user-specified number of principal components.
# FVE_threshold: a positive number between 0 and 1; It is used with the option
#                selection_k = "FVE" to select the number of principal components
#                that explain at least "FVE_threshold" of total variation.
### EXPORT ###
# stat = Functional F-statistics
{
  y = X
  res = PCA(y,tt,rr,selection_k=selection_k,FVE_threshold=FVE_threshold,verbose="off")
  n = nrow(y)
  m = ncol(y)
  mm = length(tt)
  rr1 = rep(1:mm,rr)
  mu = res$mu
  yfit = res$yfit_orig
  # if(delta=="auto"){ delta = res$sigma }
  delta = res$sigma
  ss0 = rowSums((y-matrix(rep(mu,m),ncol=m))^2,na.rm=T)
  ss1 = rowSums((y-yfit)^2,na.rm=T)
  stat = (ss0-ss1)/(ss1+delta)
  return(stat)
} 


##################################
## FUNCTION: equiv.regression() ##
##################################

equiv.regression <- function(yfd, xfd, threshold=0.01)
{	
  ## inner product matrices
  Sigma <- inprod(xfd,xfd)
  XY <- inprod(xfd,yfd)
  ## eigen decomposition for covariance matrix
  eig <- eigen(Sigma)
  cutoff <- sum(eig$values)*threshold
  l <- sum(eig$values > cutoff)
  ## equivalence
  X <- diag(sqrt(eig$values[1:l])) %*% t(eig$vectors[,1:l])
  Y <- diag(sqrt(1/eig$values[1:l])) %*% t(eig$vectors[,1:l]) %*% XY
  ## output
  return(list("y"=Y, "Xmat"=X, "cut.value"=cutoff))
}

####################################################
## FUNCTION: wMWUTest() ##
####################################################

wMWUTest <- function(test.index,statistics,weight=NULL,correlation=0,df=Inf)
# Weighted rank sum test for correlated samples.
# See also limma::rankSumTestWithCorrelation.
{
  n <- length(statistics)
  n1 <- length(test.index)
  n2 <- n-n1
  ## Mann-Whitney style ranks
  r <- vector()
  statistics0 <- setdiff(statistics, statistics[test.index])
  for (i in 1:length(test.index)){r[i] <- sum(statistics[test.index[i]] > statistics0)}
  ## weight
  if(is.null(weight)){weight <- rep(1, n1)}
  U.w <- sum(weight*r)
  n1.w <- sum(weight)
  c1 <- sum(weight^2)
  c2 <- n1.w^2 - c1
  ## mean and variance estimates of Mann-Whitney U-statistic
  mu <- n1.w*n2/2
  if(correlation==0 || n1==1) {
    sigma2 <- c1*n2*(n2+2)/12 + c2*n2/12
  } else {
    sigma2 <- asin(1)*c1*n2 + asin(0.5)*c1*n2*(n2-1) + asin(correlation/2)*c2*n2*(n2-1) + asin((correlation+1)/2)*c2*n2
    sigma2 <- sigma2/2/pi
  }
  ## potential improvment: easy to have ties, but minor effect
  TIES <- (length(r) != length(unique(r))) 
  if(TIES) {
    NTIES <- table(r)
    adjustment <- sum(NTIES*(NTIES+1)*(NTIES-1)) / (n*(n+1)*(n-1)) 
    sigma2 <- sigma2 * (1 - adjustment)
  }
  ## z-score and p-value
  zlowertail <- (U.w+0.5-mu)/sqrt(sigma2)
  zuppertail <- (U.w-0.5-mu)/sqrt(sigma2)
  pvalues <- c(less=pt(zlowertail,df=df), greater=pt(zuppertail,df=df,lower.tail=FALSE))
  pvalues	
}

############################
## FUNCTION: FUNNELtest() ## getWeightMatrix(), wMWUTest()
############################

FUNNELtest <- function(fdexpr, genesets, Fstats, rho, df, nharm=3, centerfns=FALSE, equiv.threshold=0.01, lam1=0.4, lam2=0.01)
## IMPORT ##
## fdexpr = fdobj for the expression matrix
## genesets = list of gene sets
## genenames = names of all genes in fdexpr (same annotation as in genesets)
## Fstats = vector of functional F-statistics per gene
## rho = mean of per pathway correlation
## df = degrees of freedom for wMWUTest 
## OUTPUT ##
## pvals = p-values per gene set
## weight.list = list of vectors of weights per gene set
{
  pvals <- vector()
  weight.list <- list()
  cat("Weight calculation...", "\n")
  weight.mat <- getWeightMatrix(fdexpr, genesets, nharm=nharm, centerfns=centerfns, equiv.threshold=equiv.threshold, lam1=lam1, lam2=lam2)
  cat("Gene set test...", "\n")
  for (k in 1:length(genesets)){
    testset <- genesets[k] # list of length 1
    weight.k <- weight.mat[,names(testset)]
    weight <- weight.k[!is.na(weight.k)]
    weight <- weight[unlist(testset)]
    test <- wMWUTest(test.index=unlist(testset), statistics=Fstats, weight=weight, correlation=rho, df=df)
    pvals[k] <- test["greater"]
    weight.list[[k]] <- weight
  }
  names(pvals) <- names(weight.list) <- names(genesets)
  ## output
  return(list("pvals"=pvals, "weight.list"=weight.list))
}


#############################
## FUNCTION: FUNNEL.GSEA() ##
#############################

FUNNEL.GSEA <- function(X, tt, genesets, lambda=10^-3.5, rr=rep(1,length(tt)), selection_k="FVE", FVE_threshold=0.9, nharm=3, centerfns=FALSE, equiv.threshold=0.01, lam1=0.4, lam2=0.01, alpha.level=0.05)
### INPUT ###
# X = original expression matrix
# tt = origninal tt points
# genesets = original gene sets
### OUTPUT ###
# pvals = p-values per gene set
# weight.list = list of vectors of weights per gene set
# correlation = mean pathway correlation
# sig.genesets = significant gene sets
# Fstats = functional F-statistics per gene
{
  checkInputs(X, tt, genesets)
  ## Standardize timepoint and X so that the optimum roughness/L1/L2
  ## penality parameters are applicable
  tt <- (tt - min(tt))/diff(range(tt))
  X <- t(scale(t(X)))
  ## Remove genes in predefined genesetss that are not present in X the
  ## filtered input data
  genesets <- lapply(genesets, function(z) { intersect(z, rownames(X)) })
  ## get Fstats
  Fstats <- FPCA.Fstats(X, tt, rr=rr, selection_k=selection_k, FVE_threshold=FVE_threshold)
  ## get rho.hat
  rho.hat <- getRho(X, genesets)
  ## smoothing
  fdexpr <- smoothExpr(X, tt, lambda=lambda)
  ## FUNNEL test
  FUNNEL.out <- FUNNELtest(fdexpr, genesets, Fstats=Fstats, rho=rho.hat, df=sum(rr)-1, nharm=nharm, centerfns=FALSE, equiv.threshold=equiv.threshold, lam1=lam1, lam2=lam2)
  ## significant pathways
  sig.genesets <- names(genesets)[FUNNEL.out$pvals<alpha.level]
  ## output
  return(list("pvals"=FUNNEL.out$pvals, "weight.list"=FUNNEL.out$weight.list, "correlation"=rho.hat, "sig.genesets"=sig.genesets, "Fstats"=Fstats))
}


###############################
## FUNCTION: weightPerGene() ##
###############################

weightPerGene <- function(weight.list, genesOfInterest)
### IMPORT ###
# weight.list = output from FUNNELtest() or FUNNEL.wrapper()
# genesOfInterest = vector of genes of interest
### OUTPUT ###
# weightPerGene.list = list of weights associated with each gene
{
  genesets <- lapply(weight.list,names)
  weightPerGene.list <- vector("list", length=length(genesOfInterest))
  names(weightPerGene.list) <- genesOfInterest
  for(gene.i in genesOfInterest){
    ## find pathways that contain gene.i
    index <- which(sapply(genesets, function(z){gene.i %in% z}))
    ## three cases
    if(length(index)==0){weight <- NA}
    if(length(index)==1){weight <- 1}
    if(length(index)>1){
      ## find location of gene.i in each of its pathways
      weight <- vector()
      for (kk in 1:length(index)){
        index.k <- index[kk]
        location <- which(genesets[[index.k]]==gene.i)
        weight[kk] <- weight.list[[index.k]][location]
      }
      names(weight) <- names(genesets[index])
    }
    weightPerGene.list[[gene.i]] <- weight
  }
  return(weightPerGene.list)
}


############################
## FUNCTION: plotWeight() ## 
############################

plotWeight <- function(weight.list, geneset.index, main="Weight", show_colnames=TRUE, ...)
## Plot of non-zero weights per gene set
{
  weight <- sort(weight.list[[geneset.index]], decreasing=TRUE)
  weight <- weight[weight!=0]
  weight <- as.matrix(weight)
  if (is.numeric(geneset.index)){colnames(weight) <- names(weight.list)[geneset.index]}
  if (is.character(geneset.index)){colnames(weight) <- geneset.index}
  
  genes <- rownames(weight)
  genesets <- lapply(weight.list,names)
  membership <- table(unlist(genesets))
  anno_row <- as.data.frame(membership[genes])
  names(anno_row) <- "Membership"
  anno_row[anno_row>1,] <- "Overlapping"
  anno_row[anno_row==1,] <- "Exclusive"
  
  pheatmap(weight, scale="none",
           cellwidth = 60, cellheight = 20,
           color=colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")))(10), breaks=seq(0,1,.1),
           cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE, show_colnames=show_colnames,
           annotation_row=anno_row, display_numbers=TRUE, main=main, ...)
}




