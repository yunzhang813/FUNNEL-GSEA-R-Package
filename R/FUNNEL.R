## Required packages: fda, quadrupen, pheatmap, RColorBrewer

###########################
## FUNCTION: scaleTime() ##
###########################

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

wMWUTest <- function(index,statistics,weight=NULL,correlation=0,df=Inf)
  # weighted Rank sum test as for two-sample Wilcoxon-Mann-Whitney test,
  # and allowing for correlation between members of test set.
  # Edited from limma::rankSumTestWithCorrelation by Gordon Smyth and Di Wu.
{
  n <- length(statistics)
  n1 <- length(index)
  n2 <- n-n1
  ## Mann-Whitney style ranks
  r <- vector()
  statistics0 <- setdiff(statistics, statistics[index])
  for (i in 1:length(index)){r[i] <- sum(statistics[index[i]] > statistics0)}
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
## FUNCTION: FUNNELtest() ## getWeight(), wMWUTest()
############################

FUNNELtest <- function(fdexpr, geneset, Fstats, rho, df, nharm=3, centerfns=FALSE, equiv.threshold=0.01, lam1=0.4, lam2=0.1)
  ## IMPORT ##
  ## fdexpr = fdobj for the expression matrix
  ## geneset = list of original pathway database
  ## genenames = name for all the genes in fdexpr, expressed with the same annotation as the geneset
  ## Fstats = vector of per gene F-statistic at gene level analysis
  ## rho = mean of per pathway correlation
  ## df = degree of freedom used in calculating per gene F-statistic 
  ## OUTPUT ##
  ## pvals = p-values for each pathway
  ## weight.list = list of weighting vectors for each pathway
{
  pvals <- vector()
  weight.list <- list()
  cat("Weight calculation...", "\n")
  weight.mat <- getWeightMatrix(fdexpr, geneset, nharm=nharm, centerfns=centerfns, equiv.threshold=equiv.threshold, lam1=lam1, lam2=lam2)
  cat("Gene set test...", "\n")
  for (k in 1:length(geneset)){
    testset <- geneset[k] # list of length 1
    weight.k <- weight.mat[,names(testset)]
    weight <- weight.k[!is.na(weight.k)]
    weight <- weight[unlist(testset)]
    test <- wMWUTest(index=unlist(testset), statistics=Fstats, weight=weight, correlation=rho, df=df)
    pvals[k] <- test["greater"]
    weight.list[[k]] <- weight
  }
  names(pvals) <- names(weight.list) <- names(geneset)
  ## output
  return(list("pvals"=pvals, "weight.list"=weight.list))
}

## example
# key.geneset <- c(122, 123, 124, 126)
# key.geneset <- true.pathways
# t0 <- system.time(temp0 <- FUNNELtest0(fdexpr, geneset, lam1=0.4, lam2=0.01, Fstats=Fstats, rho=rho.hat, df=15))
# temp0$pvals[key.geneset]
# 
# t <- system.time(temp <- FUNNELtest(fdexpr, geneset, lam1=0.4, lam2=0.01, Fstats=Fstats, rho=rho.hat, df=15))
# temp$pvals[key.geneset]
# 
# identical(temp$weight.list, temp0$weight.list) #TRUE


################################
## FUNCTION: FUNNEL.wrapper() ##
################################

FUNNEL.GSEA <- function(X, tt, geneset, lambda=10^-3.5, rr=rep(1,length(tt)), selection_k="FVE", FVE_threshold=0.9, nharm=3, equiv.threshold=0.01, lam1=0.4, lam2=0.1, alpha.level=0.05)
  ### INPUT ###
  # X = original expression matrix
  # tt = origninal tt points
  # geneset = original geneset database
{
  checkInputs(X,tt,geneset)
  ## Standardize timepoint and X so that the optimum roughness/L1/L2
  ## penality parameters are applicable
  tt <- (tt - min(tt))/diff(range(tt))
  X <- t(scale(t(X)))
  ## Remove genes in predefined genesets that are not present in X the
  ## filtered input data
  geneset <- lapply(geneset, function(z) { intersect(z, rownames(X)) })
  ## smoothing curve
  fdexpr <- smoothExpr(X, tt, lambda=lambda)
  ## get Fstats and rho
  Fstats <- getFstats(X, tt, rr=rr, selection_k=selection_k, FVE_threshold=FVE_threshold)
  rho.hat <- getRho(X, geneset)
  ## FUNNEL test
  FUNNEL.out <- FUNNELtest(fdexpr, geneset, Fstats=Fstats, rho=rho.hat, df=sum(rr)-1, nharm=nharm, centerfns=FALSE, equiv.threshold=equiv.threshold, lam1=lam1, lam2=lam2)
  ## significant pathways
  sig.geneset <- names(geneset)[FUNNEL.out$pvals<alpha.level]
  ## output
  return(list("pvals"=FUNNEL.out$pvals, "weight.list"=FUNNEL.out$weight.list, "correlation"=rho.hat, "sig.geneset"=sig.geneset))
}



#################################################################################################

###############################
## FUNCTION: weightPerGene() ##
###############################

weightPerGene <- function(weight.list, genesOfInterest)
  ### IMPORT ###
  # weight.list = output from FUNNELtest() or FUNNEL.wrapper()
  # genesOfInterest = vector of genes of interest
  ### OUTPUT ###
  # weightPerGene.list = list of weight for each gene in testgenes
{
  geneset <- lapply(weight.list,names)
  weightPerGene.list <- vector("list", length=length(genesOfInterest))
  names(weightPerGene.list) <- genesOfInterest
  for(gene.i in genesOfInterest){
    ## find pathways that contain gene.i
    index <- which(sapply(geneset, function(z){gene.i %in% z}))
    ## three cases
    if(length(index)==0){weight <- NA}
    if(length(index)==1){weight <- 1}
    if(length(index)>1){
      ## find location of gene.i in each of its pathways
      weight <- vector()
      for (kk in 1:length(index)){
        index.k <- index[kk]
        location <- which(geneset[[index.k]]==gene.i)
        weight[kk] <- weight.list[[index.k]][location]
      }
      names(weight) <- names(geneset[index])
    }
    weightPerGene.list[[gene.i]] <- weight
  }
  return(weightPerGene.list)
}


############################
## FUNCTION: plotWeight() ## plot of non-zero weights per geneset
############################

plotWeight <- function(weight.list, geneset.index, main="Weight", ...){
  weight <- sort(weight.list[[geneset.index]], decreasing=TRUE)
  weight <- weight[weight!=0]
  weight <- as.matrix(weight)
  colnames(weight) <- names(weight.list)[geneset.index]
  
  genes <- rownames(weight)
  geneset <- lapply(weight.list,names)
  membership <- table(unlist(geneset))
  anno_row <- as.data.frame(membership[genes])
  names(anno_row) <- "Membership"
  anno_row[anno_row>1,] <- "overlapping"
  anno_row[anno_row==1,] <- "exclusive"
  
  pheatmap(weight, scale="none",
           cellwidth = 60, cellheight = 20,
           color=colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")))(10), breaks=seq(0,1,.1),
           cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE, show_colnames=TRUE,
           annotation_row=anno_row, display_numbers=TRUE, main=main, ...)
}



############################################
## Fisher's method for combining p-values ##
############################################

fisher <- function(pvec){
  n <- length(pvec)
  chisq <- -2*sum(log(pvec))
  p <- pchisq(chisq, df=2*n, lower.tail=FALSE)
  return(p)
}
