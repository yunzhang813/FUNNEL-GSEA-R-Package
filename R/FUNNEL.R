## Required packages: fda, quadrupen, pheatmap, RColorBrewer

###########################
## FUNCTION: scaleTime() ##
###########################

scaleTime <- function(time){
  new.time <- (time - min(time)) / diff(range(time))
  round(new.time, 2)
}

#########################
## FUNCTION: scaleX() ##
#########################

scaleX <- function(X){
  X.scale <- t(scale(t(X)))
  return(X.scale)
}

###############################
## FUNCTION: modifyGeneset() ##
###############################

modifyGeneset <- function(geneset, genenames){
  newGeneset <- lapply(geneset, function(z){z<-z[z %in% genenames]})
  return(newGeneset)
}

########################
## FUNCTION: getRho() ##
########################

getRho <- function(X, geneset)
{
  rho.geneset <- NULL
  for (i in 1:length(geneset)) {
    index <- geneset[[i]]
    cor.expr.geneset <- cor(t(X[index, ]))
    m <- length(index)
    rho.geneset <- rbind(rho.geneset, (m * mean(cor.expr.geneset) - 1)/(m - 1))
  }
  rho.hat <- mean(rho.geneset)
  return(rho.hat)
}

## example
# getRho(X,mGeneset)


##########################
## FUNCTION: getBasis() ##
##########################

getBasis <- function(time){
  Nt <- length(time)
  basis <- create.bspline.basis(range(time), Nt+4-2, 4, time)
  return(basis)
}


###########################
## FUNCTION: getFstats() ##
###########################

getFstats <- function(X,time,rr=rep(1,length(time)),selection_k="FVE",FVE_threshold=0.9)
  ### INPUT ###
  # X: n*m data matrix, with missing values denoted as NA.
  # time: length mm vector, unique time points.
  # rr: number of repetitions at each unique time point.
  # selection_k: the method of choosing the number of principal components;
  #              "FVE" (fraction of variance explained) : use scree plot
  #                           approach to select number of principal
  #                           components), see "FVE_threshold" below;
  #              positive integer K: user-specified number of principal components.
  # FVE_threshold: a positive number between 0 and 1; It is used with the option
  #                selection_k = "FVE" to select the number of principal components
  #                that explain at least "FVE_threshold" of total variation.
{
  y = X
  tt = time
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

## example
# Fstats <- getFstats(X, time, rr=rep(1,length(time)))

############################
## FUNCTION: smoothExpr() ## getBasis()
############################

smoothExpr <- function(X, time, lambda=10^-3.5)
  ### IMPORT ###
  # X = gene expression matrix, with rows = genes (IMPORTANT: row names use same gene annotation as geneset), columns = time points
  # time = vector of distinct time points
  # lambda = 10^-3.5, which is the smoothing penalty mannually selected for scaled-X and scaled-time
  ### EXPORT ###
  # fdobj = fdobj for the expression matrix
{
  basis <- getBasis(time)
  par <- fdPar(basis, 2, lambda=lambda)
  fdobj <- smooth.basis(time, t(X), par)$fd
  ## output
  return(fdobj)
}

## example
# fdExpr0 <- smoothExpr(X, time)


#############################
## FUNCTION: PCA.geneset() ##
#############################

PCA.geneset <- function(fdexpr, geneset, nharm=3, centerfns=FALSE)
{
  varprop <- coef.mat <- NULL
  for(testset in geneset){
    pca <-  pca.fd(fdexpr[testset], nharm=nharm, centerfns=centerfns)  
    harms <- pca$harmonics
    coef.mat <- cbind(coef.mat, harms$coefs)
    varprop <- rbind(varprop, pca$varprop)    
  }
  colnames(varprop) <- paste0("harm",1:nharm)
  colnames(coef.mat) <- paste0(names(geneset), ".", colnames(coef.mat))
  harmonics <- fd(coef.mat, fdexpr$basis)
  ## output
  return(list("harmonics"=harmonics, "varprop"=varprop))
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

#########################
## FUNCTION: getBeta() ## PCA.geneset(), equiv.regression()
#########################

getBeta <- function(fdexpr, gene.i, geneset.i, nharm=3, centerfns=FALSE, equiv.threshold=0.01, lam1=0.4, lam2=0.1)
{
  ## equivalent regression
  yfd <- fdexpr[gene.i]
  xfd <- PCA.geneset(fdexpr, geneset.i, nharm=nharm, centerfns=centerfns)$harmonics
  equiv <- equiv.regression(yfd, xfd, threshold = equiv.threshold)
  Y <- equiv$y
  X <- equiv$Xmat
  colnames(X) <- rep(names(geneset.i), each=nharm)
  ## fitting elastic net
  en <- elastic.net(X, Y, lambda1=lam1, lambda2=lam2, intercept=FALSE, normalize=FALSE)
  beta.en <- as.numeric(attributes(en)$coef)
  names(beta.en) <- colnames(X)
  ## output
  return(beta.en)
}

## example
# gene.i <- "2645"
# index <- which(sapply(geneset, function(z){gene.i %in% z}))
# geneset.i <- geneset[index]
# getBeta(fdexpr, gene.i, geneset.i,lam1=0.5, lam2=10^-2)


#################################
## FUNCTION: getWeightMatrix() ## getBeta()
#################################
getWeightMatrix <- function(fdexpr, geneset, nharm=3, centerfns=FALSE, equiv.threshold=0.01, lam1=0.4, lam2=0.1, verbose=FALSE){
  genenames <- fdexpr$fdnames$reps
  weight <- matrix(NA, length(genenames), length(geneset))
  rownames(weight) <- genenames
  colnames(weight) <- names(geneset)
  for (gene.i in genenames){
    ## find pathways that contain gene.i
    index <- names(which(sapply(geneset, function(z){gene.i %in% z})))

    ## the actual approach
    geneset.i <- geneset[index]
    ## if non-ovelapping gene, weight = 1
    if (length(index)==1) {weight[gene.i, index] <- 1} 
    else{
      ## if FUNNEL fails, assign uniform weight
      weight[gene.i, index] <- 1/length(index)
      tryCatch({
        beta <- getBeta(fdexpr, gene.i, geneset.i, nharm=nharm, centerfns=centerfns, equiv.threshold=equiv.threshold, lam1=lam1, lam2=lam2)
      },
      error = function(err){
        if(verbose==TRUE){
          message("Uniform weighting is applied when elastic net fails.")
          message(paste("Original error message:", err)) 
        } 
      })
      for(testset in index){
        weight[gene.i, testset] <- sum(beta[names(beta)==testset]^2)/sum(beta^2) 
      }
    }
  }
  ## if not a number, then zero
  weight[is.nan(weight)] <- 0
  return(weight)
}

## example
# weight.mat <- getWeightMatrix(fdexpr, geneset)


####################################################
## FUNCTION: weightedRankSumTestWithCorrelation() ##
####################################################

weightedRankSumTestWithCorrelation <- function(index,statistics,weight=NULL,correlation=0,df=Inf)
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
## FUNCTION: FUNNELtest() ## getWeight(), weightedRankSumTestWithCorrelation()
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
    test <- weightedRankSumTestWithCorrelation(index=unlist(testset), statistics=Fstats, weight=weight, correlation=rho, df=df)
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

FUNNEL.GSEA <- function(X, time, geneset, lambda=10^-3.5, rr=rep(1,length(time)), selection_k="FVE", FVE_threshold=0.9, nharm=3, centerfns=FALSE, equiv.threshold=0.01, lam1=0.4, lam2=0.1, alpha.level=0.05)
  ### INPUT ###
  # X = original expression matrix
  # time = origninal time points
  # geneset = original geneset database
{
  checkInputs(X,time,geneset)
  ## data preparation
  time <- scaleTime(time)
  X <- scaleX(X)
  geneset <- modifyGeneset(geneset,rownames(X))
  ## smoothing curve
  fdexpr <- smoothExpr(X, time, lambda=lambda)
  ## get Fstats and rho
  Fstats <- getFstats(X, time, rr=rr, selection_k=selection_k, FVE_threshold=FVE_threshold)
  rho.hat <- getRho(X, geneset)
  ## FUNNEL test
  FUNNEL.out <- FUNNELtest(fdexpr, geneset, Fstats=Fstats, rho=rho.hat, df=sum(rr)-1, nharm=nharm, centerfns=centerfns, equiv.threshold=equiv.threshold, lam1=lam1, lam2=lam2)
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
