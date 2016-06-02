## This file includes all auxilary functions that are not exported for
## a typical end user to use.  A determined user can still get these
## functions by scaleX <- FUNNEL:::scaleX.

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


