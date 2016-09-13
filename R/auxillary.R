## This file includes all auxilary functions that are not exported for
## a typical end user to use.  A determined user can still get these
## functions by scaleX <- FUNNEL:::getFstats.


########################
## FUNCTION: getRho() ##
########################

getRho <- function(X, genesets)
{
  rho.genesets <- NULL
  for (i in 1:length(genesets)) {
    index <- genesets[[i]]
    cor.expr.genesets <- cor(t(X[index, ]))
    m <- length(index)
    rho.genesets <- rbind(rho.genesets, (m * mean(cor.expr.genesets) - 1)/(m - 1))
  }
  rho.hat <- mean(rho.genesets)
  return(rho.hat)
}


############################
## FUNCTION: smoothExpr() ##
############################

smoothExpr <- function(X, tt, lambda=10^-3.5)
### IMPORT ###
# X = gene expression matrix, with rows = genes (IMPORTANT: row names use same gene 
#     annotation as genesets), columns = time points
# tt = vector of distinct time points
# lambda = smoothing penalty. Defualt = 10^-3.5 is mannually selected for
#          standardized X and tt
### EXPORT ###
# fdobj = smoothed expression cuvres (a functional data object)
{
  basis <- create.bspline.basis(range(tt), length(tt)+4-2, 4, tt)
  par <- fdPar(basis, 2, lambda=lambda)
  fdobj <- smooth.basis(tt, t(X), par)$fd
  ## output
  return(fdobj)
}


#############################
## FUNCTION: PCA.genesets() ##
#############################

PCA.genesets <- function(fdexpr, genesets, nharm=3, centerfns=FALSE)
### IMPORT ###
# fdexpr = smoothed expression (a functional data object)
# genesets = gene sets
# nharm = number of harmonic functions.
# centerfns = if to center the functions for pca.fd()
### EXPORT ###
# harmonics = harmonic functions
# varprop = proportion of variance explained
{
  varprop <- coef.mat <- NULL
  for(testset in genesets){
    pca <-  pca.fd(fdexpr[testset], nharm=nharm, centerfns=centerfns)  
    harms <- pca$harmonics
    coef.mat <- cbind(coef.mat, harms$coefs)
    varprop <- rbind(varprop, pca$varprop)    
  }
  colnames(varprop) <- paste0("harm",1:nharm)
  colnames(coef.mat) <- paste0(names(genesets), ".", colnames(coef.mat))
  harmonics <- fd(coef.mat, fdexpr$basis)
  ## output
  return(list("harmonics"=harmonics, "varprop"=varprop))
}

#########################
## FUNCTION: getBeta() ## PCA.genesets(), equiv.regression()
#########################

getBeta <- function(fdexpr, gene.i, geneset.i, nharm=3, centerfns=FALSE, equiv.threshold=0.01, lam1=0.4, lam2=0.1)
{
  ## equivalent regression
  yfd <- fdexpr[gene.i]
  xfd <- PCA.genesets(fdexpr, geneset.i, nharm=nharm, centerfns=centerfns)$harmonics
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


#################################
## FUNCTION: getWeightMatrix() ## getBeta()
#################################

getWeightMatrix <- function(fdexpr, genesets, nharm=3, centerfns=FALSE, equiv.threshold=0.01, lam1=0.4, lam2=0.1, verbose=FALSE)
{
  genenames <- fdexpr$fdnames$reps
  weight <- matrix(NA, length(genenames), length(genesets))
  rownames(weight) <- genenames
  colnames(weight) <- names(genesets)
  for (gene.i in genenames){
    ## find pathways that contain gene.i
    index <- names(which(sapply(genesets, function(z){gene.i %in% z})))

    ## the actual approach
    geneset.i <- genesets[index]
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
          message(paste("Uniform weighting is applied because elastic-net regression failed for gene", gene.i, "in gene sets", paste(index, collapse=", ")))
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


############################################
## Fisher's method for combining p-values ##
############################################

fisher <- function(pvec){
  n <- length(pvec)
  chisq <- -2*sum(log(pvec))
  p <- pchisq(chisq, df=2*n, lower.tail=FALSE)
  return(p)
}
