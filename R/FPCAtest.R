## Required packages: MASS, akima

FPCAtest = function(y,tt,rr,H0="mean",B=1000,p.adjust.methods="BH",FDR=0.05,miss=0,
           delta="auto",bwxcov=c(0,0),ngrid=51,ngrid1=30,kern="gauss",error=1,
           selection_k="FVE",FVE_threshold=0.9,verbose="on")
### INPUT ###
# y: n*m data matrix, with missing values denoted as NA.
# tt: length mm vector, unique time points.
# rr: number of repetitions at each unique time point.
# H0: "mean" null hypothesis is that the time course data are equal to the mean;
#     "base" null hypothesis is that the time course data are equal to the baseline.
# B: number of permutations
# p.adjust.methods: multiple testing ajustment method to be used; options include:
#        c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none");
#        see ?p.adjust for more details.
# FDR: level of false discovery rate
# miss: 0-1 indicator for whether the data have missing values.
# delta: small nonnegative constant used in the modified F-statistics to stabilize
#        the variance; if is 0, the regular F statistics are used;
#        default is "auto", the program uses the estimated measurement error variance (see
#        output sigma in PCA).
# For the rest of input, please refer to function PCA
### OUTPUT ###
# PCAres: output from PCA
# p: vector of pvalues
# stat: vector of statistics
# id.sig: index for subjects that reject the null hypothesis
{
    res = PCA(y,tt,rr,miss,bwxcov,ngrid,ngrid1,kern,error,selection_k,FVE_threshold,verbose)
    if(is.null(delta)){ delta = res$sigma }
    r = permtest(res,y,tt,rr,miss,H0,B,delta)
    p = r$p
    stat = r$stat
    p.adj = p.adjust(p,method=p.adjust.methods)
    sig = which(p.adj<=FDR)
    list(PCAres=res,p=p,stat=stat,id.sig=sig)
}


permtest = function(res,y,tt,rr,miss=0,H0="mean",B=1000,delta="auto")
### INPUT ###
# res: output from PCA.
# y: n*m data matrix, with missing values denoted as NA.
# tt: length mm vector, unique time points.
# rr: number of repetitions at each unique time point.
# miss: 0-1 indicator for whether the data have missing values.
# H0: "mean" null hypothesis is that the time course data are equal to the mean;
#     "base" null hypothesis is that the time course data are equal to the baseline.
# B: number of permutations
# delta: small nonnegative constant used in the modified F-statistics to stabilize
#        the variance; if is 0, the regular F statistics are used;
#        default is "auto", the program uses the estimated measurement error variance (see
#        output sigma in PCA).
### OUTPUT ###
# p: vector of pvalues
# stat: vector of statistics
{
    n = nrow(y)
    m = ncol(y)
    mm = length(tt)
    rr1 = rep(1:mm,rr)
    mu = res$mu
    yfit = res$yfit_orig
    if(delta=="auto"){ delta = res$sigma }
    if(H0 == "mean"){
      ss0 = rowSums((y-matrix(rep(mu,m),,m))^2,na.rm=T)
      ss1 = rowSums((y-yfit)^2,na.rm=T)
      stat = (ss0-ss1)/(ss1+delta)
      ## ss = sort(stat,index.return=T)
      ss = sort(stat); oo <- order(stat)
      edge = c(min(stat)-100,ss,max(stat)+100)
      nbt = rep(0,length=length(edge)-1)
      for(i in 1:B){
        tmp = sample(1:mm)
        tmprr = rr[tmp]
        ind = sapply(1:mm,function(j) which(rr1==tmp[j]))
        samp = y[,unlist(ind)]
        sampfit = PCApred(res,samp,tt,tmprr,miss)$yfit_orig
        ss1bt = rowSums((samp-sampfit)^2,na.rm=T)
        bt = (ss0-ss1bt)/(ss1bt+delta)
        nbt = nbt+hist(bt,breaks=edge,plot=F)$count
      }
    }else{
      ss0 = rowSums((y-matrix(rep(yfit[,1],m),,m))^2,na.rm=T)
      ss1 = rowSums((y-yfit)^2,na.rm=T)
      stat = (ss0-ss1)/(ss1+delta)
      ss = sort(stat,index.return=T)
      edge = c(min(stat)-100,ss,max(stat)+100)
      nbt = rep(0,length=length(edge)-1)
      for(i in 1:B){
        tmp = sample(1:mm)
        tmprr = rr[tmp]
        ind = sapply(1:mm,function(j) which(rr1==tmp[j]))
        samp = y[,unlist(ind)]
        sampfit = PCApred(res,samp,tt,tmprr,miss)$yfit_orig
        ss0bt = rowSums((samp-matrix(rep(sampfit[,1],m),,m))^2,na.rm=T)
        ss1bt = rowSums((samp-sampfit)^2,na.rm=T)
        bt = (ss0bt-ss1bt)/(ss1bt+delta)
        nbt = nbt+hist(bt,breaks=edge,plot=F)$count
      }
    }
    p = 1-cumsum(nbt[-length(nbt)])/(n*B)
    ## p = p[sort(oo,index.return=T)$ix]
    p = p[order(oo)]

    list(p=p,stat=stat)
}


PCA = function(y,tt,rr,miss=0,bwxcov=c(0,0),ngrid=51,ngrid1=30,kern="gauss",error=1,
               selection_k="FVE",FVE_threshold=0.9,verbose="on")
### INPUT ###
# y: n*m data matrix, with missing values denoted as NA.
# tt: length mm vector, unique time points.
# rr: number of repetitions at each unique time point.
# miss: 0-1 indicator for whether the data have missing values.
# bwxcov: bandwidth for smoothing the covariance function;
#         default c(0,0), choose the bandwidth by generalized cross-validation (GCV).
# ngrid: number of support points for output time grid
# ngrid1: number of support points for the covariance surface in the GCV procedure.
# kern: a character string to define the kernel to use in the smoothing;
#       options include: "gauss" Gaussian kernel (default)
#                        "epan" Epanechnikov kernel
#                        "rect" Rectangular kernel
# error: 0-1 indicator for additional measurement errors.
# selection_k: the method of choosing the number of principal components;
#              "FVE" (fraction of variance explained) : use scree plot
#                           approach to select number of principal
#                           components), see "FVE_threshold" below;
#              positive integer K: user-specified number of principal components.
# FVE_threshold: a positive number between 0 and 1; It is used with the option
#                selection_k = "FVE" to select the number of principal components
#                that explain at least "FVE_threshold" of total variation.
# verbose: "on" display diagnostic messages; "off" suppress diagnostic messages.
### OUTPUT ###
# noeig: number of selected principal components.
# sigma: estimate of measurement error variance.
# lambda: estimated eigenvalues.
# phi: estimated eigenfunctions at original time points tt.
# eigens: estimated eigenfunctions at output time grid out21.
# xi: estimated principal component scores.
# mu: mean level of each subject (each row of y).
# xcov: smoothed covariance function, corresponding to grid out21.
# bw_xcov: bandwidth for smoothed covariance.
# xcovfit: fitted covariance function, based on truncated estimate of eigenvalues
#          ("lambda") and eigenfunctions ("eigens"), corresponding to out21.
# FVE: fraction of variance explained.
# yfit: fitted curves for each subject, corresponding to out21.
# yfit_orig: fitted curves evaluated at the same time points as original y.
# out21: output time grid.
{
  n = nrow(y)
  m = ncol(y)
  mm = length(tt)
  mu = rowMeans(y,na.rm=T)
  cy = y-matrix(rep(mu,m),nrow=n,ncol=m)
  # calculate raw covariance from centered y
  rcov = rawcov(cy,tt,rr,miss)
  # select bandwidth by GCV
  if(bwxcov[1] == 0 || bwxcov[2] == 0){
      bw_xcov = gcv_mullwlsn(tt,ngrid1,miss,error,kern,rcov,verbose)$bw_xcov
      if(any(is.na(bw_xcov))){
          cat("Error: FPCA is aborted because the observed data is too sparse to estimate the covariance function!\n")
             return(NULL)
          }
      bw_xcov = adjustBW2(kern,bw_xcov,1,0,miss,verbose)
  }else if(all(bwxcov > 0)){
      bw_xcov = bwxcov
  }else if(bwxcov[1] < 0 || bwxcov[2] < 0){
      cat("Error: Bandwidth choice for the covariance function must be positive!\n")
      return(NULL)
  }

  # smooth raw covariance
  out21 = seq(min(tt),max(tt),length=ngrid)
  rcov1 = rcov
  if(error == 1){
     tpairn = rcov1$tpairn
     tneq = tpairn[1,] != tpairn[2,]
     cyy = rcov1$cyy
     rcov1$tpairn = tpairn[,tneq]
     rcov1$cxxn = cyy[tneq]
     rcov1$win = rep(1,len = length(rcov1$cxxn))
     if(miss == 1){
        rcov1$count = rcov1$count[tneq]
     }
  }

  if(miss == 1){  #smooth raw covariance
     r = mullwlsk(bw_xcov,kern,rcov1$tpairn,rcov1$cxxn,rcov1$win,out21,out21,rcov1$count)
  }else{             #smooth raw covariance
     r = mullwlsk(bw_xcov,kern,rcov1$tpairn,rcov1$cxxn,rcov1$win,out21,out21)
  }
  invalid = r$invalid
  xcov = r$mu
  xcov = (xcov+t(xcov))/2

   # select number of components
   if(invalid == 0){
     r = no_FVE(xcov,FVE_threshold)
     no_opt = r$no_opt
     FVE = r$FVE
   }else{
     cat("FPCA is aborted because enough points to estimate the smooth covariance function!\n")
     return(NULL)
   }

   pc_options = c("FVE","user")
   if(selection_k == "FVE"){
       k_id = 1;
   }else if(is.numeric(selection_k) && selection_k > 0){
     no_opt = selection_k;
     k_id = 2;
   }else{
     cat(paste('"selection_k" must be a positive integer! Reset to "FVE" method with threshold =', FVE_threshold, '!\n'))
     k_id = 1;
   }

   if(verbose == "on"){
        cat(paste("Best number of principal components selected by ", pc_options[k_id]," : ", no_opt, ".\n", sep = ""))
        if(k_id != 1){
            cat(paste("It accounts for ", round(FVE[no_opt], digits = 4)*100, "% of total variation.\n", sep = ""))
        }else{
            cat(paste("It accounts for ", round(FVE[no_opt], digits = 4)*100, "% of total variation (threshold = ", FVE_threshold, ").\n", sep = ""))
        }
        cat(paste("FVE calculated from ", ngrid, " possible eigenvalues: \n", sep = ""))
        print(FVE)
   }

   # compute the eigenvalues, eigenfunctions and fitted covariance
   if(error == 1){
     r = pc_covE(tt,bw_xcov,ngrid,1,kern,rcov)
     invalid = r$invalid
     sigma = r$sigma

     if(invalid == 0){
         r = getEigens(xcov,tt,out21,no_opt,1)
         lambda = r$lambda
         phi = r$phi
         eigens = r$eigens
         if(length(lambda)==1){
             xcovfit = lambda*eigens%*%t(eigens)
         }else{
             xcovfit = eigens%*%diag(lambda)%*%t(eigens)
         }
         rm(r)
     }else{
         lambda = NULL; phi = NULL; eigens = NULL; xcovfit = NULL
     }
   }else{
     sigma = NULL
     r = getEigens(xcov,tt,out21,no_opt,1)
     lambda = r$lambda
     phi = r$phi
     eigens = r$eigens
     xcovfit = eigens%*%diag(lambda)%*%t(eigens)
     rm(r)
  }

  # compute the principal component scores
  phi1 = sapply(1:no_opt,function(i) rep(phi[,i],rr))
  w = ginv(t(phi1)%*%phi1)
  if(miss == 1){
    if(no_opt == 1){
      xi = matrix(sapply(1:n,function(i) cy[i,!is.na(y[i,])]%*%phi1[!is.na(y[i,]),]%*%w),nrow=n)
    }else{
      xi = t(sapply(1:n,function(i) cy[i,!is.na(y[i,])]%*%phi1[!is.na(y[i,]),]%*%w))
    }
  }else{
    xi = cy%*%phi1%*%w
  }
  yfit_orig = matrix(rep(mu,m),,m)+xi%*%t(phi1)
  yfit = matrix(rep(mu,ngrid),,ngrid)+xi%*%t(eigens)

  list(noeig=no_opt,sigma=sigma,lambda=lambda,phi=phi,eigens=eigens,xi=xi,mu=mu,
       xcov=xcov,bw_xcov=bw_xcov,xcovfit=xcovfit,FVE=FVE,yfit=yfit,
       yfit_orig=yfit_orig,out21=out21)
}


PCApred = function(res,y,tt,rr,miss)
# res: output of PCA
{
    n = nrow(y)
    m = ncol(y)
    mm = length(tt)
    rr1 = rep(1:mm,rr)

    noeig = res$noeig
    out21 = res$out21
    phi = res$phi
    eigens = res$eigens

    mu = rowMeans(y,na.rm=T)
    cy = y-matrix(rep(mu,m),,m)
    phi1 = sapply(1:noeig,function(i) rep(phi[,i],rr))
    w = ginv(t(phi1)%*%phi1)
    if(miss == 1){
      #xi = t(sapply(1:n,function(i) cy[i,!is.nan(y[i,])]%*%phi1[!is.nan(y[i,]),]%*%w))
      if(noeig == 1){
        xi = matrix(sapply(1:n,function(i) cy[i,!is.na(y[i,])]%*%phi1[!is.na(y[i,]),]%*%w),nrow=n)
      }else{
        xi = t(sapply(1:n,function(i) cy[i,!is.na(y[i,])]%*%phi1[!is.na(y[i,]),]%*%w))
      }
    }else{
      xi = cy%*%phi1%*%w
    }
    yfit_orig = matrix(rep(mu,m),,m)+xi%*%t(phi1)
    yfit = matrix(rep(mu,length(out21)),,length(out21))+xi%*%t(eigens)

    list(yfit=yfit,yfit_orig=yfit_orig,xi=xi)
}


adjustBW2 = function(kern=c("gauss","epan","rect","quar"),bw_xcov,npoly=nder+1,nder=0,miss=1,verbose="on")
{
    kern = kern[1]
    if(kern == "gauss"){
        if(miss == 0){
          bwxcov_fac = c(1.1,0.8,0.8)
        }else{
          bwxcov_fac = c(1.1, 1.2, 2)
        }
       if(nder > 2){
          facID = 3
       }else if(nder >= 0 && nder <= 2){
          facID = nder+1
       }else{
          facID = 1
       }
       bw_xcov = bw_xcov*bwxcov_fac[facID]
       if(verbose == "on")
         cat("Adjusted GCV bandwidth choice for COV function (npoly = ", npoly, "): (", bw_xcov[1],",", bw_xcov[2], ")\n")
    }
    bw_xcov
}


gcv_mullwlsn = function(tt,ngrid=30,miss=1,error=1,kern=c("gauss","epan","rect","quar"),rcov,verbose="on")
{
   kern = kern[1]
   out1 = tt
   a0 = min(out1)
   b0 = max(out1)

   h0 = getMinb(tt,miss)
   if(kern == "gauss"){
      if(is.na(h0)){
         h0 = b0
      }
      h0 = h0*0.2
   }
   if(is.na(h0)){
       cat("Error: the data is too sparse, no suitable bandwidth can be found! Try Gaussian kern instead!\n")
       return(list(bw_xcov = NA, gcv = NA))
   }

   rcovcount = rcov$count
   if(error == 1){
      tpairn = rcov$tpairn
      tneq = tpairn[1,]!=tpairn[2,]
      cyy = rcov$cyy
      tpairn = tpairn[,tneq]
      cxxn=cyy[tneq]
      win= rep(1,len = length(cxxn))
      if(miss == 1){
          rcovcount = rcovcount[tneq]
      }
   }else{
       tpairn = rcov$tpairn
       cxxn = rcov$cxxn
       win = rcov$win
   }
   rm(rcov,cyy,tneq)
   N = length(cxxn)
   r = b0-a0
   rm(out1)

   qq = (r/(4*h0))^(1/9)
   bw = sort(qq^(0:9)*h0)
   bw = matrix(rep(bw,2),,2)
   k0 = mykern(0,kern)
   out21 = seq(a0,b0,len=ngrid)

   leave = 0
   nleave = 0
   tooSparse = 0

   while(leave == 0){
       gcv = rep(Inf,nrow(bw))
       for(k in 1:nrow(bw)){
          if(miss == 1){
            xcov = mullwlsk(bw=bw[k,],kern=kern,xin=tpairn,yin=cxxn,win=win,out1=out21,out2=out21,count=rcovcount)
          }else{
            xcov = mullwlsk(bw=bw[k,],kern=kern,xin=tpairn,yin=cxxn,win=win,out1=out21,out2=out21)
          }
          invalid = xcov$invalid
          xcov = xcov$mu

          if(invalid != 1){
               o21 = expand.grid(out21,out21)
               xcov = as.vector(xcov)
               #require package akima for 2-D interpolation
               #it seems to me that only linear interpolation works
               #do not allow extrapolation
  require(akima)
               newxcov =interpp(o21[,1],o21[,2],xcov,tpairn[1,],tpairn[2,])$z

               rm(xcov)
               if(miss == 1){
                   cvsum = sum((cxxn/rcovcount-newxcov)^2)
               }else{
                   cvsum = sum((cxxn-newxcov)^2)
               }
               rm(newxcov)
               bottom = 1-(1/N)*((r*k0)/bw[k])^2

               gcv[k] = cvsum/(bottom)^2
               tmp = gcv[gcv != Inf]
               if(length(tmp) > 1 && gcv[k] > gcv[k-1]){
                    leave = 1
                    break
               }
          }
       }

       if(all(gcv == Inf)){
           if(nleave == 0 && bw[10,1] < r){
                bw_xcov = bw[10,]
                tooSparse = 1
           }else{
                cat("Error: the data is too sparse, no suitable bandwidth can be found! Try Gaussian kern instead!\n");
                return(list(bw_xcov = NA, gcv = NA))
           }
       }else{
           bw_xcov = bw[which(gcv == min(gcv))[1],]
       }

       if(bw_xcov[1] == r){
           leave = 1
           cat("data is too sparse, optimal bandwidth includes all the data!You may want to change to Gaussian kern!\n")
       }else if(bw_xcov[1] == bw[10,1] && nleave == 0){
           if((tooSparse == 1) || (sum(gcv == Inf) == 9)){
               cat("data is too sparse, retry with larger bandwidths!\n")
               h0 = bw[10,1]*1.01
           }else{
              cat("Bandwidth candidates are too small, retry with larger choices now!\n")
              h0 = bw[9,1]
           }
           newr = seq(0.5,1,by = 0.05)*r
           id = which(h0 < newr)[1]
           qq = (newr[id]/h0)^(1/9)
           bw = sort(qq^(0:9)*h0);
           bw = matrix(rep(bw,2),,2)
           if(verbose == "on"){
              cat("New bwxcov candidates:\n")
              print(bw)
           }
       }else if(bw_xcov[1] < bw[10,1] || nleave > 0){
           leave = 1
       }
       nleave = nleave+1
  }

    if(kern != "gauss" && verbose == "on")
       cat(paste("GCV bandwidth choice for COV function : (", bw_xcov[1],",",bw_xcov[2],")\n", sep = ""))

    return(list(bw_xcov=bw_xcov,gcv=gcv))
}


getEigens = function(xcov,out1,out21,noeig,disp=FALSE)
{
  r = range(out21)
  h = (r[2]-r[1])/(length(out21)-1)

  ngrid = nrow(xcov)

  r = eigen(xcov)
  eigens = r$vectors
  d = r$values
  rm(r)

  idx = which(Im(d)!=0)      #find indices for imaginary eigenvalues
  if(length(idx) >  0){
      stop(paste(length(idx),"eigenvalues are complex. The estimated auto-covariance surface is not symmetric!"))
  }
  idx = which(d <= 0)
  if(length(idx) > 0){
       if(disp)
          cat(paste(length(idx), "real eigenvalues are negative or zero and are removed!\n"))
       eigens = eigens[,d>0]
       d = d[d>0]
  }

  if(noeig > length(d)){
    noeig = length(d)
    cat(paste("At most",noeig,"number of PC can be selected!\n"))
  }

  eigens = eigens[,1:noeig]
  if(!is.matrix(eigens))
     eigens = matrix(eigens,length(out21),noeig)
  d = d[1:noeig]
  eigens = eigens/sqrt(h)
  lambda = h*d

  for(i in 1:noeig){
      eigens[,i] = eigens[,i]/sqrt(romb2(out21,eigens[,i]^2))
      if(eigens[2,i] < eigens[1,i])
        eigens[,i] = -eigens[,i]
  }

  #interpolate from the normalized the eigenfunctions
  phi = interp11(out21,eigens,out1)

  #normalize smoothed eigenfunctions
  for(i in 1:noeig){
      phi[,i] = phi[,i]/sqrt(romb2(out1,phi[,i]^2))
  }
  list(lambda=lambda,phi=phi,eigens=eigens,noeig=noeig)
}


getMinb = function(tt,miss=1,npoly=1)
{
   if(miss == 1){
      dstar = minb(tt,1+npoly)*2;
   }else{
      dstar = minb(tt,2+npoly)*1.5;
   }
   dstar
}


interp1 = function(x,y,newx = x,method="natural",nder = 0)
{
#use splinefun() in the stats package of R to perform 1-D interpolation
#default is used natual cubic spline, other possible values, see
#the description for splinefun().
    f = splinefun(x,y,method = method)
    if(is.matrix(newx)){
      res = f(newx[,1],deriv=nder)
      cols = ncol(newx)
      for(i in 2:cols){
        res = cbind(res,f(newx[,i],deriv=nder))
      }
   }else{
     res = f(newx,deriv=nder)
   }
   res
}


interp11 = function(xx,y,newx=xx,method="natural",nder = 0)
{
#xx is a vector, where all rows of y are evaluated at
#y is a matrix
   apply(y,2,function(x) interp1(xx,x,newx,method=method,nder = nder))
}


lwls = function(bw,kern=c("gauss","epan","rect","quar"),nwe=0,npoly=nder+1,nder=0,xin,yin,win,xou,bwmuLocal=0)
{
   if(npoly < nder)
      stop("Degree of polynomial should be no less than the order of derivative!")

   kern = kern[1]
   require(MASS)
   actobs = which(win!=0)
   xin = xin[actobs]
   yin = yin[actobs]
   win = win[actobs]
   invalid = 0

   aa = 1
   if(nwe == -1)
      aa = 4
   mu = numeric(length(xou))
   gap = mu

   search_ht = function(t1,h0){
     tmp = -(t1-75)^2/6000
     h0*(1+exp(tmp))
   }

   if(bw > 0){
      if(bwmuLocal){
         bw = search_ht(xou,bw)
      }else{
         bw = rep(bw,length(xou))
      }
   }else{
      stop("Bandwidth choice for mu(t) and/or its derivative must be positive!")
   }

   #LWLS with different weight functions
   for(i in 1:length(xou))
   {
       #(3-1) Locating local window
       if(kern != "gauss" && kern != "gausvar"){
          idx = xin <= xou[i] + aa*bw[i] & xin >= xou[i]-aa*bw[i]
       }else{
          idx = 1:length(xin)
       }
       lx = xin[idx]
       ly = yin[idx]
       lw = win[idx]

       if(length(unique(lx)) >= (npoly+1)){

          #Sepcify weight matrix
          llx = (lx-xou[i])/bw[i]

          if(kern == "epan"){
              w = lw*(1-llx^2)*0.75
          }else if(kern == "rect"){
              w = lw
          }else if(kern == "optp"){
              w = lw*(1-llx^2)^(nwe-1)
          }else if(kern == "gauss"){
              w = lw*dnorm(llx)
          }else if(kern == "gausvar"){
              w = lw*dnorm(llx)*(1.25-0.25*llx^2)
          }else if(kern == "quar"){
              w = lw*(15/16)*(1-llx^2)^2
          }else{
              cat("Invalid kernel, Epanechnikov kernel is used!\n")
              w = lw*(1-llx^2)*0.75
          }
          W = diag(w)
          # Define design matrix
          dx = matrix(1,length(lx),npoly+1)
          for(j in 1:npoly){
            dx[,j+1] = (xou[i]-lx)^j
          }

          p = ginv(t(dx)%*%W%*%dx)%*%t(dx)%*%W%*%ly

          #Find estimate
          mu[i] = p[(nder+1)*gamma(nder+1)*((-1)^nder)]
       }else{
          gap[i] = 1
          invalid = 1
       }
   }
   indx = which(gap == 0)
   if((length(indx)>= 0.9*length(xou))&& (length(indx)<length(xou))){
        mu1 = mu[indx]
        rr = myunique(xou[indx])
        mu = interp1(rr$out1,mu1[rr$id],xou)
   }else if(length(indx) < 0.9*length(xou)){
        mu = NULL
        cat("Too many gaps, please increase bandwidth!\n")
        invalid = 1
   }
   return(list(invalid = invalid, mu = mu))
}


minb = function(x,numPoints=2)
{
   n = length(x)
   x = sort(x)
   if(numPoints > 1){
     max(x[numPoints:n]-x[1:(n-numPoints+1)])
   }else{
     max((x[2:n]-x[1:(n-1)])/2)
   }
}


mullwlsk = function(bw,kern=c("gauss","epan","rect","quar"),xin,yin,win=rep(1,length(xin)),out1,out2,count=NULL)
{
    require(MASS) #for generalized inverse
    kern = kern[1]
    active = which(win!= 0)
    xin = xin[,active]
    yin = yin[active]
    win = win[active]
    invalid = 0
    mu = matrix(NA,length(out2),length(out1))
    for(i in 1:length(out2)){
        for(j in i:length(out1)){
           #locating local window
           if(kern != "gauss"){
               list1 = (xin[1,] >= out1[j]-bw[1]-10^(-6)) & (xin[1,] <= out1[j] + bw[1]+10^(-6))
               list2 = (xin[2,] >= out2[i]-bw[2]-10^(-6)) & (xin[2,] <= out2[i] + bw[2]+10^(-6))
               ind = list1 & list2
           }else{
               ind = !logical(dim(xin)[2])
           }
           lx = xin[,ind]
           ly = yin[ind]
           lw = win[ind]
           #computing weight matrix
           if(dim(unique(t(lx)))[1]>=3){     #at least 3 unique number of time pairs in the local window
               llx = rbind((lx[1,]-out1[j])/bw[1],(lx[2,]-out2[i])/bw[2])
               #deciding the kern used
               k = dim(llx)[2]
               if(kern == "epan"){
                   temp = lw*(1-llx[1,]^2)*(1-llx[2,]^2)*(9/16)
               }else if(kern == "rect"){
                   temp = lw*rep(1,dim(lx)[2])/4
               }else if(kern == "gauss"){
                   temp = lw*dnorm(llx[1,])*dnorm(llx[2,])
               }else if(kern == "gausvar"){
                   temp = lw*dnorm(llx[1,])*(1.25-0.25*llx[1,]^2)*(dnorm(llx[2,])*(1.5-0.5*llx[2,]^2))
               }else if(kern == "quar"){
                   temp = lw*((1-llx[1,]^2)^2)*((1-llx[2,]^2)^2)*(225/256)
               }
               W = diag(temp)

               #computing design matrix
               X = matrix(1,length(ly),3)
               X[,2] = t(lx[1,])-out1[j]
               X[,3] = t(lx[2,])-out2[i]
               if(!is.null(count)){
                   temp = temp*count[ind]
                   W1 = diag(temp)
               }else{
                   W1 = W
               }

               beta = ginv(t(X)%*%W1%*%X)%*%t(X)%*%W%*%ly
               rm(X,W,W1)
               mu[i,j] = beta[1]
           }else{
               invalid = 1
               cat("Not enough points in local window, please increase bandwidth!\n")
               return(list(invalid = invalid, mu = NULL))
           }
        }
    }

    if(!is.null(mu)){
       #mu[lower.tri(mu)] = mu[upper.tri(mu)] #assign lower triangular part of the mu matrix
                                             #to be the same as the upper triangular part
       a = matrix(0, dim(mu)[1], dim(mu)[2]);
       a[upper.tri(a)] = mu[upper.tri(mu)]
       mu = diag(mu)*diag(1, dim(mu)[1],dim(mu)[2])+a+t(a)
    }
    return(list(invalid = invalid, mu = mu))
}


mykern = function(x,kern="gauss")
{
   if(kern == "quar"){
      (15/16)*(abs(x) <= 1)*(1-x^2)^2
   }else if(kern == "epan"){
      0.75*(abs(x)<=1)*(1-x^2)
   }else if(kern == "rect"){
      0.5*(abs(x)<=1)
   }else if(kern == "gausvar"){
      (1/sqrt(2*pi))*exp(-0.5*x^2)*(1.25-0.25*x^2)
   }else if(kern == "gausvar1"){
      (1/sqrt(2*pi))*exp(-0.5*x^2)*(1.5-0.5*x^2)
   }else if(kern == "gausvar2"){
      k1 = (1/sqrt(2*pi))*exp(-0.5*x^2)*(1.25-0.25*x^2);
      k2 = (1/sqrt(2*pi))*exp(-0.5*x^2)*(1.5-0.5*x^2);
      k1*k2
   }else{
      (1/sqrt(2*pi))*exp(-0.5*x^2);
   }
}


myunique = function(xx,sorted=TRUE)
{
#get the unique value or unique sorted value of input vector "xx"
#any missing values from "xx" are ignored
#out1[id] is same as "xx" (without missing values)
   xx = xx[!is.na(xx)]
   if(sorted){
     out1 = unique(sort(xx))
   }else{
     out1 = unique(xx)
   }
   id = sapply(xx,function(x) which(out1==x))
   list(out1 = out1, id = id)
}


no_FVE = function(xcov,FVE_threshold=0.9,disp=FALSE)
{
   d = eigen(xcov,only.values=FALSE)$values
   idx = which(Im(d)!=0)      #find indices for imaginary eigenvalues
   if(length(idx) >  0){
      stop(paste(length(idx), "eigenvalues are complex. The estimated auto-covariance surface is not symmetric!"));
   }
   idx = which(d <= 0)
   if(length(idx) > 0){
       if(disp){
         cat(paste(length(idx), "real eigenvalues are negative or zero and are removed!\n"))
       }
       d = d[d > 0]
   }
   lambda = sort(d,decreasing=TRUE)
   FVE = cumsum(lambda)/sum(lambda)
   no_opt = which(FVE > FVE_threshold)[1]
   list(no_opt = no_opt,FVE = FVE,lambda = lambda)
}


pc_covE = function(tt,bw_xcov,ngrid=51,cut=1,kern=c("gauss","epan","rect","quar"),rcov,npoly=1)
{
    kern = kern[1]
    out1 = tt
    a0 = min(out1)
    b0 = max(out1)
    lint = b0-a0
    h = (b0-a0)/(ngrid-1)
    out21 = seq(a0,b0,len=ngrid)
    out22 = out21

    tpairn = rcov$tpairn
    tneq = tpairn[1,] != tpairn[2,]
    cyy = rcov$cyy

    #This is for the case when miss = 1, the raw covariance
    #matrix needs to be divided by the number of individual sums
    #for each element in the matrix. In the case of miss = 0,
    #the division is n for each of the element.
    if(!is.null(rcov$count))
       cyy = cyy/rcov$count

    cxx = cyy[tneq]
    rm(rcov)

    win1 = rep(1,length(cxx))
    #get smoothed variance function for y(t) using lwls
    teq = !tneq
    vyy = cyy[teq]
    win2 = rep(1,length(vyy))
    #yvar is the variance function

    r = lwls(bw_xcov[1],kern,1,npoly,0,tpairn[1,teq],vyy,win2,out21,0)
    yvar = r$mu
    invalid = r$invalid

    if(invalid == 0){
       #estimate variance of measurement error term
       #use quadratic form on diagonal to estimate Var(x(t))
       r = rotate_mlwls(bw_xcov,kern,tpairn[,tneq],cxx,win1,rbind(out21,out22),npoly)
       invalid = r$invalid
       xvar = r$mu
       rm(r)

       if(invalid == 0){
           if(cut == 0){
              sigma = romb(out21,yvar-xvar)/lint
           }else if(cut == 1){
              a = a0+lint*0.25
              b = a0+lint*0.75
              ind1 = out21 > a & out21 < b
              yvar1 = yvar[ind1]
              xvar1 = xvar[ind1]
              sigma = romb(out21[ind1],yvar1-xvar1)*2/lint
           }
       }
    }
    if(sigma < 0){
       cat("Estimated sigma is negative, reset to zero now!\n")
       sigma = 0
    }
    list(invalid = invalid,sigma = sigma,xvar = xvar,yvar = yvar)
}


rawcov = function(cy,tt,rr,miss=1)
# cy: centered y
{
  n = nrow(cy)
  mm = length(tt)
  rr1 = rep(1:mm,rr)
  count = NULL

  if(miss == 0){
    cy1 = matrix(NA,n,mm)
    for(i in 1:mm){
      ind = which(rr1==i)
      if(length(ind) == 1){
        cy1[,i] = cy[,ind]
      }else{
        cy1[,i] = rowMeans(cy[,ind])
      }
    }
    cyy = t(cy1)%*%cy1/n
    cyy = as.vector(cyy)
    cxxn = cyy
    tpairn = rbind(rep(tt,each=mm),rep(tt,mm))
    win = rep(1,length(cxxn))
  }else{
    cy1 = matrix(NA,n,mm)
    for(i in 1:mm){
      ind = which(rr1==i)
      if(length(ind) == 1){
        cy1[,i] = cy[,ind]
      }else{
        cy1[,i] = rowMeans(cy[,ind])
      }
    }
    ID = matrix(1,n,mm)
    for(i in 1:n){
          id = which(is.na(cy1[i,]))
          cy1[i,id] = 0
          ID[i,id] = 0
    }
    count = t(ID)%*%ID
    count = as.vector(count)
    idx = count!=0
    count = count[idx]
    cyy = t(cy1)%*%cy1
    cyy = as.vector(cyy)
    cyy = cyy[idx]
    cxxn = cyy
    tpairn = rbind(rep(tt,each=mm),rep(tt,mm))
    tpairn = tpairn[,idx]
    win = rep(1,length(cxxn))
  }
  list(tpairn=tpairn,cxxn=cxxn,win=win,cyy=cyy,count=count)
}


romb = function(x,y,decdigs=10)
{
  rom = matrix(0,2,decdigs)
  romall = numeric(2^(decdigs-1)+1)
  a = min(x)
  b = max(x)
  romall = interp1(x,y,seq(a,b,by = (b-a)/2^(decdigs-1)))
  h = b-a
  rom[1,1] = h*(romall[1]+romall[length(romall)])/2
  for(i in 2:decdigs){
      st = 2^(decdigs-i+1)
      rom[2,1] = (rom[1,1]+h*sum(romall[st/2+seq(1,2^(decdigs-1),by=st)]))/2
      for(k in 1:(i-1)){
         rom[2,k+1] = ((4^k)*rom[2,k]-rom[1,k])/((4^k)-1)
      }
      rom[1,1:i] = rom[2,1:i]
      h = h/2
  }
  res = rom[1,decdigs]
  return(res)
}


romb2 = function(x,y,decdigs = 10)
{
#perform romberg integration for each column of y
  if(is.matrix(y)){
     m = dim(y)[2]
  }else{
     res = romb(x,y,decdigs)
     return(res)
  }
  res = numeric(m)
  for(i in 1:m){
    res[i] = romb(x,y[,i],decdigs)
  }
  res
}


rotate_mlwls = function(bw,kern =c("gauss","epan","rect","quar"),xin,yin,win,d,npoly = 1)
{
   active = win!=0
   xin = xin[,active]
   yin = yin[active]
   win = win[active]

   #rotating coordinates of predictors by pi/4
   R = sqrt(2)/2*matrix(c(1,-1,1,1),2,2,byrow = TRUE)
   xn = R%*%xin
   yn = yin
   dn = R%*%d

   invalid = 0
   #minimizing local weighted least squares
   m = ncol(d)
   mu = numeric(m)
   for(i in 1:m){
      #locating local window
      if(kern != "gauss" && kern != "gausvar"){
          list1 = which(xn[1,] >= dn[1,i]-bw[1] & xn[1,] <= dn[1,i]+bw[1])
          list2 = which(xn[2,] >= dn[2,i]-bw[2] & xn[2,] <= dn[2,i]+bw[2])
          ind = intersect(list1,list2)
      }else{
          ind = 1:ncol(xn)
      }

      lx = xn[,ind]
      ly = yn[ind]
      lw = win[ind]
      #computing weight matrix
      if(length(ly) >= npoly+1){
           llx = rbind((lx[1,]-dn[1,i])/bw[1],(lx[2,]-dn[2,i])/bw[2])
           #deciding the kernel to use
           if(kern == "epan"){
              w = lw*(1-llx[1,]^2)*(1-llx[2,]^2)*(9/16)
           }else if(kern == "rect"){
              w = lw*rep(1,dim(lx)[2])/4
           }else if(kern == "gauss"){
              w = lw*dnorm(llx[1,])*dnorm(llx[2,])
           }else if(kern == "gausvar"){
              w = lw*dnorm(llx[1,])*(1.25-0.25*llx[1,]^2)*dnorm(llx[2,])*(1.5-0.5*llx[2,]^2)
           }else if(kern == "quar"){
              w = lw*((1-llx[1,]^2)^2)*((1-llx[2,]^2)^2)*(225/256)
           }
           W = diag(w)
           #computing design matrix
           X = matrix(1,length(ly),3)
           #X[,1] = rep(1,len = length(ly))
           X[,2] = (lx[1,]-dn[1,i])^2
           X[,3] = lx[2,]-dn[2,i]

           beta = ginv(t(X)%*%W%*%X)%*%t(X)%*%W%*%ly
           mu[i] = beta[1]
      }else if(length(ly) == 1){
           mu[i] = ly
      }else{
           cat("No points in local window, please increase bandwidth!\n")
           invalid = 1;
           return(list(invalid = invalid, mu = NULL))
      }
   }
   list(invalid = invalid,mu = mu)
}


xeig = function(tt,lint=c(0,1),numPhi=3)
{
   if(length(lint) == 1){
      lb = 0
      T = lint
   }else{
      lb = lint[1]
       T = lint[2]
   }

   id = which(tt < lb | tt > T)
   if(length(id) > 0){
      cat("tt must be in [",lb,",",T,"]. Invalid elements in tt are removed!\n")
      tt = tt[-id]
   }

   phi = matrix(NA,numPhi, length(tt))

   if(numPhi == 1){
      phi = -sqrt(2/T)*cos(2*pi*tt/T)
      phi = matrix(phi, 1, length(phi))

   }else if(numPhi == 2){
      phi[1,] = -sqrt(2/T)*cos(2*pi*tt/T)
      phi[2,] = sqrt(2/T)*sin(2*pi*tt/T)
   }else{
      phi = matrix(0, numPhi, length(tt))
      id = 1:numPhi
      oddID = id%%2
      oddFactor = 1:sum(oddID)
      evenID = oddID==0
      evenFactor = 1:sum(evenID)

      phiOdd = matrix(0,sum(oddID),length(tt))
      phiEven = matrix(0,sum(evenID),length(tt))
      for(i in 1:sum(oddID)){
         phiOdd[i,] = -sqrt(2/T)*cos(2*oddFactor[i]*pi*tt/T)
      }
      phi[which(oddID==1),] = phiOdd

      for(i in 1:sum(evenID)){
         phiEven[i,] = sqrt(2/T)*sin(2*evenFactor[i]*pi*tt/T)
      }
      phi[which(evenID==1),] = phiEven
   }
   phi
}
