## ------------------------------------------------------------------------
  wfNMF<-function(V,r,maxiter){
  V<-as.matrix(V)
  h<-nrow(V)
  u<-ncol(V)
  par(mar=c(1,1,1,1))
  par(mfrow=c(1,2))
  image(t(V[h:1,]),col = gray.colors(12))
  W=matrix(runif(h*r),nrow = h,ncol = r)
  H=matrix(runif(r*u),nrow = r,ncol = u)
  for (i in 1:maxiter) {
    W=W*(V%*%t(H))/(W%*%H%*%t(H))
    W=W/(matrix(rep(1,h),h,1)%*%colSums(W))
    H=H*((t(W)%*%V)/(t(W)%*%W%*%H))
  }
  img_V=W%*%H
  image(t(img_V[h:1,]),col = gray.colors(12))
    }

## ------------------------------------------------------------------------
library(SC19039)
wfNMF(V,r=20,maxiter=1000)

## ------------------------------------------------------------------------
wfNMF(V,r=40,maxiter=1000)

## ------------------------------------------------------------------------
SNMF<-function(V,r,maxiter,kong){
  V<-as.matrix(V)
  h<-nrow(V)
  u<-ncol(V)
  par(mar=c(1,1,1,1))
  par(mfrow=c(1,2))
  image(t(V[h:1,]),col = gray.colors(12))
  W<-matrix(runif(h*r),nrow = h,ncol = r)
  H<-matrix(runif(r*u),nrow = r,ncol = u)
  for (i in 1:maxiter) {
    W<-W*(V%*%t(H))/(W%*%H%*%t(H))
    W<-W/(matrix(rep(1,h),h,1)%*%colSums(W))
    H<-H*((t(W)%*%V)/(t(W)%*%W%*%H))
    y<-t(W)%*%W%*%H-t(W)%*%V
    for (k in 1:r) {
      y[k,]<-y[k,]/c(t(W[,k])%*%W[,k])
    }
    deta0<-0.5*sqrt((H-y)*(H-y))
    for (k in 1:r) {
      deta0[k,]<-deta0[k,]%*%(t(W[,k])%*%W[,k])
    }
    As<-order(deta0,decreasing = TRUE)
    deta0[As]
    H[deta0<deta0[As[kong]]]=0.00001
  }
  img_V<-W%*%H
  image(t(img_V[h:1,]),col = gray.colors(12))
  return(H[1:50])
}

## ------------------------------------------------------------------------
SNMF(V,r=20,maxiter=1000,kong=500)

## ------------------------------------------------------------------------
SNMF(V,r=40,maxiter=1000,kong=1000)

