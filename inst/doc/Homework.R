## ------------------------------------------------------------------------
x1<-c(1.2,2.6, 3.3,4.5,6.1,6.9,7.5,8.4,9.8,10.1,11.2,12.8,13.5,14.4,15.4,16.7)
x2<-c(0.98,0.87,0.85,0.80,0.74,0.69,0.65,0.60,0.52,0.49,0.44,0.39,0.28,0.19,0.15,0.08)
x3<-c(12.6,13.1,14,15.1,16.8,17.2,18.3,19.1,20.1,21.2,22.3,23.5,24.6,25.4,26.7,27.8)
y<-c(53,56,57,58.9,59.9,60.4,61.9,62.8,63.9, 65, 66.4,68.2,69.9,70.8,72.2,74)
fun1<-lm(y~I(x1^2)+x1+x2+x3)
fun1
summary(fun1)$coef 

## ----islands-------------------------------------------------------------
knitr::kable(head(islands))

## ----rock----------------------------------------------------------------
knitr::kable(head(rock))

## ----echo=FALSE----------------------------------------------------------
par(mar=c(1,1,1,1))
par(mfrow=c(2,2))
plot(fun1)

## ------------------------------------------------------------------------
x4=matrix(rnorm(100,mean=1.5,sd=2),ncol=2)
plot(x4)

## ------------------------------------------------------------------------
#Clear variable
rm(list = ls())
# Write a function to generate n random vectors from MVN(mu, Sigma)
wf1 <- function(n, mu, Sigma) {
 # dimension is inferred from mu and Sigma 
  di <- length(mu)
  #Eigenvalue decomposition for sigma
  eva <- eigen(Sigma, symmetric = TRUE)
  #Get the value of the eigenvalue
  lam <- eva$values
  #Get the vector of the eigenvalue
  Ve <- eva$vectors
  #Here we use sqrt(lambda) to make sure Z
  C <- Ve %*% diag(sqrt(lam)) %*% t(Ve)
  Y <- matrix(rnorm(n*di), nrow = n, ncol = di)
  #X is the matrix we need
  X <- Y %*% C + matrix(mu, n, di, byrow = TRUE)
}

## ------------------------------------------------------------------------
#Initial value
mu <- c(0, 0)
#Preparing to draw four graphs
par(mar=c(1,1,1,1))
par(mfrow=c(2,2))
# To get four graphs,we need a for loop.
for (i in c(1,3,8,15)){
  Sigma <- matrix(c(i, 0, 0, i), nrow = 2, ncol = 2)
  X3 <- wf1(1000, mu, Sigma)
  cov(X3);#Covariance
  a3<-X3^2;
  #calculate value
  a<-sqrt(a3[,1] +a3[,2]);
  #Set the abscissa
  break3=seq(0,20,0.3);
  #Drawing histogram
  x1=break3
  hist(a,prob=TRUE,breaks=break3,main = i,col = 'yellow');
  #Density function curve
  lines(x1,x1*exp(-x1*x1/2/i)/(i))
}

## ------------------------------------------------------------------------
n<-1e4
#Generating a normal distribution
y1 <- rnorm(n,mean=0)
y2 <- rnorm(n,mean=3)
#Assign value to p1
r <- sample(c(0,1),n,replace=TRUE,prob = c(0.25,0.75))
z1 <- r*y1+(1-r)*y2
x2=seq(-5,10,0.2)
#Draw a histogram
hist(z1,prob=TRUE,breaks = seq(-5,15,0.2),col='red')
#Draw a density function curve
lines(x2,0.75*exp(-x2*x2/2)/sqrt(2*pi)+0.25*exp(-(x2-3)*(x2-3)/2)/sqrt(2*pi))

## ------------------------------------------------------------------------
#A function of drawing picture with p1
wf2<-function(n,p1){
y1 <- rnorm(n,mean=0)
y2 <- rnorm(n,mean=3)
#Assign value to p1
r <- sample(c(0,1),n,replace=TRUE,prob = c(p1,1-p1))
z <- r*y1+(1-r)*y2}
#Prepare to draw nine pictures
par(mar=c(1,1,1,1))
par(mfrow=c(3,3))
#A loop to draw picture
for(i in seq(2,10,1)){
p2=i/20
z<-wf2(10000,p2)
x3=seq(-5,10,0.2)
#Draw a histogram
hist(z,prob=TRUE,breaks = seq(-5,10,0.3),main = p2,col='blue') 
#Density function curve
lines(x3,(1-p2)*exp(-x3*x3/2)/sqrt(2*pi)+p2*exp(-(x3-3)*(x3-3)/2)/sqrt(2*pi))}

## ------------------------------------------------------------------------
wish<-function(n,sigma1){
#Cholesky algorithm
chol1<-chol(sigma1)
#Get dimension of chol1
dim1=dim(chol1)
A=matrix(nrow=dim1[1],ncol=dim1[2])
# produce matrix A
for(i in 1:dim1[1]){
  for (j in 1:dim1[2]){
    if (i<j){A[i,j]=0}
    else if(i>j){A[i,j]=rnorm(1)}
    else{A[i,j]=sqrt(rchisq(1,n-i+1,ncp=0))}
  }
}
X=t(chol1)%*%A%*%t(A)%*%chol1}

## ------------------------------------------------------------------------
sigma1<-matrix(c(5,1,1,5),ncol = 2)
n=10
wi<-wish(n,sigma1)
wi

## ----echo=FALSE----------------------------------------------------------
wi

## ------------------------------------------------------------------------
n<-1e4;
#generate 10000  random numbers between 0 and pi/3
y<-runif(n,min=0,max=pi/3);
#take average
a1<-mean(sin(y)*pi/3);
print(c(a1,0.5))

## ----echo=FALSE----------------------------------------------------------
n<-1e6;
#generate 10000  random numbers between 0 and pi/3
y<-runif(n,min=0,max=pi/3);
#take average
a1<-mean(sin(y)*pi/3);
print(c(a1,0.5))

## ------------------------------------------------------------------------
yu<-function(m){
z<-runif(m,min=0,max = 1);
z1<-exp(-z)/(1+z^2)*1
z2<-mean(z1)
z2}

## ------------------------------------------------------------------------
anti<-function(m){
m1<-m/2
x<-runif(m1,min=0,max = 1);
x1<-exp(-x)/(1+x^2)*1
x2<-exp(-(1-x))/(1+(1-x)^2)*1
x3<-c(x1,x2)
x4<-mean(x3)
}

## ------------------------------------------------------------------------
iter=1000
aa=numeric(iter)
bb=numeric(iter)
for (i in 1:iter) {
  aa[i]=yu(10000)
  bb[i]=anti(10000)
}
print(c(sd(bb),sd(aa),sd(bb)/sd(aa)))

## ------------------------------------------------------------------------
aa1=aa[1:100]
bb1=bb[1:100]
plot(aa1,type='l',xlab='aa1 and bb1',ylab='',col='red',main='make a comparision')
lines(bb1,col='blue')


## ------------------------------------------------------------------------
m2<-100000;
k<-5
N<-50#number of iteration
sj<-numeric(m2/k)
T<-numeric(k)
mat1=numeric(N)
sd=numeric(k)
ssd=numeric(N)
g<-function(x)exp(-x)/(1+x^2)*(x>0)*(x<1)
  #{exp(-x-log(1+x^2))*(x>0)*(x<1)}
f<-function(x)
{-log(1-x*(1-exp(-1)))}
duan=c(0,0.2,0.4,0.6,0.8,1)
for(i in 1:N){
  for(j in 1:5){
    sj<-f(runif(m2/k,duan[j],duan[j+1]))#generate x
    T[j]<-mean(5*g(sj)/(5*exp(-sj)/(1-exp(-1))))#每段均值
    sd<-sd(g(sj)/(exp(-sj)/(1-exp(-1))))
    }
  mat1[i]<-mean(T)
  ssd[i]<-mean(sd)
}  
mean(mat1)
mean(ssd)

## ----echo=FALSE----------------------------------------------------------
m2<-100000;
k<-1
N<-50#number of iteration
sj<-numeric(m2/k)
T<-numeric(k)
mat1=numeric(N)
sd=numeric(k)
ssd=numeric(N)
g<-function(x)exp(-x)/(1+x^2)*(x>0)*(x<1)
  #{exp(-x-log(1+x^2))*(x>0)*(x<1)}
f<-function(x)
{-log(1-x*(1-exp(-1)))}
duan=c(0,1)
for(i in 1:N){
  for(j in 1:1){
    sj<-f(runif(m2/k,duan[j],duan[j+1]))#generate x
    T[j]<-mean(5*g(sj)/(5*exp(-sj)/(1-exp(-1))))#每段均值
    sd<-sd(g(sj)/(exp(-sj)/(1-exp(-1))))
    }
  mat1[i]<-mean(T)
  ssd[i]<-mean(sd)
}  
mean(mat1)
mean(ssd)

## ----echo=FALSE----------------------------------------------------------
m2<-100000;
k<-2
N<-50#number of iteration
sj<-numeric(m2/k)
T<-numeric(k)
mat1=numeric(N)
sd=numeric(k)
ssd=numeric(N)
g<-function(x)exp(-x)/(1+x^2)*(x>0)*(x<1)
  #{exp(-x-log(1+x^2))*(x>0)*(x<1)}
f<-function(x)
{-log(1-x*(1-exp(-1)))}
duan=c(0,0.5,1)
for(i in 1:N){
  for(j in 1:2){
    sj<-f(runif(m2/k,duan[j],duan[j+1]))#generate x
    T[j]<-mean(5*g(sj)/(5*exp(-sj)/(1-exp(-1))))#每段均值
    sd<-sd(g(sj)/(exp(-sj)/(1-exp(-1))))
    }
  mat1[i]<-mean(T)
  ssd[i]<-mean(sd)
}  
mean(mat1)
mean(ssd)

## ----echo=FALSE----------------------------------------------------------
m2<-100000;
k<-10
N<-50#number of iteration
sj<-numeric(m2/k)
T<-numeric(k)
mat1=numeric(N)
sd=numeric(k)
ssd=numeric(N)
g<-function(x)exp(-x)/(1+x^2)*(x>0)*(x<1)
  #{exp(-x-log(1+x^2))*(x>0)*(x<1)}
f<-function(x)
{-log(1-x*(1-exp(-1)))}
duan=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
for(i in 1:N){
  for(j in 1:10){
    sj<-f(runif(m2/k,duan[j],duan[j+1]))#generate x
    T[j]<-mean(5*g(sj)/(5*exp(-sj)/(1-exp(-1))))#每段均值
    sd<-sd(g(sj)/(exp(-sj)/(1-exp(-1))))
    }
  mat1[i]<-mean(T)
  ssd[i]<-mean(sd)
}
mean(mat1)
mean(ssd)


## ------------------------------------------------------------------------
k.value<-c(1,2,5,10)
mean.value<-c(0.524,0.524,0.524,0.524)
sd.value<-c(0.0968,0.0692,0.0307,0.0156)
knitr::kable(cbind(k.value,mean.value,sd.value))

## ------------------------------------------------------------------------
wf1<-function(n,alpha){
x <- rchisq(n,2)
UCL11 <- mean(x)+sd(x) *qt(alpha/2,n-1)/ sqrt(n)
UCL12<-mean(x)+ sd(x) *qt(1-alpha/2,n-1)/ sqrt(n)
c(UCL11,UCL12)
}

## ------------------------------------------------------------------------
n=10000
data1=matrix(0,nrow = 2,ncol = n)
for(i in 1:n){
  data1[,i]=wf1(20,0.05)
}
length(which(data1[1,]<2&data1[2,]>2))/n
new.t1<-mean(data1[1,])
new.t2<-mean(data1[2,])
new.t1
new.t2

## ------------------------------------------------------------------------
cv<-function(alpha){
n <- c(10, 20, 30, 50, 100, 500)#sizes of samples
cv <- qnorm(alpha, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))#crit. values for each n
round(cv,5)#Five decimal places
}

## ------------------------------------------------------------------------
#write a function to compute the sample skewness coeff.
sk <- function(x) {
xbar <- mean(x)
m3 <- mean((x - xbar)^3)#molecule
m2 <- mean((x - xbar)^2)#for calculating Denominator
return( m3 / m2^1.5 )#return skewness
}


## ------------------------------------------------------------------------
n <- c(10, 20, 30, 50, 100, 500)
var<-p.reject <- matrix(0,4,length(n)) #to store sim. results
m <- 100000 #num. repl. each sim.
num<-0
for (alpha in c(0.9875,0.975,0.525,0.5125)){
    num<-num+1
    cv1<-cv(alpha)
    for (i in 1:length(n)) {
        sktests <- numeric(m) #test decisions
        for (j in 1:m) {
            x <- rnorm(n[i])
            #test decision is 1 (reject) or 0
            sktests[j] <- sum(abs(sk(x)) >= cv1[i] )
         }
         p.reject[num,i] <- mean(sktests) #proportion  rejected
         sigma<-sqrt(6*(n[i]-2) / ((n[i]+1)*(n[i]+3)))
         var[num,i]<-alpha*(1-alpha)/(n[i]*(exp(-cv1[i]^2/2/sigma^2)/(sigma*sqrt(2*pi)))^2)
    }
}
p.reject
var

## ------------------------------------------------------------------------
par(mar=c(1,1,1,1))
par(mfrow=c(1,2))
break1 = seq(0,1,0.02)

hist(rbeta(3000,10,10),prob=TRUE,break1,col = 'yellow')
lines(break1,dbeta(break1,10,10))

hist(rbeta(3000,2,6),prob=TRUE,breaks = seq(0,1,0.02),col = 'yellow')
lines(break1,dbeta(break1,2,6))

## ------------------------------------------------------------------------
sk <- function(x) {
  #computes the sample skewness coeff.
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}

## ------------------------------------------------------------------------
alpha <- 0.05
n <- 30
m <- 1000
seq1 <-seq(2, 10, 0.2)
N <- length(seq1)
pwr <- numeric(N) #Store power value
#critical value for the skewness test
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { 
  sktests <- numeric(m)
  for (i in 1:m) { #for each replicate
    x <- rbeta(30, 2, seq1[j])
    sktests[i] <- as.integer(abs(sk(x)) >= cv)
  }
  pwr[j] <- mean(sktests)
}
#plot power vs seq1
plot(seq1, pwr, type = "b",
     xlab = bquote(seq1), ylim = c(0,1),col='blue')
abline(h = 0.05, lty = 4) #Increase the confidence level of 0.05
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(seq1, pwr+se, lty = 3,col='violet')
lines(seq1, pwr-se, lty = 3,col='violet')

## ------------------------------------------------------------------------
alpha <- 0.05
n <- 30
m <- 1000
seq1 <-seq(2, 10, .2)
N <- length(seq1)
pwr <- numeric(N)
#critical value for the skewness test
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { 
  sktests <- numeric(m)
  for (i in 1:m) { #for each replicate
    x <- rbeta(30, seq1[j], seq1[j])
    sktests[i] <- as.integer(abs(sk(x)) >= cv)
  }
  pwr[j] <- mean(sktests)
}
#plot power vs epsilon
plot(seq1, pwr, type = "b",
     xlab = bquote(seq1), ylim = c(0,0.2),col='blue')
abline(h = .05, lty = 4) #Increase the confidence level of 0.05
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(seq1, pwr+se, lty = 3,col='violet')
lines(seq1, pwr-se, lty = 3,col='violet')

## ------------------------------------------------------------------------
alpha <- 0.05
n <- 30
m <- 1000
seq2 <- c(seq(0, .15, .01), seq(.15, 1, .05))
N <- length(seq2)
pwr <- numeric(N)
#critical value for the skewness test
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { #for each epsilon
  e <- seq2[j]
  sktests <- numeric(m)
  for (i in 1:m) { #for each replicate
    par <- sample(c(2,10), replace = TRUE,
                    size = n, prob = c(1-e, e))
    x <- rbeta(30, par, par)
    sktests[i] <- as.integer(abs(sk(x)) >= cv)
  }
  pwr[j] <- mean(sktests)
}
#plot power vs epsilon
plot(seq2, pwr, type = "b",
     xlab = bquote(seq2),col='blue')
abline(h = 0.05, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(seq2, pwr+se, lty = 3,col='violet')
lines(seq2, pwr-se, lty = 3,col='violet')


## ------------------------------------------------------------------------
alpha <- .05
n <- 50
m <- 1000
seq1 <-seq(2, 50, 1)
N <- length(seq1)
pwr <- numeric(N)
#critical value for the skewness test
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { 
  #print(j)
  sktests <- numeric(m)
  for (i in 1:m) { #for each replicate
    x <- rt(n,seq1[j])
    sktests[i] <- as.integer(abs(sk(x)) >= cv)
  }
  pwr[j] <- mean(sktests)
}
#plot power vs epsilon
plot(seq1, pwr, type = "b",
     xlab = bquote(seq1),ylim = c(0,1),col='blue')
abline(h = .05, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(seq1, pwr+se, lty = 3,col='violet')
lines(seq1, pwr-se, lty = 3,col='violet')

## ------------------------------------------------------------------------
alpha <- 0.05
n <- 30
m <- 1000
seq2 <- c(seq(0, .15, .01), seq(.15, 1, .05))
N <- length(seq2)
pwr <- numeric(N)
#critical value for the skewness test
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { #for each epsilon
  e <- seq2[j]
  sktests <- numeric(m)
  for (i in 1:m) { #for each replicate
    par <- sample(c(2,80), replace = TRUE,
                    size = n, prob = c(1-e, e))
    x <- rt(30, par)
    sktests[i] <- as.integer(abs(sk(x)) >= cv)
  }
  pwr[j] <- mean(sktests)
}
#plot power vs epsilon
plot(seq2, pwr, type = "b",
     xlab = bquote(seq2),col='green')
abline(h = 0.05, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(seq2, pwr+se, lty = 3,col='yellow')
lines(seq2, pwr-se, lty = 3,col='yellow')


## ------------------------------------------------------------------------
n <- 1000
mu0 <- 1
m <- 10000 #number of replicates
p <- numeric(m) #storage for p-values
chisq.hat<- numeric(3)
fw<-c(0.025,0.05,0.1)
for(i in 1:3){
   for (j in 1:m) {
      x <- rchisq(n, 1)
      ttest <- t.test(x, mu = mu0)
      p[j] <- ttest$p.value
   }
 chisq.hat[i]<- mean(p < fw[i])
}  
se.hat1 <- sqrt(chisq.hat * (1 - chisq.hat) / m)
print(c(chisq.hat, se.hat1))

## ------------------------------------------------------------------------
n <- 1000
alpha <- .05
mu0 <- 1
m <- 10000 #number of replicates
p <- numeric(m) #storage for p-values
unif.hat<- numeric(3)
fw<-c(0.025,0.05,0.1)
for(i in 1:3){
    for (j in 1:m) {
         x <- runif(n, 0,2)
         ttest <- t.test(x, mu = mu0)
         p[j] <- ttest$p.value
    }
  unif.hat[i]<- mean(p < fw[i])
}
se.hat2 <- sqrt(unif.hat * (1 - unif.hat) / m)
print(c(unif.hat, se.hat2))

## ------------------------------------------------------------------------
n <- 1000
alpha <- .05
mu0 <- 1
m <- 10000 #number of replicates
p <- numeric(m) #storage for p-values
exp.hat<- numeric(3)
fw<-c(0.025,0.05,0.1)
for(i in 1:3){
    for (j in 1:m) {
        x <- rexp(n, 1)
        ttest <- t.test(x, mu = mu0)
         p[j] <- ttest$p.value
    }
    exp.hat[i] <- mean(p < fw[i])
}

se.hat3 <- sqrt(exp.hat * (1 - exp.hat) / m)
print(c(exp.hat, se.hat3))

## ------------------------------------------------------------------------
data.frame(fw,chisq.hat,unif.hat,exp.hat)


## ------------------------------------------------------------------------
library(bootstrap) #for the score data
cor(scor)
#show a part of data
scor[1:5,]

## ------------------------------------------------------------------------
par(mar=c(1,1,1,1))
par(mfrow=c(2,5))
for(i in 1:5){
  for(j in 1:5){
    if(i!=j & i<j){
    plot(scor[,i],scor[,j],col='red')}
  }
}

## ------------------------------------------------------------------------
N=300 # Initialize the numbers of replicates
hat.rho12<-numeric(N)
n<-nrow(scor)
for(i in 1:N){
  samp<-sample(1:n,n,replace = TRUE)
  hat.rho12[i]<-cor(scor$mec[samp],scor$vec[samp])
}
# Get Standard deviation of Correlation coefficient
hat.sd12<-mean(hat.rho12)
print(hat.sd12)
# real Correlation coefficient
real.rho12<-cor(scor$mec,scor$vec)
print(real.rho12)
#plot and compare with real Correlation coefficient
plot(hat.rho12,ylim = c(0,1),col='green')
abline(h = real.rho12, lty = 4,col='red')
#Draw a histogram and a true correlation coefficient vertical line
hist(hat.rho12)
abline(v=real.rho12,col='red',lwd=2)

## ------------------------------------------------------------------------
hat.rho34<-numeric(N)
for(i in 1:N){
  samp<-sample(1:n,n,replace = TRUE)
  hat.rho34[i]<-cor(scor$alg[samp],scor$ana[samp])
}
# Get Standard deviation of Correlation coefficient
hat.sd34<-mean(hat.rho34)
print(hat.sd34)
# real Correlation coefficient
real.rho34<-cor(scor$alg,scor$ana)
print(real.rho34)
#plot and compare with real Correlation coefficient
plot(hat.rho34,ylim = c(0,1),col='green')
abline(h = real.rho34, lty = 4,col='red')
#Draw a histogram and a true correlation coefficient vertical line
hist(hat.rho34)
abline(v=real.rho34,col='red',lwd=2)

## ------------------------------------------------------------------------
hat.rho35<-numeric(N)
for(i in 1:N){
  samp<-sample(1:n,n,replace = TRUE)
  hat.rho35[i]<-cor(scor$alg[samp],scor$sta[samp])
}
# Get Standard deviation of Correlation coefficient
hat.sd35<-mean(hat.rho35)
print(hat.sd35)
# real Correlation coefficient
real.rho35<-cor(scor$alg,scor$sta)
print(real.rho35)
#plot and compare with real Correlation coefficient
plot(hat.rho35,ylim = c(0,1),col='green')
abline(h = real.rho35, lty = 4,col='red')
#Draw a histogram and a true correlation coefficient vertical line
hist(hat.rho35)
abline(v=real.rho35,col='red',lwd=2)

## ------------------------------------------------------------------------
hat.rho45<-numeric(N)
for(i in 1:N){
  samp<-sample(1:n,n,replace = TRUE)
  hat.rho45[i]<-cor(scor$ana[samp],scor$sta[samp])
}
# Get Standard deviation of Correlation coefficient
hat.sd45<-mean(hat.rho45)
print(hat.sd45)
# real Correlation coefficient
real.rho45<-cor(scor$ana,scor$sta)
print(real.rho45)
#plot and compare with real Correlation coefficient
plot(hat.rho45,ylim = c(0,1),col='green')
abline(h = real.rho45, lty = 4,col='red')
#Draw a histogram and a true correlation coefficient vertical line
hist(hat.rho45)
abline(v=real.rho45,col='red',lwd=2)

## ------------------------------------------------------------------------
sk <- function(x,j) {
  #computes the sample skewness coeff.
  xbar <- mean(x[j])
  m3 <- mean((x[j] - xbar)^3)
  m2 <- mean((x[j] - xbar)^2)
  return( m3 / m2^1.5 )
}

## ------------------------------------------------------------------------
me<-0
n<-20
m<-1000
set.seed(12345)
library(boot)
nor.norm<-nor.basic<-nor.perc<-matrix(NA,m,2)
for (i in 1:m) {
  data.nor<-rnorm(n,0,2)
  nor.ske<-boot(data.nor,statistic=sk,R=1000)
  nor<- boot.ci(nor.ske,type=c("norm","basic","perc"))
  nor.norm[i,]<-nor$norm[2:3];
  nor.basic[i,]<-nor$basic[4:5];
  nor.perc[i,]<-nor$percent[4:5];
}
#Calculate the coverage probability of a normal distribution
cat('norm =',mean(nor.norm[,1]<=me & nor.norm[,2]>=me),
    'basic =',mean(nor.basic[,1]<=me & nor.basic[,2]>=me),
    'perc =',mean(nor.perc[,1]<=me & nor.perc[,2]>=me))
#Calculate the probability of the left side of the normal distribution
cat('norm.left=',mean(nor.norm[,1]>=me ),
    'basic.left =',mean(nor.basic[,1]>=me ),
    'perc.left =',mean(nor.perc[,1]>=me))
#Calculate the right side probability of a normal distribution
cat('norm.right=',mean(nor.norm[,2]<=me ),
    'basic.right =',mean(nor.basic[,2]<=me ),
    'perc.right =',mean(nor.perc[,2]<=me))

## ------------------------------------------------------------------------
pian.real<-mean(replicate(1000,expr = {
  x0<-rchisq(1000,5)
  mean((x0-5)^3)/10^1.5
}))
print(pian.real)

## ------------------------------------------------------------------------
me<-pian.real
n<-20
m<-1000
set.seed(12345)
library(boot)
chi.norm<-chi.basic<-chi.perc<-matrix(NA,m,2)
for (i in 1:m) {
  data.chisq<-rchisq(n,5)
  chisq.ske<-boot(data.chisq,statistic=sk,R=1000)
  chi<- boot.ci(chisq.ske,type=c("norm","basic","perc"))
  chi.norm[i,]<-chi$norm[2:3];
  chi.basic[i,]<-chi$basic[4:5];
  chi.perc[i,]<-chi$percent[4:5];
}
#Calculate the coverage probability of the chi-square distribution
cat('norm =',mean(chi.norm[,1]<=me & chi.norm[,2]>=me),
    'basic =',mean(chi.basic[,1]<=me & chi.basic[,2]>=me),
    'perc =',mean(chi.perc[,1]<=me & chi.perc[,2]>=me))
#Calculate the probability of the left side of the chi-square distribution
cat('norm.left=',mean(chi.norm[,1]>=me ),
    'basic.left =',mean(chi.basic[,1]>=me ),
    'perc.left =',mean(chi.perc[,1]>=me))
#Calculate the right side probability of the chi-square distribution
cat('norm.right=',mean(chi.norm[,2]<=me ),
    'basic.right =',mean(chi.basic[,2]<=me ),
    'perc.right =',mean(chi.perc[,2]<=me))

## ------------------------------------------------------------------------
library(bootstrap) #for the score data
cov(scor)
#show a part of data
scor[1:5,]

## ------------------------------------------------------------------------
mcov <- function(x,j) cov(x[j,])

## ------------------------------------------------------------------------
#Computing data dimension
n<-lengths(scor)[1]
#Calculate the real covariance matrix
real.cov <- mcov(scor,1:n)
#Calculate the real theta value
the.hat<-eigen(real.cov)$values[1]/sum(eigen(real.cov)$values)
#Define variables to store values
the.jack <- numeric(n)
#Loop n times
for(i in 1:n){
the.cov <- mcov(scor,(1:n)[-i])
#Get matrix eigenvalues
eig1<-eigen(the.cov)$values
#Find the estimate of lambda
the.jack[i]<-eig1[1]/sum(eig1)
}
#Calculation bias
bias1 <- (n-1)*(mean(the.jack)-the.hat)
#Calculating  standard deviation
se1 <- sqrt((n-1)*mean((the.jack-the.hat)^2))
round(c(original=the.hat,estimate=mean(the.jack),bias1=bias1,
se1=se1),3)

## ------------------------------------------------------------------------
r2<-function(x,y,p){
  a<-(n-2)*(sum((y-mean(x))^2)/sum((x-mean(x))^2))/(n-2-p)
}

## ------------------------------------------------------------------------
library(DAAG);
attach(ironslag);

## ------------------------------------------------------------------------
n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n)
# for n-fold cross validation
# fit models on leave-one-out samples
yhat1<-yhat2<-yhat3<-yhat4<-numeric(n)
for (k in 1:n) {
y <- magnetic[-k]
x <- chemical[-k]
#Linear model
J1 <- lm(y ~ x)
yhat1[k] <- J1$coef[1] + J1$coef[2] * chemical[k]
e1[k] <- magnetic[k] - yhat1[k]

#Quadratic model
J2 <- lm(y ~ x + I(x^2))
yhat2[k] <- J2$coef[1] + J2$coef[2] * chemical[k] +
J2$coef[3] * chemical[k]^2
e2[k] <- magnetic[k] - yhat2[k]

#Exponential model
J3 <- lm(log(y) ~ x)
logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
yhat3[k] <- exp(logyhat3)
e3[k] <- magnetic[k] - yhat3[k]

# cubic polynomial model
J4 <- lm(y ~ x+I(x^2)+I(x^3))
yhat4[k]<- J4$coef[1] + J4$coef[2]* chemical[k] +J4$coef[3]* chemical[k]^2+J4$coef[4]* chemical[k]^3
e4[k]<- magnetic[k] - yhat4[k]

}
#calculate R^2 of four methods
rsq1<-r2( magnetic,yhat1,1)
rsq2<-r2( magnetic,yhat2,2)
rsq3<-r2( magnetic,yhat3,1)
rsq4<-r2( magnetic,yhat4,3)

## ------------------------------------------------------------------------
c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))
c(rsq1,rsq2,rsq3,rsq4)

## ------------------------------------------------------------------------
rm(list = ls())
count5test <- function(z,ix,sizes) {
n11<-sizes[1]
n22<-sizes[2]
n0<-n11+n22
z<-z[ix]
x1<-z[1:n11]
y1<-z[(n11+1):n0]
X <- x1 - mean(x1)
Y <- y1- mean(y1)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0)
return(max(c(outx, outy)))
}

## ------------------------------------------------------------------------
library(boot)
n1<-30
n2<-c(30,40,50,60,70,80)
u1<-u2<-0
sig1<-1
sig2<-1
p.value<-numeric(6)
a=1
for (i in n2) {
   n3<-i
   size1<-c(n1,n3)
   x<-rnorm(n1,u1,sig1)
   y<-rnorm(n3,u2,sig2)
   z<-c(x,y)
   set.seed(4*i);
   boot.obj <- boot(data = z, statistic = count5test, R = 999, sim = "permutation",sizes=size1)
   ts <- c(boot.obj$t0,boot.obj$t)
   p.value[a]<- mean(ts>=ts[1])
   a=a+1
}
data.frame(n2,p.value)

## ------------------------------------------------------------------------
rm(list = ls())
count5test <- function(z,ix,sizes,kk) {
n11<-sizes[1]
n22<-sizes[2]
n0<-n11+n22
z<-z[ix]
x1<-z[1:n11]
y1<-z[(n11+1):n0]
X <- x1 - mean(x1)
Y <- y1- mean(y1)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0)
return(as.integer(max(c(outx, outy))>kk))
}

## ------------------------------------------------------------------------
library(boot)
n1<-30
n2<-c(35,40,50,60,70,80)
u1<-u2<-0
sig1<-1
sig2<-1
p.value2<-numeric(6)
a=1
for (i in n2) {
   n3<-i
   size1<-c(n1,n3)
   x<-rnorm(n1,u1,sig1)
   y<-rnorm(n3,u2,sig2)
   z<-c(x,y)
   set.seed(i);
   boot.obj <- boot(data = z, statistic = count5test, R = 999, sim = "permutation",sizes=size1,kk=a+4)
   ts <- c(boot.obj$t0,boot.obj$t)
   p.value2[a] <- mean(ts!=ts[1])
   a<-a+1
}
data.frame(n2,p.value2)
plot(n2,p.value2,type = 'o',col='red',ylim = c(0,1))

## ------------------------------------------------------------------------
dCov <- function(x, y) {
    x <- as.matrix(x); y <- as.matrix(y)
    n <- nrow(x); m <- nrow(y)
    if (n != m || n < 2) stop("Sample sizes must agree")
    if (! (all(is.finite(c(x, y)))))
    stop("Data contains missing or infinite values")
    Akl <- function(x) {
        d <- as.matrix(dist(x))
        m <- rowMeans(d); M <- mean(d)
        a <- sweep(d, 1, m); b <- sweep(a, 2, m)
        b + M
    }
    A <- Akl(x); B <- Akl(y)
    sqrt(mean(A * B))
}


## ----cars----------------------------------------------------------------
library(MASS)
library(Ball)
#options(digits=3)
alpha=0.1
mean <- c(0,0)
sigma <- matrix(c(1,0,0,1),nrow=2,ncol=2)
cn<-c(10,20,30,40,50,60,80,100)
p.value1 <- matrix(NA,8,4)
a=1

ndCov2 <- function(z, ix, dims) {
    #dims contains dimensions of x and y
    p <- dims[1]
    q <- dims[2]
    d <- p + q
    x0 <- z[, 1:p] #leave x as is
    y0 <- z[ix, -(1:p)] #permute rows of y
    return(nrow(z) * dCov(x0, y0)^2)
}

for (i in cn) { 
    p.cor1<-numeric(50)
    p.cor2<-numeric(50)
    p.cor3<-numeric(50)
    p.cor4<-numeric(50)
  for (j in 1:50) {
      
    xn<-mvrnorm(n=i,mean,sigma)
    en<-mvrnorm(n=i,mean,sigma)
    yn1<-xn/4+en
    yn2<-xn/4*en
    set.seed(12345)
    
    z1 <- cbind(xn,yn1)
    boot.obj1 <- boot(data = z1, statistic = ndCov2, R = 199,sim = "permutation", dims = c(2, 2))
    tb1 <- c(boot.obj1$t0, boot.obj1$t)
    p.cor1[j]<- mean(tb1>=tb1[1])
    p.cor3[j]<-bcov.test(x = xn, y = yn1, R=999,seed =i*1234*j)$p.value
    z2 <- cbind(xn,yn2)
    boot.obj2<- boot(data = z2, statistic = ndCov2, R = 199,sim = "permutation", dims = c(2, 2))
    tb2 <- c(boot.obj2$t0, boot.obj2$t)
    p.cor2[j] <- mean(tb2>=tb2[1])
    p.cor4[j] <-bcov.test(x = xn,y =yn2,R=999,seed = i*1233*j)$p.value
  }
  p.value1[a,1]<-mean(p.cor1<alpha)
  p.value1[a,2]<-mean(p.cor3<alpha)
  p.value1[a,3]<-mean(p.cor2<alpha)
  p.value1[a,4]<-mean(p.cor4<alpha)
  a=a+1
  print(a)
}
par(mar=c(1,1,1,1))
par(mfrow=c(1,2))
    plot(cn,p.value1[,1],type='l',col='red',ylim = c(0,1))
    lines(cn,p.value1[,2],type='l',col='yellow')
    plot(cn,p.value1[,3],type='l',col='red',ylim = c(0,1))
    lines(cn,p.value1[,4],type='l',col='yellow')

## ------------------------------------------------------------------------
rp.Metropolis <- function(sigma, x0, N) {
x <- numeric(N)
x[1] <- x0
u <- runif(N)
k <- 0
for (i in 2:N) {
y <- rnorm(1, x[i-1], sigma)
#Replace probability density function with Laplace distribution
if (u[i] <= (0.5*exp(-abs(y))) / (0.5*exp(-abs(x[i-1]))))
    x[i] <- y 
else {
    x[i] <- x[i-1]
    k <- k + 1
}
}
return(list(x=x, k=k))
}

## ------------------------------------------------------------------------
N <- 2000
sigma <- c(.2, .5, 2, 16)
x0 <- 25
#Get the return value of different variances
rp1 <- rp.Metropolis(sigma[1], x0, N)
rp2 <- rp.Metropolis(sigma[2], x0, N)
rp3 <- rp.Metropolis(sigma[3], x0, N)
rp4 <- rp.Metropolis(sigma[4], x0, N)
#Calculate acceptance probability
p1<-1-rp1$k/N
p2<-1-rp2$k/N
p3<-1-rp3$k/N
p4<-1-rp4$k/N
data.frame(p1,p2,p3,p4)

## ------------------------------------------------------------------------
   
    par(mar=c(1,1,1,1))
    par(mfrow=c(2,2))  #display 4 graphs together
    refline <- c(log(0.05),-log(0.05))
    rp <- cbind(rp1$x, rp2$x, rp3$x,  rp4$x)
    for (j in 1:4) {
        plot(rp[,j], type="l",
             xlab=bquote(sigma == .(round(sigma[j],3))),
             ylab="X", ylim=c(-5,30),col='blue')
        abline(h=refline,col='red')
    }
    

## ------------------------------------------------------------------------
seq1<-seq(1,100,1)
tr1<-tr2<-numeric(length(seq1))
wf<-function(i){
return(exp(log(i))==log(exp(i)))}
dx<-function(i){
  isTRUE(all.equal(exp(log(i)),log(exp(i))))
}
for (i in seq1){
  tr1[i]<-wf(i)
  tr2[i]<-dx(i)
}
mean(tr1)
mean(tr2)

## ------------------------------------------------------------------------
rm(list=ls())
k=4
ck<-function(k,a){
  return(sqrt(a^2*k/(k+1-a^2)))
  }
wf1 <- function(u) {
(1+u*u/(k-1))^(-k/2)
}
wf2<- function(u) {
(1+u*u/k)^(-(k+1)/2)
}
wf3<-function(a){
  2*gamma(k/2)/(sqrt(pi*(k-1))*gamma((k-1)/2))*integrate(wf1,0,ck(k-1,a))$value-2*gamma((k+1)/2)/(sqrt(pi*k)*gamma(k/2))*integrate(wf2,0,ck(k,a))$value
}

## ------------------------------------------------------------------------
kt<-c(seq(4,25,1),100,200)
solution11.5<-solution11.4<-numeric(length(kt))
js=1
for (i in kt) {
  k<-kt[js]
  solution11.5[js]<-uniroot(wf3,c(0.001,sqrt(k)/2+1))$root
  js=js+1
}

## ------------------------------------------------------------------------
js=1
wf4<-function(a){pt(sqrt(a^2*(k-1)/(k-a^2)),k-1)-pt(sqrt(a^2*k/(k+1-a^2)),k)}
for (i in kt) {
   k<-kt[js]
  solution11.4[js]<-uniroot(wf4,c(0.00001,sqrt(k)-0.0001))$root
  js=js+1
}
data.frame(kt,solution11.4,solution11.5)

## ------------------------------------------------------------------------
formulas1<- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)
j=1:4
li=list()
#Use lapply
la<-lapply(j,function(j)lm(formulas1[[j]],data = mtcars))
la

#Or use a for loop
for (j in 1:4) {
  li[[j]]<-lm(formulas1[[j]],data = mtcars)
}
li

## ------------------------------------------------------------------------
rm(list = ls())
formulas <- list(
l1<-function(mtcars) glm(mtcars$mpg ~ mtcars$disp),l2<-function(mtcars) glm(mtcars$mpg ~ I(1 / mtcars$disp)),l3<-function(mtcars) glm(mtcars$mpg ~ mtcars$disp + mtcars$wt),l4<-function(mtcars) glm(mtcars$mpg ~ I(1 / mtcars$disp) +mtcars$wt)
)
#Use lapply
la1<-lapply(formulas, function(f) f(mtcars))

#Or use a for loop
li1=list()
for (i in 1:4){
  #li1<-list(li1,formulas[[i]](mtcars))
  li1[[i]]<-formulas[[i]](mtcars)
}

## ------------------------------------------------------------------------
formulas1<- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)
  bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
  })
li0=list()
#Use lapply
j=1:10
la0<-lapply(j,function(j)lm(formulas1[[1]],data = bootstraps[[j]]))
la0
#Or use a for loop
for (j in 1:10) {
  li0[[j]]<-lm(formulas1[[1]],data = bootstraps[[j]])
}
li0

## ------------------------------------------------------------------------
#Use lapply
la2<-lapply(bootstraps,formulas[[1]])
li2<-list()
#Or use a for loop
for (i in 1:10) {
  li2[[i]]<-formulas[[1]](bootstraps[[i]])
}

## ------------------------------------------------------------------------
formulas <- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)
wf <- numeric(4)
rsq <- function(mod) summary(mod)$r.squared
j = 1:4
wf <- lapply(j,function(j) {rsq(lm(formulas[[j]], data = mtcars))})
wf

## ------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
  })
wf <- numeric(10)
rsq <- function(mod) summary(mod)$r.squared
j = 1:10
wf <- lapply(j,function(j) {rsq(lm(formulas[[1]], data = bootstraps[[j]]))})
wf


## ------------------------------------------------------------------------
trials2 <- replicate(
100,
t.test(rpois(10, 10), rpois(7, 10)),
simplify = FALSE
)
tq<-function(a)a$p.value
sapply(trials2,tq)

## ------------------------------------------------------------------------
sapply(trials2,"[[",3)

## ------------------------------------------------------------------------
library(parallel)
formulas <- list(
l1<-function(mtcars) glm(mtcars$mpg ~ mtcars$disp),l2<-function(mtcars) glm(mtcars$mpg ~ I(1 / mtcars$disp)),l3<-function(mtcars) glm(mtcars$mpg ~ mtcars$disp + mtcars$wt),l4<-function(mtcars) glm(mtcars$mpg ~ I(1 / mtcars$disp) +mtcars$wt)
)
cl<-makeCluster(4)
bootstraps2 <- lapply(1:3000, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]})
system.time(sapply(bootstraps2,formulas[[1]]))
system.time(parSapply(cl,bootstraps2,formulas[[1]]))

## ------------------------------------------------------------------------
rp.Metropolis <- function(sigma, x0, N) {
x <- numeric(N)
x[1] <- x0
u <- runif(N)
k <- 0
for (i in 2:N) {
   y <- rnorm(1, x[i-1], sigma)
   #Replace probability density function with Laplace distribution
   if (u[i] <= (0.5*exp(-abs(y))) / (0.5*exp(-abs(x[i-1]))))
       x[i] <- y 
   else {
       x[i] <- x[i-1]
       k <- k + 1
        }
     }
return(list(x=x, k=k))
}

## ------------------------------------------------------------------------
library(Rcpp)
cppFunction('NumericVector rp(double sigma,double xx,int N){
#include<Rcpp.h>
NumericVector x(N+1);
x[0]=xx;
int k=0;
for(int i=1;i<N;++i){
  srand(i);
  double u=rand()/(RAND_MAX+1.0);
  double y=as<double>(rnorm(1, x[i-1], sigma));
    if(u<=(0.5*exp(-abs(y)))/(0.5*exp(-abs(x[i-1]))))
        {x[i]=y;}
    else{
        x[i]=x[i-1];
            k=k+1;
            }
}
x[N]=k;
return x;
}
            ')

## ------------------------------------------------------------------------
N <- 2000
sigma <- c(3, 5, 10, 16)
library(car)
x0 <-5
library(microbenchmark)

pic<-function(i){
ts <- microbenchmark(rp1 <- rp.Metropolis(sigma[i], x0, N),rpp1 <- rp(sigma[i], x0, N))
tsr<-summary(ts)[,c(1,3,5,6)]
#rpp1
par(mar=c(1,1,1,1))
par(mfrow=c(2,2))
#qqnorm(rp1$x,col='green')
#qqline(rp1$x,col='red')
#qqnorm(rpp1[1:N],col='green')
#qqline(rpp1[1:N],col='red')
qqplot(rp1$x,rpp1[1:N],col='green')
qqline(rpp1[1:N],col='red')
return(list(rpp1[N+1],tsr))
}

## ------------------------------------------------------------------------
par(mar=c(1,1,1,1))
par(mfrow=c(2,2))
for (i in 1:4) {
  res<-pic(i)
  print("Probability of rejection")
  print(res[1])
  print("Comparison of the two methods")
  print(res[2])
  
}

