library(mvtnorm)
library(MASS)
library(boot)
library(quantreg)
library(lpSolve)
library(extraDistr)

#tau=.5
p=4
n=20
set.seed(1)

#beta
beta<-rnorm(p+1,-1,1)

xMat<-rmvnorm(n,rep(10,p+1),diag(1,p+1))
xMat[,1]<-1

#Error distribution
error<-rcauchy(n=n,location = 0,scale = 5)#rt(n = 20,df = 1)
y<-xMat%*%beta+error

#LP formulation
a=c(rep(0,2*(p+1)),rep(1,n))
coMat<-matrix(c(1,-1,-1,1),2,2) #tackle the -ve coefficients
A=cbind(xMat%x%coMat,diag(1,n,n)%x%c(-1,-1))
b=y%x%c(1,-1)

#Deriving the LP
const_type <- rep("<=",2*n)
linprog <- lp("min",a,A,const_type,b)
linprog$solution[seq(1,2*(p+1),by=2)]-linprog$sol[seq(2,2*(p+1),by=2)]

#Quantreg Package
xData<-data.frame(xMat[,-1],y)
quantreg::rq(y~., tau=.5, data=xData)


#Function 
lpMedReg<-function(n,p,e_distr,err_p1=0,err_p2=5,seed=1){
  results=list()
  p=4
  n=20
  set.seed(seed)
  
  #beta
  beta<-rnorm(p+1,-1,1)
  
  #X matrix
  xMat<-rmvnorm(n,rep(10,p+1),diag(1,p+1))
  xMat[,1]<-1
  
  #Error distribution
  if (e_distr=='cauchy'){
    error<-rcauchy(n=n,location = err_p1,scale = err_p2)
  } else if(e_distr=='uniform'){
    error<-runif(n=n,min = err_p1,max= err_p2)
  } else if(e_distr=='laplace'){
    error<-rlaplace(n, mu = err_p1, sigma = err_p2)
  }else stop("Not compatible error function")
  
  #error<-rcauchy(n=n,location = 0,scale = 5)#rt(n = 20,df = 1)
  #y
  y<-xMat%*%beta+error
  
  #LP formulation
  a=c(rep(0,2*(p+1)),rep(1,n))
  coMat<-matrix(c(1,-1,-1,1),2,2)#tackle the -ve coefficients
  A=cbind(xMat%x%coMat,diag(1,n,n)%x%c(-1,-1))
  b=y%x%c(1,-1)
  
  #Deriving the LP
  const_type <- rep("<=",2*n)
  linprog <- lp("min",a,A,const_type,b)
  lpResult<-linprog$solution[seq(1,2*(p+1),by=2)]-linprog$sol[seq(2,2*(p+1),by=2)]
  
  #Quantreg Package
  xData<-data.frame(xMat[,-1],y)
  qResult<-quantreg::rq(y~., tau=.5, data=xData)
  qResult$coefficients
  
  #Results
  results[['trueValues']]=beta
  results[['quantregResults']]=qResult$coefficients
  results[['lpResults']]=lpResult
  
  return(results)
}
  
lpMedReg(n=20,p=4,e_distr='cauchy',err_p1=0,err_p2=5)
# $trueValues
# [1] -1.6264538 -0.8163567 -1.8356286  0.5952808 -0.6704922
# 
# $quantregResults
# (Intercept)          X1          X2          X3          X4 
# 2.582429   -2.980454   -2.712336    1.086279    1.514030 
# 
# $lpResults
# [1]  2.582429 -2.980454 -2.712336  1.086279  1.514030

lpMedReg(n=20,p=4,e_distr='laplace',err_p1=0,err_p2=5)
# $trueValues
# [1] -1.6264538 -0.8163567 -1.8356286  0.5952808 -0.6704922
# 
# $quantregResults
# (Intercept)          X1          X2          X3          X4 
# 21.6889349  -1.4884689  -4.2098585   1.5880037  -0.7657169 
# 
# $lpResults
# [1] 21.6889349 -1.4884689 -4.2098585  1.5880037 -0.7657169

lpMedReg(n=20,p=4,e_distr='uniform',err_p1=-10,err_p2=10)
# $trueValues
# [1] -1.6264538 -0.8163567 -1.8356286  0.5952808 -0.6704922
# 
# $quantregResults
# (Intercept)          X1          X2          X3          X4 
# -31.8703083   1.5128665  -0.4850726   2.6725677  -3.0509969 
# 
# $lpResults
# [1] -31.8703083   1.5128665  -0.4850726   2.6725677  -3.0509969