library(mvtnorm)
library(MASS)
library(boot)
library(quantreg)
library(lpSolve)

tau=.5
p=4
n=20
set.seed(1)
beta<-rnorm(5,-1,1)
xCov<-abs(rmvnorm(20,rep(10,4),diag(1,4)))
xMat<-cbind(1,xCov)
error<-rcauchy(n=20,location = 0,scale = 5)#rt(n = 20,df = 1)
y<-xMat%*%beta+error

a=c(rep(0,2*(p+1)),rep(1,n))
A=cbind((xMat%x%c(1,-1))%x%t(c(1,-1)),diag(1,n,n)%x%c(-1,-1))
b=y%x%c(1,-1)

simplex(a,A1=A,b1=b)
xData=data.frame(xCov,y)
quantreg::rq(y~., tau=.5, data=xData)

const_type <- rep("<=",2*n)
linprog <- lp("min",a,A,const_type,b)
linprog$solution[seq(1,2*(p+1),by=2)]-linprog$sol[seq(2,2*(p+1),by=2)]