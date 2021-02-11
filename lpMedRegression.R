library(MASS)
library(quantreg)
library(lpSolve)
library(tidyverse)


norm <- function(x) (x - min(x))/(max(x) - min(x))

p=1
n=10
set.seed(1)
beta<-c(sample(1:6,p+1,replace = T))

data("Pima.tr")

xMat<-Pima.tr %>% 
  mutate(cons=1,.before = everything())%>%
  slice(1:n) %>% 
  dplyr::select(cons,bp) %>%
  unlist()%>% 
  matrix(ncol = p+1, byrow = F)

set.seed(2)
error<-rt(n = n,df = 1)
y<-xMat%*%beta+error

xData<-Pima.tr %>% 
  slice(1:n) %>% 
  dplyr::select(bp) %>%
  #lapply(norm) %>% 
  data.frame()%>%
  add_column(y = y) 

a=c(rep(0,p+1),rep(1,n))
A=cbind(xMat%x%c(1,-1),diag(1,n,n)%x%rep(-1,2))
b=y%x%c(1,-1)

simplex(a,A1=A,b1=b)

quantreg::rq(y~., tau=.5, data=xData)

const_type <- rep("<=",2*n)
linprog <- lp("min",a,A,const_type,b)
linprog$solution[1:(p+1)]




##Quantile Regression
base=read.table("http://freakonometrics.free.fr/rent98_00.txt",header=TRUE)
attach(base)
tau <- 0.5
base<-base[1:30,]
# Problem (1) only one covariate
X <- cbind(1,base$area)
K <- ncol(X)
N <- nrow(X)

A <- cbind(X,-X,diag(N),-diag(N))
c <- c(rep(0,2*ncol(X)),tau*rep(1,N),(1-tau)*rep(1,N))
b <- base$rent_euro
const_type <- rep("=",N)

linprog <- lp("min",c,A,const_type,b)
beta <- linprog$sol[1:K] -  linprog$sol[(1:K+K)]
beta
rq(rent_euro~area, tau=tau, data=base)


# Problem (2) with 2 covariates
X <- cbind(1,base$area,base$yearc)
K <- ncol(X)
N <- nrow(X)

A <- cbind(X,-X,diag(N),-diag(N))
c <- c(rep(0,2*ncol(X)),tau*rep(1,N),(1-tau)*rep(1,N))
b <- base$rent_euro
const_type <- rep("=",N)

linprog <- lp("min",c,A,const_type,b)
beta <- linprog$sol[1:K] -  linprog$sol[(1:K+K)]
beta
rq(rent_euro~ area + yearc, tau=tau, data=base)



