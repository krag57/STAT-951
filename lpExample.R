library(lpSolve)
library(boot)
a=c(50,30,40)
A1=matrix(c(2,3,5,5,2,4),nrow=2, byrow=TRUE)
b1=c(100,80)
simplex(a,A1,b1,maxi = TRUE)


const_type <- rep("<=",2)
linprog <- lp("max",a,A1,const_type,b1)
linprog$solution
