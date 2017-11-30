#########################################################
## Stat 202A - Homework 3
## Author: Ruoxuan Jia
## Date : 10/22/2017
## Description: This script implements QR decomposition,
## linear regression, and eigen decomposition / PCA
## based on QR.
#########################################################

#############################################################
## INSTRUCTIONS: Please fill in the missing lines of code
## only where specified. Do not change function names,
## function inputs or outputs. You can add examples at the
## end of the script (in the "Optional examples" section) to
## double-check your work, but MAKE SURE TO COMMENT OUT ALL
## OF YOUR EXAMPLES BEFORE SUBMITTING.
##
## Very important: Do not use the function "setwd" anywhere
## in your code. If you do, I will be unable to grade your
## work since R will attempt to change my working directory
## to one that does not exist.
#############################################################

##################################
## Function 1: QR decomposition ##
##################################

myQR <- function(A){

  ## Perform QR decomposition on the matrix A
  ## Input:
  ## A, an n x m matrix

  ########################
  ## FILL IN CODE BELOW ##
  ########################

  n<-dim(A)[1]
  m<-dim(A)[2]
  Q<-diag(n)

  for (k in 1:(m-1)){
    X<-matrix(0,nrow=n,ncol=1)
    X[(k:n),1]=A[(k:n),k]
    V=X
    V[k,1]=X[k,1]+sign(X[k,1])*norm(X,type="F")
    s=norm(V,type="F")
    u=V/s   #unit vector
    A=A-2*(u%*%(t(u)%*%A))
    Q=Q-2*(u%*%(t(u)%*%Q))

  }


  ## Function should output a list with Q.transpose and R
  ## Q is an orthogonal n x n matrix
  ## R is an upper triangular n x m matrix
  ## Q and R satisfy the equation: A = Q %*% R
  return(list("Q" = t(Q), "R" = A))

}


# now fitting a polynomial


# using QR decomposition directly
n    <- 100
p    <- 3

## Simulate data from our assumed model.
## We can assume that the true intercept is 0
X    <- matrix(rnorm(n * p), nrow = n)
beta <- matrix(1:p, nrow = p)
Y    <- X %*% beta + rnorm(n)
#myQR(a)
###############################################
## Function 2: Linear regression based on QR ##
###############################################

myLM <- function(X, Y){

  ## Perform the linear regression of Y on X
  ## Input:
  ## X is an n x p matrix of explanatory variables
  ## Y is an n dimensional vector of responses
  ## Do NOT simulate data in this function. n and p
  ## should be determined by X.
  ## Use myQR inside of this function

  ########################
  ## FILL IN CODE BELOW ##
  ########################

  ## Simulate data from our assumed model.
  ## We can assume that the true intercept is 0

  n=dim(X)[1]
  p=dim(X)[2]
  Z=cbind(rep(1,n),X,Y)
  result<-myQR(Z)
  R_m<-result$R
  R_1=R_m[1:(p+1) ,1:(p+1)]
  Y_1=R_m[(1:(p+1)),(p+2)]
  beta_ls=solve(R_1,Y_1)
  X=cbind(rep(1,n),X)
  sigma_square<-sum((Y- X %*% beta_ls)^2)/(n-(p+1))
  variance<-sigma_square*diag(solve(t(X)%*%X))
  sd_error<-sqrt(variance)
  R_std  <- coef(summary(lm(Y ~ X)))[,2]

  ## Function returns the 1 x (p + 1) vector beta_ls,
  ## the least squares solution vector
  return(list("beta_hat"=beta_ls,"standard_error"=sd_error))

}

a=myLM(X,Y)
print (a)
R_std  <- coef(summary(lm(Y ~ X)))[,2]
print(R_std)

#d <- qr(a)

## Simulate data from our assumed model.
## We can assume that the true intercept is 0

# using QR decomposition directly




##################################
## Function 3: PCA based on QR  ##
##################################



## Save R's linear regression coefficients
#R_coef  <- coef(lm(Y ~ X))

## Save our linear regression coefficients
#my_coef <- myLM(X, Y)



myEigen_QR <- function(A, numIter = 1000){

  ## Perform PCA on matrix A using your QR function, myQRC.
  ## Input:
  ## A: Square matrix
  ## numIter: Number of iterations

  ########################
  ## FILL IN CODE BELOW ##
  ########################

  r=dim(A)[1]
  c=dim(A)[2]

  V=matrix(sample(rnorm(n), r*r, TRUE),nrow=r,ncol=r)
  for (i in 1:numIter){
    Q<-myQR(V)["Q"]
    Q_ele<-(unlist(Q))
    Q_m<-matrix(Q_ele,ncol=r,byrow=FALSE)


    R<-myQR(V)["R"]
    R_ele<-(unlist(R))
    R_m<-matrix(R_ele,nrow=r,byrow=FALSE)
    A
    Q_m
    V=A %*% Q_m

  }
  ## Function should output a list with D and V
  ## D is a vector of eigenvalues of A
  ## V is the matrix of eigenvectors of A (in the
  ## same order as the eigenvalues in D.)
  return(list("D" = diag(R_m), "V" = Q_m))

}

testing_Linear_Regression <- function(){

  ## This function is not graded; you can use it to
  ## test out the 'myLinearRegression' function

  ## Define parameters
  n    <- 100
  p    <- 3

  ## Simulate data from our assumed model.
  ## We can assume that the true intercept is 0
  X    <- matrix(rnorm(n * p), nrow = n)
  beta <- matrix(1:p, nrow = p)
  Y    <- X %*% beta + rnorm(n)

  ## Save R's linear regression coefficients
  R_coef  <- coef(lm(Y ~ X))

  ## Save our linear regression coefficients
  my_coef <- myLM(X, Y)

  ## Are these two vectors different?
  sum_square_diff <- sum((R_coef - my_coef)^2)
  if(sum_square_diff <= 0.001){
    return('Both results are identical')
  }else{
    return('There seems to be a problem...')
  }

}
print(testing_Linear_Regression())
