/*
#########################################################
## Stat 202A - Homework 3
## Author: Ruoxuan Jia
## Date : 10/26/2017
## Description: This script implements QR decomposition,
## linear regression, and eigen decomposition / PCA 
## based on QR.
#########################################################
 
###########################################################
## INSTRUCTIONS: Please fill in the missing lines of code
## only where specified. Do not change function names, 
## function inputs or outputs. MAKE SURE TO COMMENT OUT ALL 
## OF YOUR EXAMPLES BEFORE SUBMITTING.
##
## Very important: Do not change your working directory
## anywhere inside of your code. If you do, I will be unable 
## to grade your work since R will attempt to change my 
## working directory to one that does not exist.
###########################################################
 
*/ 


# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


/* ~~~~~~~~~~~~~~~~~~~~~~~~~ 
 Sign function for later use 
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

// [[Rcpp::export()]]
double signC(double d){
  return d<0?-1:d>0? 1:0;
}



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   Problem 1: QR decomposition 
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */  
  

// [[Rcpp::export()]]
List myQRC(const mat A){ 
  
  /*
  Perform QR decomposition on the matrix A
  Input: 
  A, an n x m matrix (mat)

  #############################################
  ## FILL IN THE BODY OF THIS FUNCTION BELOW ##
  #############################################
  
  */ 
  
  
  
  // Function should output a List 'output', with 
  // Q.transpose and R
  // Q is an orthogonal n x n matrix
  // R is an upper triangular n x m matrix
  // Q and R satisfy the equation: A = Q %*% R
  int n = A.n_rows;
  int m = A.n_cols; 
  mat X; 
  mat Q;
  mat B=A; 
  mat W;
  Q.eye(n,n);
  
  for (int k=0; k<(m-1);k++){
    mat X=zeros(n,1);
    X(span(k,(n-1)),0)=B(span(k,(n-1)),k);
    mat V=X;
    V(k,0)=X(k,0)+signC(X(k,0))*(norm(X));
    double s=norm(V);
    mat u=V/s;
    B = B-2*(u*u.t()*B); 
    Q=Q-2*(u*u.t()*Q); 
  }
  List output;
  output["Q"] = Q.t();
  output["R"] = B;
  return(output);
  

}
  
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   Problem 2: Linear regression using QR 
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  
   
// [[Rcpp::export()]]
mat myLinearRegressionC(const mat X, const mat Y){
  mat beta_ls; 
  int n=X.n_rows;
  int p=X.n_cols;
  
  mat Q=ones(n,1);
  mat Z1=join_rows(Q,X); //cbind
  mat Z=join_rows(Z1,Y);
 
  mat R=myQRC(Z)["R"];
  mat R_1=R(span(0,p),span(0,p));
  mat Y_1=R(span(0,p),p+1);
  beta_ls=R_1.i()*Y_1;
  /*  
  Perform the linear regression of Y on X
  Input: 
  X is an n x p matrix of explanatory variables
  Y is an n dimensional vector of responses
  Do NOT simulate data in this function. n and p
  should be determined by X.
  Use myQRC inside of this function
  
  #############################################
  ## FILL IN THE BODY OF THIS FUNCTION BELOW ##
  #############################################
  
  */  
  
  
  // Function returns the 'p+1' by '1' matrix 
  // beta_ls of regression coefficient estimates
  return(beta_ls.t());
  
}  

/* ~~~~~~~~~~~~~~~~~~~~~~~~ 
 Problem 3: PCA based on QR 
 ~~~~~~~~~~~~~~~~~~~~~~~~~~ */


// [[Rcpp::export()]]
List myEigen_QRC(const mat A, const int numIter = 1000){
  
  
  
  
  /*
  Perform PCA on matrix A using your QR function, myQRC.
  Input:
  A: Square matrix
  numIter: Number of iterations
   
  #############################################
  ## FILL IN THE BODY OF THIS FUNCTION BELOW ##
  #############################################
   
 */
  //mat B=A; 
  mat D;
  int n=A.n_rows;
  int p=A.n_cols;
  mat V (n, p, fill::randu);
  mat Q=myQRC(V)["Q"];
  for (int i=0;i<numIter;i++){
     mat R=myQRC(V)["R"];
     mat Q=myQRC(V)["Q"];
     V=A*Q; 
  }
  mat R=myQRC(V)["R"];
  V=A*Q;
  D=R.diag(); 
  // Function should output a list with D and V
  // D is a vector of eigenvalues of A
  // V is the matrix of eigenvectors of A (in the 
  // same order as the eigenvalues in D.)
  List output;
  output["D"] = D;
  output["V"] = Q;
  return(output);

}


 

  
  // [[Rcpp::export()]]
  void test_eig (){
    mat B (5, 5, fill::randu);
    mat C = B.t() * B;
    vec eigval;
    mat eigvec;
    eig_sym(eigval, eigvec, C);
    cout << "rcpp_eigen_value: " << eigval;
    cout << "rcpp_eigen_vec: " << eigvec;
    
    List my_eval = myEigen_QRC(C);
    mat D = my_eval["D"];
    mat V = my_eval["V"];
    cout << "my_eigen_value: " << D;
    cout << "my_eigen_vec: " << V;
  }
  

