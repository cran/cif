#' Wald test for Multivariate Least-Squares regression
#'
#' @param mlsresults output of mls, mlsresults<-mls(y, x, df_flag)
#' @author P. Paruolo
#' @description  Wald test for multivariate Least-Squares regression
#' @importFrom stats pchisq
# @usage wald<-Wald.mls(mlsresults)
#' @references Berta et al. 2020
#' @return wald table of Wald tests on significance of single regressors and pvalues based on chi square distribution
#' @export

Wald.mls <- function(mlsresults)
{ # Wald test on significance of single x variables in
  # y     =  x      \beta   +  u
  # (n x g) (n x k) (k x g) (n x g)
  # returns wald  = trace(c'(x'x)^{-1}c)^{-1}c'\beta\Omega{-1}\beta'c) and pval
  # where c selects each x, one at the time.
  mB<-t(mlsresults$coeffs); XXinv<-mlsresults$XXinv; # mB =t (\beta), x'x^{-1}
  Oinv<-solve(mlsresults$Omega);                     # Inv Omega
  ny<-nrow(Oinv); nx<-nrow(XXinv);                   # ny nx
  Waldtable<-matrix(nrow=2,ncol=nx);                 # Table of results. 1st row: Wald stat, 2nd row: pval
  colnames(Waldtable)<-colnames(mB);                 # x var names
  rownames(Waldtable)<-c("Wald stat", "p value")     # row names
  for(j in (1:nx)){
  b<-as.matrix(mB[,j])
  Wj<-t(b)%*%Oinv%*%b/XXinv[j,j]
  pvj<-pchisq(Wj,ny,lower.tail = FALSE)
  Waldtable[1,j]<-Wj; Waldtable[2,j]<-pvj}
  # write output
  out <- Waldtable
  # result:
  return(out)
}
