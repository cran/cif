#' computes companion matrix of the VAR
#'
#' @param est is the output of ec.EG1.R
#' @param p (positive integer) is the dimension of the VAR
#' @param nlag  (positive integer) is the number of lags in the VAR
# @usage  mA<-ec.companion(est,p,nlag)
#' @author P. Berta, P. Paruolo, S. Verzillo, PG. Lovaglio
#' @description builds the companion matrix of the VAR
#' @references Berta et al. 2020
#' @return mA companion matrix
#' @export

ec.companion <- function(est,p=2,nlag=4)
{ # state vector: (Dx_t',x_{t-1}',Dx_{t-1}',...,Dx_{t-nlag+2}')'
  # indices     : (i1   ,i2      ,i3                          )
  twop<-2*p;
  if(nlag==1){s<-twop}else{s<-p*nlag}         # dim of the state vector
  i1<-1:p; i2<-(p+1):twop ;i3<-(twop+1):s     # index of3 blocks
  mA <- matrix(0,s,s)                         # initialize companion matrix
  eyep <- diag(rep(1,p))                      # I_p
# ---- first 2x2 blocks (nlag=2) --------------
  if(nlag==1){mA[i1,i1] <- est$mPi            # block (1,1)
  }else{mA[i1,i1] <- est$mPhi[,i1]+est$mPi}
  mA[i1,i2] <- est$mPi;                       # block (1,2)
  mA[i2,i1] <- eyep; mA[i2,i2] <- eyep;       # blocks (2,1), (2,2)
# ----- blocks for nlags > 2 ------------------
  if(nlag>2){
  mA[(twop+1):(3*p),i1] <- eyep                # block I_p in (3,1)
  mA[i1,i3] <- est$mPhi[,-i1]                  # block Gamma(2:k-1) in (1,3)
# -----  blocks for nlags > 3 ------------------
  if(nlag>3){                                  # block I_{p(k-3)} in (3,3)
  mA[(3*p+1):s,(twop+1):(s-p)] <- diag(rep(1,p*(nlag-3)))}
  }
# result:
return(mA)
}

