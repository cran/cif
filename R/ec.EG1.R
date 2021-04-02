#' estimates the VECM with the 2-stage procedure of Engle & Granger
#'
#' @param det1 deterministic matrix of constant(s) and trend(s)
#' @param det2 deterministic matrix of seasonals and point dummies
#' @param ymat matrix of lags
#' @param npl  n, p, nlag
#' @param befpn begtrim, endtrim, nforecast, npred
#' @param ndet order of the model d(i,j)
#' @param drop1 selection of det1 regressors in first stage to drop
#' @param drop2 selection of det1 regressors in second stage to drop
#' @author P. Berta, P. Paruolo, S. Verzillo, P.G. Lovaglio
#' @description Estimates the EC with EG. Cointegration rank fixed at 1
# @usage ec.EG1(det1,det2,ymat,npl,ndet,drop1,drop2)
#' @references Berta et al. 2020
#' @return out a list with estimates
#' @export

ec.EG1 <- function(det1,det2,ymat,npl,befpn,ndet,drop1=NA,drop2=NA){
  n<-npl[1];p<-npl[2];nlag<-npl[3]                 # n, p, nlag
  begtrim <- befpn[1]; endtrim <- befpn[2]         # begtrim, endtrim
  nforecast<-befpn[3]; npred<-befpn[4]             # nforecast, npred,
  twoqp1 <-ncol(det1); qp1<-twoqp1/2; q <-qp1-1    # 2q+1,q+1, q
  ntot <-nrow(det1); twop <-2*p; pp1<-p+1          # ntot, 2p, p+1
  boodet2<-anyNA(det2)                             # boolean for missing det2
# ------ data
  longdata <- as.matrix(ymat);                     # total set of obervations
  if(dim(longdata)[1]!=n){cat("wrong number of lines in ymat")}
  ndet1<-ndet[1]; ndet2<-ndet[2];                  # det terms in EG1 and EG2
  mynames<-colnames(longdata[,pp1:twop])           # extracts names of lagged y
  ylev<-longdata[,pp1:twop]                        # y in levels at t-1
# -------------- Engle Granger 1st stage regression ------------------
yEG <- ylev[,1,drop=FALSE]                         # l.h.s. y
if(ndet1==1){ivec <- (1:qp1)}                      # selection vector with 1_vec
if(ndet1==2){ivec <- (1:(2*qp1))}                  # selection vector with 1_vec, t_vec
if(anyNA(drop1)==F){ivec <-ivec[-drop1]}           # deletes determ. components from first stage
if(ndet1==0){xEG <- as.matrix(ylev[,-1,drop=FALSE])
  ivec<-NA}else{
  xEG <- cbind(ylev[,-1,drop=FALSE],det1[1:n,ivec])}
init <- max(begtrim,nlag,2)                        # min trim is 2
estsample <- (init+1):(n-endtrim)                  # est sample
yEGest <- yEG[estsample]                           # trims y to est sample
xEGest <- xEG[estsample,]                          # trims x to est sample
outEG <- mls(yEGest,xEGest)                        # EG First stage regression
afEG <- cbind(yEGest,outEG$fit)                    # EG 1st stage actual and fitted in est sample
vEG1step <- as.matrix(outEG$coeffs)                # coeffs
colnames(afEG)[2] <- paste0("fit EG1stage")        # names
colnames(vEG1step) <- paste0("EG1st eq. ",colnames(yEG),sep="")
ibeta<-1:(p-1); beta <- matrix(c(1,-vEG1step[ibeta])) # beta
longafEG <- cbind(yEG,xEG%*%vEG1step)              # EG 1st stage actual and fitted in full sample
colnames(longafEG)[2] <- paste0("f+f EG1stage")    # names
# myplot(longafEG)
if(anyNA(ivec)==T){betad<-NA}else{                 # betad
  betad<-matrix(0,(twoqp1),1); betad[ivec]<--vEG1step[-ibeta]
  rownames(betad)<-colnames(det1)}
ecmlong <- as.matrix(yEG-xEG%*%vEG1step)           # computes ecm for the whole sample
colnames(ecmlong)<-"ecm_1";                        # ecm
# -------------- Engle Granger 2nd stage regression ----------------------
Dy <- as.matrix(longdata[estsample,1:p])           # Second stage EG
eqnames <- paste("eq",colnames(Dy),sep=". ")       # dependent variables
colnames(Dy) <- eqnames                            # names
rwsigma <- (colMeans(Dy^2))^.5                     # sigma_i for RW model for all columns
rwabsmean <- colMeans(abs(Dy))                     # mean of abs(Dy)  all columns
xreg <- as.matrix(ecmlong[estsample])              # build regressors starting from ecm_1
colnames(xreg)<-"ecm_1"                            # names
if(nlag>1){xreg<-cbind(xreg,longdata[estsample,(twop+1):((nlag+1)*p)])} # lags
iposini<- (p*(nlag-1))+1;                          # number of positions to skip before det1 in 2nd EG regression
if(ndet2==0){ivec2 <- NA; iseldet12<-iposini}else{ # no det1 in second stage
  if(ndet2==1){ivec2 <- (1:qp1)}                   # only 1_vec in second stage
  if(ndet2==2){ivec2 <- (1:(2*qp1))}               # both 1_vec and t_vec in second stage
  if(anyNA(drop2)==F){ivec2 <-ivec2[-drop2]}       # drop components of det1 from second stage
  xreg<-cbind(xreg,det1[estsample,ivec2])          # adds det1 to xreg
  iseldet12<-iposini+ivec2                         # iseldet12 is the vector of positions of det1[,ivec2] in xreg
  colnames(xreg)[((nlag-1)*p+1+ivec2)]<-colnames(det1)[ivec2]  # names of cols of det1[,ivec2] in xreg
}
if(boodet2==F){xreg<-cbind(xreg,det2[estsample,])  # adds det2 to xreg
   det2col <-ncol(det2)}
outEC <- mls(Dy,xreg)                              # EG second stage regression
acoef <- t(outEC$coeffs)                           # all coefficients
wald <- Wald.mls(outEC)                            # test of joint significance of regressors
DafEG2<-cbind(Dy,outEC$fit)                        # actual and fitted in second stage
if(nlag>1){mPhi <- acoef[,2:iposini]               # mPhi
    }else{mPhi<-matrix(0,p,p)}
alpha <- acoef[,1,drop=FALSE]                      # alpha
# --------- xi zeta nu ------------------
zeta<-matrix(0,p,(twoqp1))                         # initialize zeta
colnames(zeta)<-colnames(det1)                     # names cols zeta
rownames(zeta)<-rownames(acoef)                    # names rows zeta
xi<-zeta;                                          # initialize xi
if(anyNA(ivec2)==F){zeta[,ivec2]<-acoef[,iseldet12]} # compute zeta
if(ndet1==0&&ndet2!=0){xi<-zeta}                   # xi = zeta
if(ndet1!=0&&ndet2==0){xi<-alpha%*%t(betad)}       # xi = alpha betad'
if(ndet1!=0&&ndet2!=0){xi<-alpha%*%t(betad)+zeta}  # xi = alpha betad'+zeta
if(boodet2==T){nu<-NA}else{                        # nu
  nu<-acoef[,((iseldet12[length(iseldet12)]+1):ncol(acoef))]}
# --------- write output -------------------
out <-list();
# --------- model specification
out$ndet <- ndet; out$det1names<-colnames(det1)    # d(i,j) choice ndet, names of cols in det1
out$det2names<-colnames(det2); out$drop1<-drop1    # names of cols in det2, selection of det1 in 1st stage
out$drop2 <- drop2; out$estsample <- estsample     # selection of det1 in 2nd stage, obs in estimantion sample
out$npl<-npl; out$befpn<-befpn                     # n,p,lag,begtrim,entrim,nforecast,npred,nhstar
# --------- first stage  --------------
out$beta <- beta; out$betad <- betad               # beta and betad
out$vEG1step <- vEG1step                           # first-stage-coeff
out$afEG <- afEG; out$longafEG <- longafEG         # actual and fitted 1st stage EG
out$ecm <- ecmlong                                 # ecm
# --------- second stage --------------
out$alpha <- alpha; out$mPhi <- mPhi;              # alpha, Phi
out$mPi <- alpha %*% t(beta)                       # Pi
out$acoef <- acoef; out$acoef_se <- t(outEC$semat) # 2nd stage coefficients
out$acoef_t <- t(outEC$tmat)                       # 2nd stage t ratios
out$acoef_pval <- t(outEC$pmat)                    # 2nd stage p values
out$res <- outEC$residuals;out$DafEG2 <- DafEG2    # 2nd stage residuals and actual and fitted
out$Omega <- outEC$Omega                           # Omega
out$nu <-nu; out$zeta<-zeta; out$xi<-xi            # nu, zeta, xi
out$wald<-wald;                                    # wald
# --------- additional output -------------------------
out$rwsigma <- rwsigma                             # sigma in RW model for all variables
out$rwabsmean <- rwabsmean                         # mean of abs(Dy) for all variables
out$EG1 <- outEG; out$EG2 <- outEC                 # 1st and 2nd stage regressions output
# result:
return(out)
}


