#' produces predictions for the VECM via its VAR companion form
#'
#' @param est output from estimation by ec.EG1.R
#' @param det1 deterministic matrix of constant(s) and trend(s)
#' @param det2 deterministic matrix of seasonals and point dummies
#' @param ymat matrix of lags
#' @param npl  n, p, nlag
#' @param befpn begtrim, endtrim, nforecast, npred, nhstar
#' @param ndet order of the model d(i,j)
#' @param cal calendar, should match the number of rows in ymat
#' @param kval how many se to use, default kval= 1.959964
# @usage ec.predict(est,det1,det2,ymat,npl,befp,ndet,cal)
#' @author P. Berta, P. Paruolo, S. Verzillo, PG. Lovaglio
#' @description Predicts both in-sample (1 step ahead) and out-of-sample (1 step ahead and dynamic forecasts)
#' @references Berta et al. 2020
#' @return list with contains:
#' afl (actual and 1 step ahead fitted levels)
#' afd (actual and 1 step ahead fitted differences)
#' fit (1 step ahead fit)
#' dynpred (dynamic predictions)
#' mAt mB (companion matrix and selection of it)
#' Sigmah (Sigmah for dyn forecasts)
#' forstartdate (starting date for dyn forecast)
#' outcal  (dates for the prediction)
#' h1star (h1star)
#' cspred (table with change in sign of pred for Dx_1)
#' indexfa (indices of forecast accuracy)
#' @export

ec.predict <- function(est,det1,det2,ymat,npl,befpn,ndet,cal,kval=1.959964){
  n<-npl[1];p<-npl[2];nlag<-npl[3]                # n, p, nlag
  names(befpn)<-rep(NULL,6);
  begtrim <- befpn[1]; endtrim <- befpn[2]        # begtrim, endtrim
  nforecast<-befpn[3];npred<-befpn[4];nhstar<-befpn[5] # nforecast, npred, nhstar
  twoqp1 <-ncol(det1); qp1<-twoqp1/2; q <-qp1-1   # 2q+1,q+1, q
  ntot <-nrow(det1); twop <-2*p; pp1<-p+1         # ntot, 2p, p+1
  Omega <- est$Omega                              # Omega
  est_se<-as.matrix((diag(Omega)^.5))             # in-sample se
  rownames(est_se)<-paste("se",rownames(est_se))  # in-sample se names
  nump <- endtrim+npred                           # nump: total number of predictions
if(nump <= 0){stop("endtrim+npred<=0\n")}           # check if nump > 0
if(nforecast > endtrim){stop("nforecast > endtrim")} # check  npred > nforecast
  obsdata <- as.matrix(ymat);                     # total set of obervations
  longafdl<-longafdd<-longaffore<-matrix(NA,n+npred,(3*p)); # initialize longafdd, longafdl, longaffore
  ncheck <-length(cal); # number of time period
if(n!= ncheck){stop("rows in ymat not equal to cal\n")} # check if nrows don't match
  longcal <- seq(from=cal[1], by="days",length.out = ntot)
  mynames <- colnames(obsdata[,(p+1):twop])        # extracts names of lagged y
  mynames<- vapply(strsplit(mynames,"_"), `[`, 1, FUN.VALUE=character(1))
  ylev <- obsdata[,1:p]+obsdata[,(p+1):twop];      # y in levels
  colnames(ylev) <- mynames                        # assigns names to y in lev
  mAt <- t(ec.companion(est,p=p,nlag=nlag))        # builds transpose companion matrix
  roots <- eigen(mAt,symmetric = FALSE,only.values = TRUE) # roots of the companion matrix
  s <- ncol(mAt)                                   # dim state vector
# -------------  1-step ahead fitted and predicted over the entire sample
  mB <- mAt[,1:twop]                               # predictor for Dx_t and x_{t-1}
  if(anyNA(det2)==T){detlong<-det1;detest<-detlong[1:n,] # deterministics, no det2
  mut<-t(est$xi)}else{
    detlong<-cbind(det1,det2);detest<-detlong[1:n,]# deterministics, long and short (estimation sample)
    mut<-t(cbind(est$xi,est$nu))                   # mu transposed
  }
  eyep <- diag(rep(1,p)); msum <- rbind(eyep,eyep) # I_p and (I_p,I_p)'
  mK <- cbind(rbind(eyep,matrix(0,p,p)),msum)      # K matrix
  mypre<-obsdata[,1:s,drop=FALSE] %*% mB + cbind(detest %*% mut,matrix(0,n,p))
  fit <- lagn(mypre,1) # fitted over the whole sample
  colnames(fit)<-paste0(colnames(obsdata[,1:twop]),"_pred1",sep="")
  fl <- fit%*%msum; colnames(fl)<-paste0(mynames,"_pred1",sep="")     # fitted levels
  fd <- fit[,1:p]; colnames(fd)<-paste0("D",mynames,"_pred1",sep="")  # fitted differences
  act <- obsdata[,1:twop]                        # actual
  al <- act%*%msum; colnames(al)<-paste0(mynames,"_actual",sep="")    # actual levels
  ad <- act[,1:p]; colnames(ad)<-paste0("D",mynames,"_actual",sep="") # actual differences
  vse <- matrix(1,n,1)%*%t(est_se)                # matrix of est se
  colnames(vse)<-rownames(est_se)                 # names for se
  afl <- cbind(al,fl,vse)                         # actual and fitted levels with se
  afd <- cbind(ad,fd,vse)                         # actual and fitted differences with se
  if(npred==0){longafd<-afd; longafl<-afl}else{   # longafl, longafd
    myfill<-matrix(data=NA,npred,3*p)
    longafl<-rbind(afl,myfill); longafd<-rbind(afd,myfill)}
# -------------  multi step ahead forecast out of sample (after t0=n-entrim)
  t0 <- n-endtrim; forstartdate <-longcal[t0]     # t0 is the last point before prediction
  outcal <- longcal[t0:(t0+nump)]                 # calendar for multi step ahead pred
  lastvalue <- obsdata[t0,1:s,drop=FALSE]         # last available data vector & obs to forecast
  dynpred <- matrix(0,nump+1,s)                   # matrix of predictions
  fcse <- matrix(0,nump+1,twop)                   # matrix of forecast standard errors
  Sigmah <- array(0, dim = c(twop,twop,nump+2))   # initialize Sigmah
  colnames(dynpred) <- colnames(lastvalue)        # variable names
  mAtj<-diag(1,s,s)                               # initialize t(A^{h-1}) j<-1
  dynpred[1,]<-lastvalue[1,]                      # initialize dynpred
  filler<-matrix(0,1,s-p)                         # filler of 0s
  for(j in 1:nump){                               # multi-step prediction
  dynpred[j+1,]<-dynpred[j,]%*%mAt+cbind(detlong[t0+j,] %*% mut,filler)
  mVjt <- mAtj[1:p,1:twop] %*% mK
  Sigmah[,,j+1]<-Sigmah[,,j]+t(mVjt)%*%Omega%*%mVjt     # Sigma_h
  fcse[j+1,]<-(diag(Sigmah[,,j+1]))^.5            # fcse for Dx and x
  mAtj<-mAtj %*% mAt                              # t(A^{h-1})
  }
  fdd <- dynpred[,1:p]                            # dynamic forecast for differences
  colnames(fdd)<-paste0("D",mynames,"_dynfor",sep="")   # names
  fdl <- dynpred[,1:twop] %*% msum                # dynamic forecast for levels
  colnames(fdl)<-paste0(mynames,"_dynfor",sep="") # names
  colnames(fcse)<- c(paste0("D",mynames,"_fcse",sep=""),paste0(mynames,"_fcse",sep=""))
  if(npred>0){ald <- rbind(al[t0:n,],matrix(data=NA,npred,p)) # creates extended series of actual lev
  add <- rbind(ad[t0:n,],matrix(data=NA,npred,p)) # creates extended series of actual diff
  }else{ald <- al[t0:n,]; add <- ad[t0:n,]}
# ------ reporting --------------
  afdd <- cbind(add,fdd,fcse[,1:p])               # actual and forecast dynamic differences with fcse
  afdl <- cbind(ald,fdl,fcse[,(p+1):twop])        # actual and forecast dynamic levels with fcse
  affore <- afdl[1:(nforecast+1),];               # affore
  longaffore[t0:(t0+nforecast),]<-afdl[1:(nforecast+1),]   # longaffore
  longafdl[t0:(n+npred),]<-afdl[1:(nump+1),]      # longafdl
  longafdd[t0:(n+npred),]<-afdd[1:(nump+1),]      # longafdd
  colnames(longafdl)<-colnames(longaffore)<-colnames(afdl) # colnames
  colnames(longafdd)<-colnames(afdd)              # colnames
  cspred <- crossing(afdd[,p+1])                  # at what horizons forecasts of Dx_1 change sign
  negseq  <- as.matrix(cspred[,cspred[2,] < 0])   # select neg sign
  posneg  <- negseq[1,negseq[3,] > nhstar]        # at what horizons forecasts of Dx_1 change sign
  if(length(posneg)==0){h1star<-NA}else{h1star<-outcal[posneg[1]]}    # h1star
# --------- calculate indices of forecast performance -------------------
  rwsigma<-est$rwsigma; rwabsmean<-est$rwabsmean
  if(anyNA(affore)==T||nforecast==0){tablefap <- NA
    cat("Warning: indices of forecast accuracy cannot be calculated\n")
    }else{
    tablefap <- ec.ifp(affore,rwsigma,rwabsmean,kval)} # indices of forecast accuracy, 95% for Gaussian
# --------- write output
out <-list()
out$afl<-afl; out$afd<-afd;                       # 1 step ahead prediction with se
out$longafl<-longafl;   out$longafd<-longafd;     # long versions of afl and afd with se
out$affore<-affore;  out$longaffore<-longaffore   # afdl[1:(nforecast+1),] of lengthnforecast+1
out$afdl <- afdl;  out$afdd <- afdd               # actual fitted dynamic prediction with fcse
out$longafdl <-longafdl;out$longafdd <-longafdd;  # longafdd & longafdl, of lenght n+npred
out$fit <- fit;                                   # 1 step ahead prediction for levels and differences
out$dynpred <- dynpred;                           # raw dynamic output
out$mAt <- mAt; out$mB <- mB; out$Sigmah <- Sigmah # raw output
out$forstartdate <- forstartdate                  # starting date for dyn forecast
out$outcal <- outcal                              # dates for the prediction
out$h1star <- h1star                              # h1star
out$cspred <- cspred                              # table with change in sign of pred for Dx_1
out$tablefap <- tablefap                          # table with indices of forecast accuracy
out$affore <- affore                              #
out$t0 <-t0                                       # start date of dyn pred is t0+1
# result:
return(out)
}

