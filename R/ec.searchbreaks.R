#' search for breaks dates for given q (=1,2,3,4)
#'
#' @param qse q: number of (additional) breaks, s: start date for search, e: end date for search
#' @param ymat matrix of lags
#' @param npl  n, p, nlag
#' @param befpn begtrim, endtrim, nforecast, npred
#' @param ndet order of the model d(i,j)
#' @param gfillmin gfill value
#' @param fixed vector of breaks to be taken as fixed (not between s=start and e=end)
#' @author P. Paruolo
#' @description Search for location of break points in 1st-stage of Engle-Granger
# @usage out<-ec.searchbreaks(qse,ymat,npl,befpn,ndet,gfillmin=11,fixed=NA)
#' @references Berta et al. 2020
#' @return out list with break dates and values of regression average sum of squares
#' @export

ec.searchbreaks <- function(qse,ymat,npl,befpn,ndet,gfillmin=10,fixed=NA){
  qadd<-qse[1]; start<-qse[2];  enddate<-qse[3]    # qadd, start, enddate
  if(anyNA(fixed)==T){qfix<-0}else{
    qfix<-length(fixed)}                           #fixed
  n<-npl[1]; p<-npl[2]; nlag<-npl[3]               # n, p, nlag
  begtrim <- befpn[1]; endtrim <- befpn[2]         # begtrim, endtrim
  nforecast<-befpn[3]; npred<-befpn[4]             # nforecast, npred,
  twoq<-2*(qadd+qfix); twoqp1 <-twoq+1;            # 2q+1,q+1, q
  ntot <-n+npred; twop <-2*p; pp1<-p+1             # ntot, 2p, p+1
  gfill<-max(gfillmin,twoq+2)                      # choice of gfill REM: g \geq 2(q+1)+1
# ------ data
  longdata <- as.matrix(ymat);                     # total set of obervations
  if(dim(longdata)[1]!=n){cat("wrong number of lines in ymat")}
  ndet1<-ndet[1]                                   # det terms in EG1
  mynames<-colnames(longdata[,pp1:twop])           # extracts names of lagged y
  ylev<-longdata[,pp1:twop]                        # y in levels at t-1
# -------------- Engle Granger 1st stage regression preparation ------------------
  yEG <- as.matrix(ylev[,1,drop=FALSE])              # l.h.s. y
  xEG <- as.matrix(ylev[,-1,drop=FALSE])             # start of r.h.s x
  init <- max(begtrim,nlag,2)                        # min trim is 2
  eini<-init+1; eend<-n-endtrim;                     # est initial and est last dates
  estsample <- (eini:eend)                           # est sample
  nini<-max(0,start-init); nend<-max(eend-init,enddate-init)     # search initial and search last dates within estsample
  yE<- yEG[estsample]                                # trims y to est sample
  myT<- n-endtrim-init                               # sample size
  mycons<- rep(1,n); mytrend<-(1:n)                  # const, trend
# ------------ breaks ------
  if(ndet1==1){xEG<-cbind(xEG,mycons); colnames(xEG)[2]<-"const"}
  if(ndet1==2){xEG<-cbind(xEG,mycons,mytrend);
  colnames(xEG)[2:3]<-c("const","trend")}
  if(qfix>0){
    myCON <-matrix(1,n,qfix); myTRE <-matrix(mytrend,n,qfix)   # constant and trend with breaks
    for(i in (1:qfix)){aa<-(1:(fixed[i]-1)); myTRE[aa,i]<-myCON[aa,i]<-0}
    if(ndet1==1){xEG<-cbind(xEG,myCON); colnames(xEG)[3:(2+qfix)]<-paste0("const",fixed)}
    if(ndet1==2){xEG<-cbind(xEG,myCON,myTRE);
    colnames(xEG)[4:(3+2*qfix)]<-c(paste0("const",fixed),paste0("trend",fixed))}
  }
  nreg0<-ncol(xEG); #number of fixed regressors
# --- additional breaks ------------------
  brstart<-nini+gfill+1                        # initial break date for search
  brend<-rev(nend-(1:qadd)*(gfill+1))        # final break dates search
  myCON <-matrix(1,n,qadd); myTRE <-matrix(mytrend,n,qadd)   # constant and trend with breaks
  colnames(myCON)<-paste0("const_",(1:qadd))    # names constants
  colnames(myTRE)<-paste0("trend_",(1:qadd))    # names trends
  if(ndet1==1){xEG<-cbind(xEG,myCON)}
  if(ndet1==2){xEG<-cbind(xEG,myCON,myTRE)}
  cpos<-nreg0+(1:qadd)                          # column position of broken constants
  tpos<-nreg0+((qadd+1):(2*qadd))               # column position of broken trends
  xE<-xEG[estsample,];
  mytrend<-mytrend[estsample];
# ------ ssmat -----------
  vsize<-listsize(nend,gfill,nini)              # number of models for q=1,2,3,4
  numcases<-vsize[qadd]
  ssmat<-matrix(nrow=numcases+1,ncol=qadd+1)   # matrix with break dates and ave RSS ## (adding one extra row)
  colnames(ssmat)<-make.names(c(paste0("break_",(1:qadd)),"ave SSR"))
# ============= ndet1 = 2 =====================================
# ------------- ndet1 = 2 and qadd=1 --------------------------
if(qadd==1&&ndet1==2){mycount<-0;
      for(i1 in (brstart:brend[1])){mycount<-mycount+1;
        cvec<-c(rep(0,i1),rep(1,myT-i1));tvec<-cvec*mytrend;
        xE[,cpos[1]]<-cvec;xE[,tpos[1]]<-tvec;
        outEG <- mls(yE,xE); fval<-outEG$Omega
        ssmat[mycount,]<-c(i1+init,fval)}
}
# ------------- ndet1 = 2 and qadd=2 --------------------------
if(qadd==2&&ndet1==2){mycount<-0;
  for(i2 in (brstart:brend[1])){
    cvec<-c(rep(0,i2),rep(1,myT-i2));tvec<-cvec*mytrend;
    xE[,cpos[1]]<-cvec;xE[,tpos[1]]<-tvec;
    for(i1 in ((i2+gfill+1):brend[2])){mycount<-mycount+1;
    cvec<-c(rep(0,i1),rep(1,myT-i1));tvec<-cvec*mytrend;
    xE[,cpos[2]]<-cvec;xE[,tpos[2]]<-tvec;
    outEG <- mls(yE,xE); fval<-outEG$Omega
    ssmat[mycount,]<-c(i2+init,i1+init,fval)}}
}
# ------------- ndet1 = 2 and qadd=3 --------------------------
if(qadd==3&&ndet1==2){mycount<-0;
      for(i3 in (brstart:brend[1])){
        cvec<-c(rep(0,i3),rep(1,myT-i3));tvec<-cvec*mytrend;
        xE[,cpos[1]]<-cvec;xE[,tpos[1]]<-tvec;
        for(i2 in ((i3+gfill+1):brend[2])){
          cvec<-c(rep(0,i2),rep(1,myT-i2));tvec<-cvec*mytrend;
          xE[,cpos[2]]<-cvec;xE[,tpos[2]]<-tvec;
          for(i1 in ((i2+gfill+1):brend[3])){mycount<-mycount+1;
            cvec<-c(rep(0,i1),rep(1,myT-i1));tvec<-cvec*mytrend;
            xE[,cpos[3]]<-cvec;xE[,tpos[3]]<-tvec;
            outEG <- mls(yE,xE); fval<-outEG$Omega
            ssmat[mycount,]<-c(i3+init,i2+init,i1+init,fval)}}}
}
# ------------- ndet1 = 2 and qadd=4 --------------------------
if(qadd==4&&ndet1==2){mycount<-0;
    for(i4 in (brstart:brend[1])){
      cvec<-c(rep(0,i4),rep(1,myT-i4));tvec<-cvec*mytrend;
      xE[,cpos[1]]<-cvec;xE[,tpos[1]]<-tvec;
    for(i3 in ((i4+gfill+1):brend[2])){
      cvec<-c(rep(0,i3),rep(1,myT-i3));tvec<-cvec*mytrend;
      xE[,cpos[2]]<-cvec;xE[,tpos[2]]<-tvec;
      for(i2 in ((i3+gfill+1):brend[3])){
        cvec<-c(rep(0,i2),rep(1,myT-i2));tvec<-cvec*mytrend;
        xE[,cpos[3]]<-cvec;xE[,tpos[3]]<-tvec;
        for(i1 in ((i2+gfill+1):brend[4])){mycount<-mycount+1;
          cvec<-c(rep(0,i1),rep(1,myT-i1));tvec<-cvec*mytrend;
          xE[,cpos[4]]<-cvec;xE[,tpos[4]]<-tvec;
          outEG <- mls(yE,xE); fval<-outEG$Omega             #
          ssmat[mycount,]<-c(i4+init,i3+init,i2+init,i1+init,fval)}}}}
}
# ============= ndet1 = 1 =====================================
  # ------------- ndet1 = 1 and qadd=1 --------------------------
  if(qadd==1&&ndet1==1){mycount<-0;
  for(i1 in (brstart:brend[1])){mycount<-mycount+1;
  cvec<-c(rep(0,i1),rep(1,myT-i1));
  xE[,cpos[1]]<-cvec;
  outEG <- mls(yE,xE); fval<-outEG$Omega
  ssmat[mycount,]<-c(i1+init,fval)}
  }
  # ------------- ndet1 = 1 and qadd=2 --------------------------
  if(qadd==2&&ndet1==1){mycount<-0;
  for(i2 in (brstart:brend[1])){
    cvec<-c(rep(0,i2),rep(1,myT-i2));
    xE[,cpos[1]]<-cvec;
    for(i1 in ((i2+gfill+1):brend[2])){mycount<-mycount+1;
    cvec<-c(rep(0,i1),rep(1,myT-i1));
    xE[,cpos[2]]<-cvec;
    outEG <- mls(yE,xE); fval<-outEG$Omega
    ssmat[mycount,]<-c(i2+init,i1+init,fval)}}
  }
  # ------------- ndet1 = 1 and qadd=3 --------------------------
  if(qadd==3&&ndet1==1){mycount<-0;
  for(i3 in (brstart:brend[1])){
    cvec<-c(rep(0,i3),rep(1,myT-i3));
    xE[,cpos[1]]<-cvec;
    for(i2 in ((i3+gfill+1):brend[2])){
      cvec<-c(rep(0,i2),rep(1,myT-i2));
      xE[,cpos[2]]<-cvec;
      for(i1 in ((i2+gfill+1):brend[3])){mycount<-mycount+1;
      cvec<-c(rep(0,i1),rep(1,myT-i1));
      xE[,cpos[3]]<-cvec;
      outEG <- mls(yE,xE); fval<-outEG$Omega
      ssmat[mycount,]<-c(i3+init,i2+init,i1+init,fval)}}}
  }
  # ------------- ndet1 = 1 and qadd=4 --------------------------
  if(qadd==4&&ndet1==1){mycount<-0;
  for(i4 in (brstart:brend[1])){
    cvec<-c(rep(0,i4),rep(1,myT-i4));
    xE[,cpos[1]]<-cvec;
    for(i3 in ((i4+gfill+1):brend[2])){
      cvec<-c(rep(0,i3),rep(1,myT-i3));
      xE[,cpos[2]]<-cvec;
      for(i2 in ((i3+gfill+1):brend[3])){
        cvec<-c(rep(0,i2),rep(1,myT-i2));
        xE[,cpos[3]]<-cvec;
        for(i1 in ((i2+gfill+1):brend[4])){mycount<-mycount+1;
        cvec<-c(rep(0,i1),rep(1,myT-i1));
        xE[,cpos[4]]<-cvec;
        outEG <- mls(yE,xE); fval<-outEG$Omega             #
        ssmat[mycount,]<-c(i4+init,i3+init,i2+init,i1+init,fval)}}}}
  }
# ============= choose min ====================
brdates<-ssmat[which.min(ssmat[,qadd+1]),(1:qadd)]
# ---------- results -------------
out<-list(); out$brdates<-brdates;          # dates of additional breaks
out$fixed<-fixed;                           # dates of fixed breaks
out$ssmat<-ssmat; out$nbreaks<-qadd+qfix;   # matrix of results, num of breaks (total)
out$qadd<-qadd; out$qfix<-qfix;             # num of additional breaks, num of fixed breaks
# ---------- return -------------
return(out)
}


