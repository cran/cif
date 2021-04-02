#' Forecast with Vector Error Correction Model
#'
#' @param y  matrix with time across rows and variables in columns
#' @param ndet vector of lenght 3, (i,j,q): i for EG1-st stage, j for EG-2nd stage, q number of breaks
#'             i,j=0 no deterministics
#'             i,j=1 constant
#'             i,j=2 constant and trend
#' @param nlag number of lags in the VAR
#' @param befpn begtrim, endtrim, nforecast, npred
#' @param breaks vector with observation numbers for T1,T2,...
#' @param booseas boolean =T if seasonal dummies, =F otherwise
#' @param pntdates vector with observation numbers for point dummies
#' @param drop1 selection of det1 regressors in first stage to drop
#' @param drop2 selection of det1 regressors in second stage to drop
#' @param cal calendar for the y matrix
#' @param kval how many se to use, default kval=1.959964
#' @author P. Berta, P. Paruolo, S. Verzillo, PG. Lovaglio
#' @description This function estimate VECM model. Selects begtrim and entrim period, define lag and run.
# @usage ec.main(y,ndet,nlag,befpn,breaks,booseas,pntdates,cal)
#' @references Berta et al. 2020
#' @return results Output contains the a set of estimates and forecasting results.
#' @export
#'

ec.main <- function(y,ndet=c(2,1),nlag,befpn,breaks=NA,booseas=NA,
                    pntdates=NA,drop1=NA,drop2=NA,cal,kval=1.959964){
  y<-as.matrix(y)               # y is taken as a matrix
  p<-ncol(y); n<-nrow(y)        # n p
  npl<-c(n,p,nlag);
  begtrim <- befpn[1]; endtrim <- befpn[2]          # begtrim, endtrim
  nforecast<-befpn[3]; npred<-befpn[4]              # nforecast, npred,
  ymat <- ec.datalag(y,nlag)                        # ymat
  onethirdss <- floor((n-endtrim-begtrim)/3)        # sample size / 3
  # ---------- optimize breaks -------------
  inidates<-floor(c(begtrim+onethirdss,begtrim+2*onethirdss))
  if(length(ndet)==3&&(ndet[3] %in% c(1,2,3,4))){q<-ndet[3]; # read q
  qse<-c(ndet[3],begtrim,n-endtrim); ndet12<-ndet[1:2]
  odates<-ec.searchbreaks(qse,ymat,npl,befpn,ndet12)# optimal dates
  breaks<-odates$brdates[1:q]                       # new break dates
  }
  # ---------- checks -------------
  if(anyNA(breaks)==F){breaks<- breaks[(breaks <n-endtrim) & (breaks>max(nlag,2,begtrim))]}
  if(anyNA(pntdates)==F){pntdates<- pntdates[(pntdates <n-endtrim) & (pntdates>max(nlag,2,begtrim))]}
  if(sum(ndet)==0){breaks<-NA}
  # ---------- data -------------
  det1 <- ec.datadet1(n,befpn,breaks)              # det1
  if(booseas==T||anyNA(pntdates)==F){det2 <- ec.datadet2(det1,booseas,pntdates)
  }else{det2<-NA}                                  # det2
  # ---------- estimation -------------
  est<-ec.EG1(det1,det2,ymat,npl,befpn,ndet,drop1,drop2)
  # ---------- forecasting -------------
  outB<-ec.predict(est,det1,det2,ymat,npl,befpn,ndet,cal,kval)
  # ---------- results -------------
  results<-list(); results$ymat<-ymat; results$breaks<-breaks;
  results$det1<-det1; results$det2<-det2;
  results<-append(results,est); results<-append(results,outB)
  results$cal<-cal;
  # ---------- return -------------
  class(results) <- "presize"
  return(results)
}

