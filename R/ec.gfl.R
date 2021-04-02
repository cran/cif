#' plots level forecasts with confidence bars
#'
#' @param obj output of ec.main
#' @param whichseries series number
#' @param nsigma how many standard deviations in confidence bars
#' @param xvec vector of dates to place on x axis
#' @param yvec vector of exp(y) values to display on y axis
#' @param cal calendar vector
#' @param lar length of arrows in error bars
#' @param ... other plot parameters
#' @author P. Paruolo
# @usage ec.gfl(obj,whichseries=1,nsigma=3,xvec,yvec,cal,lar=0.025,...)
#' @description plots level forecasts with confidence bars
#' @importFrom graphics plot lines abline arrows axis text
#' @importFrom lubridate month day as_date
#' @return does not return output, just creates a graph
#' @export

ec.gfl <- function(obj,whichseries=1,nsigma=3,xvec,yvec,cal,lar=0.025,...)
{ p<-ncol(obj$Omega); n<-nrow(obj$ymat); nfor <-obj$befpn[3];  # p,n,nforecast
  endtrim <-obj$befpn[2]; n_eff<-n-endtrim+nfor; mysample<-(1:n_eff)# mysample
  x <- obj$cal[mysample]; vpos <-c(0,p,2*p)+whichseries;       # cal, positions of series
  pp12p <- (p+1):(2*p)                                         # second block of columns
  # ----------------- afdl = levels -------------------
  shortafdl<-obj$afdl[1:(nfor+1),];                            # forecasts
  yact <-obj$ymat[mysample,(1:p)]+obj$ymat[mysample,((p+1):(2*p))]; # actual levels
  mynames <- colnames(obj$ymat[,pp12p])                        # extracts names of lagged y
  colnames(yact)<- vapply(strsplit(mynames,"_"), `[`, 1, FUN.VALUE=character(1))
  longfor<-rbind(matrix(NA,n_eff-nfor-1,3*p),shortafdl)            # long forecast lev
  longfor[(1:(n_eff-nfor)),pp12p]<-obj$afl[(1:(n_eff-nfor)),pp12p]     # assign fit lev in sample
  yg <- cbind(yact,longfor)[,-pp12p]; yf <- yg[,vpos]          # select series
  yf[n_eff-nfor,3]<-NA                                             # se from 0 to NA at t
  yf.min<-yf[,2]-nsigma*yf[,3]; yf.max<-yf[,2]+nsigma*yf[,3]   # min max forecast interval lev
  # -------------- differences ----------
  labc<-c("blue","red");labtype<-rep(1,2); small<-0.05;
  mysc<-cbind(yf[,c(1,2)],yf.min,yf.max)
  slmin<-min(min(mysc,na.rm = TRUE)); slmax<-max(max(mysc,na.rm = TRUE))
  sclmax<-(1+small*sign(slmax))*slmax; sclmin<-(1-small*sign(slmin))*slmin
  myylim <-c(sclmin,sclmax)
  xtick<-as_date(x[xvec])
  xticklab<-paste0(month(xtick,label=T)," ",day(xtick))
  ytick<-c(log(yvec)); yticklab<-yvec
  # ---- plot begin ----
  plot(x,yf[,2],type="l",col=labc[2],ylim=myylim,lty=labtype[2],ylab="",xlab="",xaxt="n", yaxt="n")
  axis(side=1, at=xtick, labels = FALSE)
  text(x=xtick,  par("usr")[3],
       labels = xticklab, srt = 0, pos = 1, offset=1, xpd = TRUE)
  axis(side=2, at=ytick, labels = FALSE)
  text(par("usr")[1], y=ytick,
       labels = yticklab, srt = 0, pos = 2, xpd = TRUE)
  lines(x,yf[,1],type="l",lty=labtype[1], col=labc[1])
  abline(v=c(cal[min(obj$estsample)],
             cal[max(obj$estsample)]), col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1))
  # abline(v=cal[max(obj$estsample)], col="blue", lty=2, lwd=1)
  arrows(x0=x, y0=yf.min, x1=x, y1=yf.max, code=3, angle=90, length=lar)
  # ---- plot end ----
}

