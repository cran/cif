#' plots forecasts of difference with confidence bars
#'
#' @param obj output of ec.main
#' @param whichseries series number
#' @param nsigma how many standard deviations in confidence bars
#' @param xvec vector of dates to place on x axis
#' @param yvec vector of exp(y) values to display on y axis
#' @param cal calendar vector
#' @param lar length of arrows in error bars
#' @param ... other plot parameters
#' @author P. Paruolo,
# @usage ec.gfd(obj,whichseries=1,nsigma=3,xvec,yvec,cal,lar=0.025,...)
#' @description plots forecasts of difference with confidence bars
#' @importFrom graphics plot lines abline arrows axis text
#' @importFrom lubridate month day as_date
#' @return does not return output, just creates a graph
#' @export

ec.gfd <- function(obj,whichseries=1,nsigma=3,xvec,yvec,cal,lar=0.025,...)
{p<-ncol(obj$Omega); n<-nrow(obj$ymat); nfor <-obj$befpn[3];  # p,n,nforecast
endtrim <-obj$befpn[2]; n_eff<-n-endtrim+nfor; mysample<-(1:n_eff)# mysample
x <- obj$cal[mysample]; vpos <-c(0,p,2*p)+whichseries;       # cal, positions of series
pp12p <- (p+1):(2*p)                                         # second block of columns
larr<-0.025                                                  # length of arrows
# -------------- differences ----------
shortafdd<-obj$afdd[1:(nfor+1),]; yactd <-obj$ymat[mysample,(1:p)];  # forecasts diff
longford<-rbind(matrix(NA,n_eff-nfor-1,3*p),shortafdd)           # long forecast diff
longford[(1:(n_eff-nfor)),pp12p]<-obj$afd[(1:(n_eff-nfor)),pp12p]    # assign fit diff in sample
ygd <- cbind(yactd,longford)[,-pp12p]                        # drop second block of cols
yfd <- ygd[,vpos]; yfd[n_eff-nfor,3]<-NA                         # select series
yfd.min<-yfd[,2]-nsigma*yfd[,3]; yfd.max<-yfd[,2]+nsigma*yfd[,3]  # min max forecast interval diff
# --------------- plot settings ----------------
labc<-c("blue","red");labtype<-rep(1,2); small<-0.05;
xtick<-as_date(x[xvec])
xticklab<-paste0(month(xtick,label=T)," ",day(xtick))
ytick<-c(log(yvec)); yticklab<-yvec
# ------------------- plot forecast differences ---------------------
myscale<-cbind(yfd[,c(1,2)],yfd.min,yfd.max)
sc_min<-min(min(myscale,na.rm = TRUE))
sc_max<-max(max(myscale,na.rm = TRUE))
smax<-(1+small*sign(sc_max))*sc_max
smin<-(1-small*sign(sc_min))*sc_min
myylimd <-c(smin,smax)
# ------------ plot 2 begins ------------
plot(x,yfd[,2],type="l",col=labc[2],ylim=myylimd,lty=labtype[2],ylab="",xlab="",xaxt="n")
axis(side=1, at=xtick, labels = FALSE)
text(x=xtick,  par("usr")[3],
     labels = xticklab, srt = 0, pos = 1, offset=1, xpd = TRUE)
lines(x,yfd[,1],type="l",lty=labtype[1], col=labc[1])
abline(v=c(cal[min(obj$estsample)],
           cal[max(obj$estsample)]), col=c("blue", "blue"), lty=c(2,2), lwd=c(1, 1))
arrows(x0=x, y0=yfd.min, x1=x, y1=yfd.max, code=3, angle=90, length=lar)
}

