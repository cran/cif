library(cif)

#Load data
qtr2<-iculomb
ssize <-nrow(qtr2)
cal <- seq(from=as.Date("2020-02-24"), by="days",length.out = ssize)

ic = log(qtr2$terapia_intensiva)
hws=log(qtr2$ricoverati_con_sintomi)

y <- cbind(ic,hws)
colnames(y) <- make.names(c("ic","hws"))
ssize <-nrow(qtr2)
cal <- seq(from=as.Date("2020-02-24"), by="days",length.out = ssize)

n<- nrow(y);  nlag=4; nhstar=0; booseas <- F;
# pntdates <-c(75,76,101);
pntdates <-NA; drop1<-NA; drop2<-NA;
npl<-c(nrow(y), ncol(y), nlag)           # collected n, p, nlag
breaks <- NA

ndet=c(2,1); begtrim=7;  endtrim=7;  nforecast=7; npred=0;
befpn<-c(begtrim,endtrim,nforecast,npred,nhstar) # collected specs for endtrim etc

longcal<-seq(from=cal[1], by="days",length.out = n+npred) # total calendar
# boosop  <- booseas||anyNA(pntdates)
det1 <- ec.datadet1(n,befpn,breaks)
det2 <- ec.datadet2(det1,booseas,pntdates)
ymat <- ec.datalag(y,nlag)
# odates<-ec.bfind(det1,ymat,npl,befpn,ndet)
oua1<-ec.main(y=y,ndet=ndet,nlag,befpn,breaks,booseas,
              pntdates,drop1=NA,drop2=NA,cal=cal)

class(oua1)

summary(oua1)



ix <- t(as.matrix(cbind(oua1$tablefap[c(1,3,4,6,9,10,13,14),3])))


barplot(ix,
        main = "Performance of the VEC model relative to the Random Walk",
        xlab = "",
        col = c("red","green","blue"),
        beside = TRUE
)
abline(h=1, col="blue", lty=2, lwd=1)


ix <- t(as.matrix(cbind(oua1$tablefap[c(1,3,4,6,9,10,13,14),3])))


rownames(ix)<-c("Model1")

myylim<-c(min(ix),max(ix)); yvec<-c(0.1,0.3,1,3,10,30)
plot(y=ix,x=seq(1:8),type="p",ylim=myylim,main="VEC model compared to the Random Walk",font.main=1,pch=19,col="blue",ylab="",xlab="",xaxt="n", yaxt="n");

xtick<-seq(1:8)
xticklab<-colnames(ix)
text(x=xtick, par("usr")[3],
     labels = xticklab, srt = 0, pos = 1, offset=1, xpd = TRUE)

ytick<- round(seq(min(ix),max(ix),length.out = 5),digits = 2)
yticklab<-round(seq(min(ix),max(ix),length.out = 5),digits = 2)
axis(side=2, at=ytick, labels = FALSE)
text(par("usr")[1], y=ytick,
     labels = yticklab, srt = 0, pos = 2, xpd = TRUE)

ytick1<- 1
yticklab1<-1
axis(side=2, at=ytick1, labels = FALSE)
text(par("usr")[1], y=ytick1,
     labels = yticklab1, srt = 0, pos = 2, xpd = TRUE)


abline(h=1, col="black", lty=2, lwd=1)




