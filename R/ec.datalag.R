#' prepares Dy y_1 Dy_1 ... Dy_{nlag-1} for estimation
#'
#' @param y is the data matrix of variables in the VAR
#' @param nlag is the number of lags in the VAR (min = 2)
# @usage ymat <- ec.datalag(y,nlag=4)
#' @author P. Berta, P. Paruolo, S. Verzillo, PG. Lovaglio
#' @description Prepares data for estimation
#' @references Berta et al. 2020
#' @return ymat contains the folloiwing columns {Dy y_1 Dy_1 ... Dy_{nlag-1}}
#' @export

ec.datalag <- function(y,nlag=4)
{if(nlag<=1){nlag<-2}
a<-dim(as.matrix(y)); n<-a[1]; p<-a[2]
Dy <- diffe(y)                              # differences
ylag <- lagn(y,1)                           # lagged levels
ymat<-cbind(Dy,ylag)                        # start
for(j in 1:(nlag-1))                        # lags of the differences
{ymat<-cbind(ymat,lagn(Dy,j))}
# result:
return(ymat)
}


