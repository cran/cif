#' prepares deterministics D^(1)
#'
#' @param n is the number of obs in available data
#' @param befpn is a vector with (begtrim,endtrim,nforecast,npred,nhstar)
#' @param breaks is a vector of integers where the trend breaks should be
# @usage matdet1<-ec.datadet1(n,befpn,breaks)
#' @author P. Berta, P. Paruolo, S. Verzillo, PG. Lovaglio
#' @description Prepares deterministic data
#' @references Berta et al. 2020
#' @return matdet1 a matrix with the following columns (1_vec, t_vec) and (n+npred) rows
#' @export

ec.datadet1 <- function(n,befpn,breaks){
begtrim <- befpn[1]; endtrim <- befpn[2]           # begtrim, endtrim
nforecast<-befpn[3];npred<-befpn[4];nhstar<-befpn[5] # nforecast, npred, nhstar
nump<-endtrim+npred; ntot<-n+npred                 # nump: total number of predictions
if(anyNA(breaks)==F){nbr1<-length(breaks)+1}else{nbr1<-1} # nbr1: number of breaks plus 1
t<-matrix(0,ntot,nbr1); const <-t;                 # initialize
t[,1]<-matrix((1:ntot),ntot,1)                     # trend
const[,1]<-matrix(1,ntot,1)                        # const
if(nbr1>1){for(i in (2:nbr1)){                     # breaks
    aa<-ifelse(t[,1]>breaks[i-1],1,0); bb<-cumsum(aa)
    const[,i]<-aa; t[,i]<-bb}
colnames(const)<-paste0("const",c(0,breaks))       # 1_vec names
colnames(t)<-paste0("trend",c(0,breaks))}else{     # t_vec names
  colnames(const)<-paste0("const")                 # 1 name
  colnames(t)<-paste0("trend")                     # t name
}
matdet1<-cbind(const,t)
# result:
return(matdet1)
}

