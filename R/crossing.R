#' computes at which observation a vector y crosses ref for the first time
#'
#' @param yfor yfor is either a vector and a matrix
#' @param ref ref is the refence value
#' @author P. Berta, P. Paruolo, S. Verzillo, PG. Lovaglio
#' @description Computes at which observation vector yfor crosses yref for the first time if it is not crossed, then 0 is returned
# @usage whensign<-crossing(yfor, ref=0)
#' @references Berta et al. 2020
#' @return whensign, a matrix with observation number at which there is crossing
#' @export

crossing <- function(yfor, ref=0)
{ # computes at which observation
  # vector yfor crosses ref for the first time yfor<-afdd[,p+1]
  # if it is not crossed, then ref<-0 is returned
  nobs <- length(yfor); obseq <- 1:(nobs-1)            # index of observations
  cdif <- sign(yfor-ref); changesign <- diff(cdif)     # change in signs
  vind <- changesign != 0; vdate <- obseq[vind]+1      # obs num of sign changes
  n1 <- length(vdate)+1; whensign <- matrix(0,3,n1+1)  # output by rows
  whensign[1:2,1] <- c(1,cdif[1]);
  if(n1>1){
    whensign[1,2:n1] <- vdate
    whensign[2,2:n1] <- cdif[vdate]
    mcnames<-c("t",paste0("change",1:(n1-1)),"last")
  }else{mcnames<-c("t","last")}
  whensign[1,n1+1] <- nobs
  whensign[2,n1+1] <- cdif[nobs]
  whensign[3,1:n1] <- diff(whensign[1,])
  rownames(whensign)<- c("position","sign","duration")
  colnames(whensign)<- mcnames
# result:
return(whensign)
}
