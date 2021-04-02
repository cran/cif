#' lag j of matrix or vector y
#'
#' @param y column vector or matrix
#' @param j number of lags
#' @param fill value to be used to fill the missing values at the beginning, default = NA
#' @author P. Berta, P. Paruolo, S. Verzillo, PG. Lovaglio
# @usage x <- lagn(y,j,fill=NA)
#' @description lagn(y,j,fill=NA) produces lag j of matrix or vector y, with fill in missing j cells
#'              REM: alternative to "lead-lag" {dplyr} which applies to vector y
#' @references Berta et al. 2020
#' @return y lagged j cells, with fill in the missing j positions
#' @export

lagn <- function(y,j,fill=NA)
{ my <- as.matrix(y);
  a <- dim(my); n <-a[1]; p<-a[2];
  ylag <- rbind(matrix(fill,j,p),as.matrix(my[1:(n-j),]))
  rownames(ylag)<-rownames(y)
  myname <- colnames(y)
  if(is.null(myname)){colnames(ylag)<- paste("var",1:p,"_",j,sep="")}else
  {colnames(ylag)<-paste(myname,j,sep="_")}
  # result:
  return(ylag)
}
