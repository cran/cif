#' appends NA at beginning of diff(y)
#'
#' @param y either a vector and a matrix
#' @author P. Paruolo
#' @description appends NA at beginning of diff(y) and creates column names accordingly when y is either a vector and a matrix
# @usage Dy<-diffe(y)
#' @references Berta et al. 2020
#' @return Dy contains the differences of y, with NA appended at the start
#' @export

diffe <- function(y)
{ y<-as.matrix(y); p <- ncol(y)                  # p is number of columns
  Dy <- rbind(matrix(NA,nrow=1,ncol=p),diff(y))  # creates diff
  rownames(Dy)<-rownames(y)
  myname <- colnames(y)                          #  names
  if(is.null(myname)){colnames(Dy)<- paste("Dvar",1:p,sep="")}else
  {colnames(Dy)<-paste("D",myname,sep="")}
  return(Dy)
}
