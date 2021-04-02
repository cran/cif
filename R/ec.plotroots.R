#' Companion matrix of the VAR
#'
#' @param roots are the roots of the companion matrix, see ec.companion.R
# @usage  ec.plotroots(roots)
#' @author P. Berta, P. Paruolo, S. Verzillo, PG. Lovaglio
#' @description plots roots and the unit circle
#' @importFrom graphics plot lines
#' @references Berta et al. 2020
#' @return does not return output, just creates a graph
#' @export

ec.plotroots <- function(roots)
{ val<-c(seq(from = 0.02, to = 0.7, by =0.02),seq(from = 0.71, to = 1, by =0.01))
  x<-c(rev(-val),0,val)
  z<-roots$values
  plot(z,ylim=c(-1,1),xlim=c(-1,1),asp=1,col="blue")
  lines(x,y=(1-x^2)^.5,type="l",col="red")
  lines(x,y=-(1-x^2)^.5,type="l",col="red")
}
