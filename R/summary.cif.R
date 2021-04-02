#' summary function for cif
#'
#' @param object is the name of the cif object created by cif
#' @param ... other parameters
#' @param digits integer indicating the number of decimal places (round) or significant digits (signif) to be used.
#' @author P. Berta, P. Paruolo, S. Verzillo, PG. Lovaglio
# @usage summary(object,...,digits)
#' @description Summary function for presize
#' @references Berta et al. 2020
#' @return returns summary output from model estimation and forecasting
#' @export

summary.cif <- function(object, ... , digits=4) {
  stopifnot(inherits(object, "cif"))
  # ----------- stars -------------------
  pval<-object$acoef_pval; myd<-dim(pval); nreg<-myd[2]; neq<-myd[1]
  acoef<-object$acoef; acoef_se<-object$acoef_se; acoef_t<-object$acoef_t

  symbols=c("***", "**", "*", ""); cutpoints = c(0, 0.01, 0.05, 0.1, 1)
  vecp <- as.vector(t(pval)); mystars <- cut(vecp,breaks=cutpoints,labels=symbols)

  # --- Eq Cointegration -
  EqC <- round((object$vEG1step),digits=digits)
  cat("\t\n", sprintf("Cointegrating Equation:\n"))
  print(EqC)
  # ----------- equations -------------------
  for(i in (1:neq)){
  otab<-round(cbind(acoef[i,],acoef_se[i,],acoef_t[i,],pval[i,]),digits=digits)
  regOutput <- as.data.frame(cbind(otab,as.character(mystars[((i-1)*nreg+1):(i*nreg)])))
  colnames(regOutput) <- c("Estimate", "Std. Error", "T-Value", "P-Value", "Stars")
  cat("\t\n", sprintf("Results for Equation"),i,":\n")
  print(regOutput)
  }
  # ----------- Wald stat -------------------
  cat("\t\n", sprintf("Joint significance of predictors:\n"))
  mywald<-t(object$wald)
  mystars <- cut(mywald[,2],breaks=cutpoints,labels=symbols)
  WOutput <- as.data.frame(cbind(round(mywald,digits=2),as.character(mystars)))
  colnames(WOutput)[3]<-"Stars"
  print(WOutput,quote=FALSE)
  # print(t(WOutput),quote=FALSE)
  # ----------- Forecast -------------------
  cat("\t\n", sprintf("Indices of forecast accuracy in equation 1:\n"))
  print(round(cbind(object$tablefap),digits=digits))

  NextMethod("summary")
}
