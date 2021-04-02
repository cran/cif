#' Multivariate Least-Squares regression
#'
#' @param y left hand side data matrix (one or more columns)
#' @param x right hand side data matrix (one or more columns)
#' @param df_flag flag = TRUE for degrees of freedom correction for the variance
#' @author P. Berta, P. Paruolo, S. Verzillo, PG. Lovaglio
#' @description Multivariate Least-Squares regression y = x beta + u
# @usage out<-mls(y, x, df_flag=FALSE)
#' @importFrom stats pchisq
#' @references Berta et al. 2020
#' @return out regression coefficients and related statistics
#' @export

mls <- function(y, x, df_flag=FALSE)
{ # Multivariate Least-Squares regression
  # y     =  x      \beta   +  u
  # (n x g) (n x k) (k x g) (n x g)
  # returns \hat{\beta}  = (x'x)^{-1}x'y and related stats
  a <- dim(as.matrix(x))
  n <- a[1]                                     # number of obs n
  k <- a[2]                                     # number of regressors
  if(df_flag==TRUE){df <- n-k}else{df <- n}     # degrees of freedom
  qx <- qr(x)                                   # QR decomposition
  beta <- solve.qr(qx, y)                       # betahat
  fit <- qr.fitted(qx,y)                        # fitted values
  xtxinv <- chol2inv(qx$qr)                     # (X'X)^{-1}
  residuals <- y - fit                          # residuals
  Omega <- (t(residuals) %*% residuals)/df      # estimate of Omega = E(u_i u_i')
  semat <- (diag(xtxinv) %*% t(diag(Omega)))^.5 # matrix of standard errors
  tmat <- beta / semat                          # matrix of t values
  pmat <- pchisq(tmat^2,1,lower.tail = FALSE)   # matrix of p values
  # write output
  out <- list(); out$n <- n; out$k <- k; out$df <- df
  out$coeffs <- beta; out$fit <- fit; out$residuals <- residuals
  out$XXinv <- xtxinv; out$Omega <- Omega; out$semat <- semat
  out$tmat <- tmat; out$pmat <- pmat
  # result:
  return(out)
}
