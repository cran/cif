#'plots forecasts
#'
#' @param y actual values and forecasts (point forecast, lower bound, upper bound)
#' @param x time calendar
#' @param lcolact color actual (scalar)
#' @param lcolfor color forecasts
#' @param ltypefor type forecasts
#' @param polycol color polygons
#'              if one wishes to have different lcolfor,ltypefor,polycol
#'              by week > make linecol, linetype, polycol vectors, indexed by week
#' @param myylim vector with min and max for y axis
#' @param ...  other plot parameters
#' @author P. Paruolo
#' @description plot actual and forecast intervals
# @usage myplotfor(y,x,lcolact,lcolfor,ltypefor,polycol,...)
#' @references Berta et al. 2020
#' @importFrom graphics plot polygon lines
#' @return does not return output, just creates a graph
#' @export

ec.plotfor<-function(y,x=NA,lcolact=NA,lcolfor=NA,
                    ltypefor=NA,polycol=NA,myylim=NA,...){
  p<-ncol(as.matrix(y)); n<-nrow(as.matrix(y))
  act<-y[,1]; yfor<-y[,-1]; nweeks<-(p-1)/3;                      # actual, forecasts, number of weeks
  if(anyNA(x)){x<-1:n} # if x missing, create sequence
  if(anyNA(myylim)==T){
  myylim <-c(min(min(y,na.rm = TRUE)),max(max(y,na.rm = TRUE)))}  # range of plot
  if(length(lcolfor)==1){lcolfor<-matrix(lcolfor,nweeks,1)}       # duplicate lcolfor
  if(length(ltypefor)==1){ltypefor<-matrix(ltypefor,nweeks,1)}    # duplicate ltypefor
  if(length(polycol)==1){polycol<-matrix(polycol,nweeks,1)}       # duplicate polycol
  plot(x,act,type="l",col=lcolact,ylim=myylim,lty=1, ...)         # graph actual once
  for(i in (1:nweeks)){                                           # cycle over weeks
    icol<-(i-1)*3;  iweek <- !is.na(yfor[,icol+1])                # choose days of forecast
    x1<-c(x[iweek],rev(x[iweek]))                                 # x coord polygon
    y1<-c(yfor[iweek,icol+2],rev(yfor[iweek,icol+3]))             # y coord polygon
    polygon(x1, y1, col=polycol[i],border = NA)                   # graph polygon
    lines(x,yfor[,icol+1],type="l",lty=ltypefor[i],col=lcolfor[i])# graph point forecast
  }
  lines(x,act,type="l",col=lcolact,ylim=myylim,lty=1)             # graph actual again
}



