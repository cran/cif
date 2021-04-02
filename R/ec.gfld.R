#' ec.gfld plots forecasts of levels and difference with confidence bars
#'
#' @param obj output of ec.main
#' @param whichseries series number
#' @param nsigma how many standard deviations in confidence bars
#' @param jointboo boolean: TRUE if 1x2 graph, FALSE otherwise
#' @param epsboo boolean: TRUE eps graph, FALSE pdf graph
#' @param filename string, name of the file (no extension)
#' @param xvec vector of dates to place on x axis
#' @param yvec vector of exp(y) values to display on y axis
#' @param cal calendar vector
#' @param lar length of arrows in error bars
#' @param ... other plot parameters
#' @author P. Paruolo
# @usage ec.gfld(obj, whichseries=1, nsigma=3,jointboo=TRUE, epsboo=TRUE,filename="whatever",xvec,yvec,cal,lar=0.025,...)
#' @description plots forecasts of levels and difference with confidence bars
#' @importFrom grDevices dev.off pdf postscript setEPS
#' @importFrom graphics par
#' @return does not return output, just creates a double graph
#' @export

ec.gfld <- function(
  obj, whichseries=1, nsigma=3,jointboo=TRUE, epsboo=TRUE,
  filename="whatever",xvec,yvec,cal,lar=0.025,...)
{ # --------------- begin plot ----------------
  if(jointboo==TRUE){mywidth<-2*5.83}else{mywidth<-5.83}       # graph width
  if(epsboo==TRUE){myfilename<-paste0(filename,".eps");        # filename eps
    setEPS(); postscript(file=myfilename,width=mywidth, height=4.13) # graph begin eps
  }else{
    myfilename<-paste0(filename,".pdf")                        # filename pdf
    pdf(file=myfilename,width=mywidth, height=4.13)}           # graph begin pdf
  if(jointboo==TRUE){par(mfcol=c(1,2))}                        # split graph if joint
  # ------------------- plot forecast levels ---------------------
  ec.gfl(obj,whichseries,nsigma,xvec,yvec,cal,lar,...)
  # ------------------- end plot forecast levels -----------------
  if(jointboo==FALSE){dev.off();myfilename2<-paste0("D",myfilename);
    if(epsboo==TRUE){setEPS();
      postscript(file=myfilename2,width=mywidth, height=4.13)
    }else{
      pdf(file=myfilename2,width=mywidth, height=4.13)}}
  # ------------------- plot forecast differences -----------------
  ec.gfd(obj,whichseries,nsigma,xvec,yvec,cal,lar,...)
  # ------------------- end plot forecast differences -------------
  dev.off()
}

