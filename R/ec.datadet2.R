#' prepares deterministics D^(2)
#'
#' @param det1 is the det term with constant and trend created by ec.datadet1.R
#' @param booseas is a boolean for daily seasonal dummies
#' @param pntdates is a vector of integers where the point dummies should be
# @usage det2mat <- ec.datadet2(det1,booseas,pntdates)
#' @author P. Berta, P. Paruolo, S. Verzillo, PG. Lovaglio
#' @description Prepares deterministic dummies for de-meaned daily seasonal and difference point dummies
#' @references Berta et al. 2020
#' @return det2mat a matrix with the following columns (daily_seas, point_dummies) and n+npred rows
#' @export

ec.datadet2 <- function(det1,booseas=NA,pntdates=NA)  {
  ntot<-nrow(det1)                              # number of rows
  if(booseas==T){nseas<-6}else{nseas<-0}        # number of daily seasonals
  boodates<-anyNA(pntdates)==F                  # boolean for point dummies
  if(booseas==F&&boodates==F){matdet2<-NA;
    return(matdet2);stop("no seas or point dummies")} # check
  if(boodates==T){ndates<-length(pntdates)      # number of point dummies
       }else{ndates<-0}
  matdet2<-matrix(0,ntot,nseas+ndates)          # empty matdet2
  # ------------- daily seasonal dummies ---------------
  if(booseas==T){sixo7<-6/7; moneo7<- -1/7;
    myvec<-c(sixo7,moneo7,moneo7,moneo7,moneo7,moneo7,moneo7)
    mymult<-ceiling(ntot/7); mycol<-matrix(myvec,mymult*7,1)
    matdet2[,1]<-mycol[1:ntot];                 # first daily seasonal
  for(i in (2:6)){matdet2[,i]<-lagn(matdet2[,1],i-1,moneo7)}# lags of first daily seasonal
  }
  # ------------- point dummies ---------------
  if(boodates==T){for(i in (1:ndates)){         # <--- point dummies
    matdet2[pntdates[i],nseas+i]<-1;            # assigns 1
    matdet2[pntdates[i]+1,nseas+i]<--1}         # assigns -1
  }
  # ------------- colnames ---------------
  if(booseas==T&&boodates==T){mynames<-c(paste0("seas",(1:6)),paste0("day",pntdates))}
  if(booseas==T&&boodates==F){mynames<-paste0("seas",(1:6))}
  if(booseas==F&&boodates==T){mynames<-paste0("day",pntdates)}
  colnames(matdet2)<-mynames
# result:
return(matdet2)
}

