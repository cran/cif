#' Computes Indices of Forecast Performance
#'
#' @param afdlin actual + forecast values + fcse
#' @param rwsigma standard deviation of Random Walk in sample
#' @param rwabsmean mean absolute deviation of Random Walk in sample
#' @param kval how many se to use, default kval = 1.959964
#' @author P. Berta, P. Paruolo, S. Verzillo, PG. Lovaglio
#' @description indices of forecast performance
#' @references Berta et al. 2020
#' @importFrom stats pnorm
# @usage ec.ifp <- function(afdlin,rwsigma,rwabsmean,kval)
#' @return list of indices of forecast performance
#'  1: index for model forecast
#'  0: index for Random Walk forecast
#' @export

ec.ifp <- function(afdlin,rwsigma,rwabsmean,kval=1.959964)
{# -------------  indices for point forecasts 1st variable -----------------------
p<-length(rwsigma) # p
chebyalpha <- 1/(kval^2)                                           # alpha using Chebyshev
gaussalpha <- 2*(1-pnorm(kval,0,1))                                # alpha using Gaussian
table<-matrix(0,14,3*p)                                             # initialise table
mycolnames <- vector(mode='character',3*p)
for(ivar in (1:p)){                                                 # build table for all equations
fe1<-afdlin[-1,ivar]-afdlin[-1,(p+ivar)]
afe1<-abs(fe1); mae1<-mean(afe1)                                    # MAE1
fe0<-afdlin[-1,ivar]-afdlin[1,ivar]
afe0<-abs(fe0); mae0<-mean(afe0); baaser<-maser<-maer<-mae1/mae0    # MAE0, MAER, MASER, BAASER
fee1<-exp(afdlin[-1,ivar])-exp(afdlin[-1,(p+ivar)])
afee1<-abs(fee1);maee1<-mean(afee1)                                 # MAEexp1
fee0<-exp(afdlin[-1,ivar])-exp(afdlin[1,ivar])
afee0<-abs(fee0);maee0<-mean(afee0); maeer <- maee1/maee0           # MAEexp0, MAEexpR
smapeden1 <- 0.5*(abs(afdlin[-1,ivar])+abs(afdlin[-1,(p+ivar)]))    # den smape1
smape1 <- mean(afe1/smapeden1)                                      # sMAPE1
smapeden0 <- 0.5*(abs(afdlin[-1,ivar])+abs(afdlin[1,ivar]))         # den smape0
smape0 <- mean(afe0/smapeden0)                                      # sMAPE0
smaper <- smape1/smape0                                             # sMAPER
rfe1 <- afe1/abs(afdlin[-1,ivar]); mape1 <- mean(rfe1)              # MAPE1
rfe0 <- afe0/abs(afdlin[1,ivar]);  mape0 <- mean(rfe0)              # MAPE0
maper <- mape1/mape0                                                # MAPER
rfee1 <- afee1/abs(exp(afdlin[-1,ivar])); mapee1 <- mean(rfee1)     # MAPEexp1
rfee0 <- afee0/abs(exp(afdlin[1,ivar]));  mapee0 <- mean(rfee0)     # MAPEexp0
mapeer <- mapee1/mapee0                                             # MAPEexpR
maape1 <- mean(atan(rfe1)); maape0 <- mean(atan(rfe0))              # MAAPE1, MAAPE0
maaper <- maape1/maape0                                             # MAAPER
baase1 <- mean(afe1/afe0)                                           # BAASE (may be +Inf)
msisden <- as.numeric(rwabsmean[ivar]);
mase0 <- mae0/msisden; mase1 <- mae1/msisden;                       # MASE0, MASE1
# -------------  scoring for interval forecasts ----------------------
obs <- afdlin[-1,ivar]; n1 <- length(obs)              # observed
    #--- model i scores ----
l1 <- afdlin[-1,(p+ivar)]-kval*afdlin[-1,(2*p+ivar)]  # lower point of forecast int
u1 <- afdlin[-1,(p+ivar)]+kval*afdlin[-1,(2*p+ivar)]  # upper point of forecast int
aa1 <- ifelse(obs<l1,l1-obs,0)+ifelse(obs>u1,obs-u1,0) # penalty
s1c <- u1-l1 + (kval/chebyalpha)*aa1                  # score using Chebyshev alpha
s1g <- u1-l1 + (kval/gaussalpha)*aa1                  # score using Gaussian alpha
    #--- RW scores ----
rwfse <- ((1:n1)^.5)*rwsigma[ivar]                     # interval 0 for RW model
l0 <- afdlin[1,ivar]-kval*rwfse                       # lower point of forecast int
u0 <- afdlin[1,ivar]+kval*rwfse                       # upper point of forecast int
aa0 <- ifelse(obs<l0,l0-obs,0)+ifelse(obs>u0,obs-u0,0) # penalty
s0c <- u0-l0 + (kval/chebyalpha)*aa0                  # score using Chebyshev alpha
s0g <- u0-l0 + (kval/gaussalpha)*aa0                  # score using Gaussian alpha
# -------------  indices for interval forecasts ----------------------
aisc1<-mean(s1c);aisc0<-mean(s0c);msiscr<-aiscr<-aisc1/aisc0    # AIS and AISR Chebyshev MSIScR
aisg1<-mean(s1g);aisg0<-mean(s0g);msisgr<-aisgr<-aisg1/aisg0    # AIS and AISR Gaussian MSISgR
msisc1 <- aisc1/msisden; msisg1 <- aisg1/msisden;               # MSIS1
msisc0 <- aisc0/msisden; msisg0 <- aisg0/msisden;               # MSIS0
basisc1 <- mean(s1c/s0c); basisg1 <- mean(s1g/s0g);             # BASIS
# -------------  build table ----------------------
mae<-c(mae1,mae0,maer)
maee<-c(maee1,maee0,maeer)
smape<-c(smape1,smape0,smaper)
mape<-c(mape1,mape0,maper)
mapee<-c(mapee1,mapee0,mapeer)
maape<-c(maape1,maape0,maaper)
baase<-c(baase1,1,baaser)
mase<-c(mase1,mase0,maser)
aisc<-c(aisc1,aisc0,aiscr)
aisg<-c(aisg1,aisg0,aisgr)
msisc<-c(msisc1,msisc0,msiscr)
msisg<-c(msisg1,msisg0,msisgr)
basisc<-c(basisc1,1,basisc1)
basisg<-c(basisg1,1,basisg1)
table1<-rbind(mae,maee,smape,mape,mapee,maape,baase,mase,aisc,aisg,msisc,msisg,basisc,basisg)
isel<-((ivar-1)*3+1):((ivar-1)*3+3)
mycolnames[isel]<-paste0("eq",ivar,c(":model",":RW",":ratio"))
table[,isel]<-table1
}
indnames <-c("MAE","MAEexp","sMAPE","MAPE","MAPEexp","MAAPE","BAASE","MASE","AISc","AISg",
             "MSISc","MSISg","BASISc","BASISg")
rownames(table)<-indnames; colnames(table)<-mycolnames
# result:
return(table)
}

