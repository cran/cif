#' listsize number of terms in the search for 1,2,3,4 number of breaks
#'
#' @param myT sample size
#' @param gfill number of gap periods
#' @param start beginning
#' @author P. Paruolo
#' @description computes length-4 vector with number of terms in the search for 1,2,3,4 number of breaks
# @usage sizev<-listsize(myT,gfill,start)
#' @return a vector of 4 elements, with the number of candidate models for 1,2,3,4 breaks
#' @export

listsize <- function(myT,gfill,start)
{vvec<- c(myT-gfill-1,myT-2*(gfill+1),myT-3*(gfill+1),myT-4*(gfill+1))
uvec<- vvec -gfill*rep(1,4)
# ------------ c2i ----------------
c20<- uvec[1]*uvec[2]-0.5*vvec[2]*(vvec[2]+1)+0.5*gfill*(gfill+1)
check <- (myT-3*gfill)*(myT-3*gfill-1)/2
c21<- -uvec[1]+(2*gfill+1)/2; c22<-0.5
c2vec<- c(c20,c21,c22) # c2 vector = (c20,c21,c22)
# ------------ c3i ----------------
c30<- uvec[3]*c2vec[1]+c2vec[2]*vvec[3]*(vvec[3]+1)/2-
c2vec[2]*gfill*(gfill+1)/2 + c2vec[3]*(2*(vvec[3]^3)+
3*(vvec[3]^2)+vvec[3]-gfill-2*(gfill^3)-3*(gfill^2))/6
c31<- -c2vec[1]-0.5*c2vec[2]*(2*gfill+1)-c2vec[3]*(gfill^2)-c2vec[3]*gfill-c2vec[3]/6
c32<- -0.5*c2vec[2]-c2vec[3]*(6*gfill+3)/6
c33<- -c2vec[3]/3
c3vec<- c(c30,c31,c32,c33) # c3 vector = (c30,c31,c32,c33)
# ------------ c2s ----------------
svec<-as.vector(c(1,start,start^2));cvec<-as.vector(c2vec);
c2s<-t(svec) %*% cvec;
# ------------ c3s ----------------
svec<-as.vector(c(1,start,start^2,start^3));
cvec<-as.vector(c3vec);
c3s<-t(svec) %*% cvec;
# ------------ c4s ----------------
h<-start+gfill;
c4s<- c3vec[1]*(vvec[4]-gfill-start)+0.5*c3vec[2]*(vvec[4]*(vvec[4]+1)-h*(h+1))+
  c3vec[3]*((2*vvec[4]+1)*vvec[4]*(vvec[4]+1)-(2*h+1)*h*(h+1))/6+
  c3vec[4]*((vvec[4]^2)*(vvec[4]+1)^2-(h^2)*((h+1)^2))/4
c1s<-vvec[1]-gfill-start
out<-c(c1s,c2s,c3s,c4s)
return(out)}
