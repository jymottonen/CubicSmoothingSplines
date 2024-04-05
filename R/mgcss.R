#' Testing of spline growth curve model
#'
#' \code{mgcss} is used to test the mean growth curves using cubic smoothing splines. 
#'
#' @param Y a qxn matrix of response vectors. The ith column contains the q-variate response vector
#' of the ith individual. The jth row contains the responses of the individuals in the jth time point.
#' @param A an nxm between-individual design matrix. You should give either \code{A} or \code{grp}, not both.
#' @param grp an n-dimensional vector giving the m treatment groups of the observations. 
#' The nxm between-individual design matrix \code{A} is constructed using the values of this vector:
#' The ith row of \code{A} is a unit vector with 1 in position \code{grp[i]} and zeros elsewhere.
#' @param model covariance structure. The default "unif" assumes the uniform covariance structure
#' \eqn{R = d^2 1_q1_q'+I_q}. 
#' @details 
#' Here are the details of the function...
#' @return A list containing the following components:
#' \describe{
#' \item{F.value}{the observed value of the F test statistic.}
#' \item{p.value}{approximate p-value.}
#' \item{p.perm}{permutation test p-value.}
#' \item{alpha}{the value of alpha.}
#' \item{c}{the value of c.}
#' \item{Ghat}{the spline fit.}
#' \item{Gtilde}{the approximated spline fit.}
#' }
#' @importFrom stats pf
#' @export
mgcss<-function(Y,A=NULL,grp=NULL,model="unif")
{
  if(is.null(A)&is.null(grp))
    stop("You should give A or grp")
  if(!is.null(A)&!is.null(grp))
    stop("You should give either A or grp, not both")
  if(!is.null(grp))
  {
    groups<-unique(grp)
    m<-length(groups)
    n<-ncol(Y)
    A<-matrix(0,n,m)
    for(j in 1:m)
    {
      A[grp==groups[j],j]<-1
    }
  }
  
  q<-nrow(Y)
  n<-ncol(Y)
  m<-ncol(A)
  tt<-1:q
  K<-roughness(tt)
  
  y<-c(Y)
  Pa<-A%*%solve(t(A)%*%A)%*%t(A)
  alphat<-seq(0.01,0.5,length=100)
  GCV<-NULL
  for(i in 1:length(alphat))
  {
    alpha<-alphat[i]
    S<-solve(diag(q)+alpha*K)
    edf<-sum(diag(S))
    yhat<-kronecker(Pa,S)%*%y
    GCV[i]<-mean((y-yhat)^2)/(1-m*edf/(n*q))^2
  }
  alpha<-alphat[which.min(GCV)]
  print(alpha)
  
  gcv.res<-gcv.alpha(Y,tt,A,alpha.min=0.01,alpha.max=0.5,len=100)
  print(gcv.res$alpha.hat)
  
  
  
  S<-solve(diag(q)+alpha*K)
  Ghat<-S%*%Y%*%A%*%solve(t(A)%*%A)
  
  M<-eigen(S)$vectors
  eigen(S)$values
  m1<-rep(1,q)/sqrt(q)
  St<-sqrt(sum((tt-mean(tt))^2))
  m2<-(tt-mean(tt)*rep(1,q))/St
  M[,1]<-m1
  M[,2]<-m2
  
  cct<-1:12
  GCV2<-NULL
  for(i in 1:length(cct))
  {
    cc<-cct[i]
    Pm<-M[,1:cc]%*%t(M[,1:cc])
    yhat<-kronecker(Pa,Pm)%*%y
    GCV2[i]<-mean((y-yhat)^2)/(1-m*cc/(n*q))^2
  }
  cc<-cct[which.min(GCV2)]
  
  Mstar<-M[,1:cc]
  Pm<-Mstar%*%t(Mstar)
  Gtilde<-Pm%*%Y%*%A%*%solve(t(A)%*%A)
  
  Omegahat<-t(Mstar)%*%Y%*%A%*%solve(t(A)%*%A)
  C<-cbind(rep(0,cc-1),diag(cc-1))
  D<-matrix(c(1,-1),2,1)
  nu<-nrow(C)
  g<-ncol(D)
  COD<-C%*%Omegahat%*%D
  Qstar<-sum(diag(COD%*%solve(t(D)%*%solve(t(A)%*%A)%*%D)%*%t(COD)))
  sigmahat2<-sum(diag(t(Y)%*%(diag(q)-Pm)%*%Y))/(n*(q-cc)); sigmahat2
  FF<-(Qstar/(nu*g))/sigmahat2
  FF<-c(FF)
  p.value<-1-pf(FF,nu*g,n*(q-cc))
  
  nperm<-100000
  FF2<-NULL
  for(i in 1:nperm)
  {
    smp<-sample(1:n,n)
    A2<-A[smp,]
    Omegahat2<-t(Mstar)%*%Y%*%A2%*%solve(t(A2)%*%A2)
    COD2<-C%*%Omegahat2%*%D
    Qstar2<-sum(diag(COD2%*%solve(t(D)%*%solve(t(A2)%*%A2)%*%D)%*%t(COD2)))
    FF2[i]<-(Qstar2/(nu*g))/sigmahat2
  }
  FF2
  p.perm<-mean(FF2>=FF)
  res<-list(F.value=FF,p.value=p.value,p.perm=p.perm,alpha=alpha,c=cc,Ghat=Ghat,Gtilde=Gtilde)
  class(res) <- "mgcss"
  return(res)
}
