#' Testing of multivariate spline growth curve model
#'
#' \code{mgcss} is used to test the mean growth curves using multivariate cubic smoothing splines. 
#'
#' @param Y a \eqn{qs\times n} matrix of the vectors of measurements of \eqn{s} responses. 
#' @param s number of responses.
#' @param A an \eqn{n\times m} between-individual design matrix. You should give either \eqn{A} or \eqn{grp}, not both.
#' @param grp an \eqn{n}-dimensional vector giving the \eqn{m} treatment groups of the observations. 
#' The \eqn{n\times m} between-individual design matrix \eqn{A} is constructed using the values of this vector:
#' The \eqn{i}th row of \eqn{A} is a unit vector with 1 in position \eqn{grp[i]} and zeros elsewhere.
#' @param t a \eqn{q}-dimensional vector of the time points.
#' @param alpha an \eqn{s}-dimensional vector of fixed values of \eqn{\alpha}. If it is not given, \eqn{\alpha}s are estimated using gcv criteria. 
#' The minimum of gcv criteria is found using grid search.
#' @param alpha.min the lower bound of grid of \eqn{\alpha} values when estimating \eqn{\alpha} using gcv criteria.
#' @param alpha.max the upper bound of grid of \eqn{\alpha} values when estimating \eqn{\alpha} using gcv criteria.
#' @param model covariance structure. The default "unif" assumes the multivariate version of the 
#' uniform covariance structure \eqn{R = (I_s\otimes 1_q)D(I_s\otimes 1_q)'+I_{qs}}. 
#' @details 
#' Here are the details of the function...
#' @return A list containing the following components:
#' \describe{
#' \item{F.value}{the observed value of the F test statistic.}
#' \item{df1}{the degrees of freedom of the numerator of the F test statistic.}
#' \item{df2}{the degrees of freedom of the denominator of the F test statistic.}
#' \item{p.value}{approximate p-value.}
#' \item{sigmahat2}{estimate of \eqn{\sigma^2}.}
#' \item{alpha}{the value of \eqn{\alpha}.}
#' \item{c}{the value of the number of eigenvectors \eqn{c}.}
#' \item{Gtilde}{the approximated spline fit.}
#' }
#' @importFrom stats pf
#' @export
mgcss<-function(Y, s, A=NULL, grp=NULL, t=1:(dim(Y)[1]/s), 
                alpha=NULL, alpha.min=0.1, alpha.max=1000, model="unif")
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
  
  qs<-nrow(Y)
  q<-qs/s
  n<-ncol(Y)
  m<-ncol(A)
  K<-roughness(t)
  
  #Estimation of the smoothing parameters
  S<-vector("list",s)
  if(is.null(alpha))
  {
    for(j in 1:s)
    {
      Yj<-Y[((j-1)*q+1):(j*q),]
      gcv.res<-gcv.alpha(Yj,t,A,alpha.min=alpha.min,alpha.max=alpha.max,len=100)
      if(gcv.res$minpoint==(-1))
        cat("The minimum gcv of response",j,"was found at the point alpha.min=",alpha.min,"\n")
      if(gcv.res$minpoint==1)
        cat("The minimum gcv of response",j,"was found at the point alpha.max=",alpha.max,"\n")
      alpha[j]<-gcv.res$alpha.hat
      S[[j]]<-solve(diag(q)+alpha[j]*K)
    }
  }
  else
  {
    for(j in 1:s)
    {
      Yj<-Y[((j-1)*q+1):(j*q),]
      S[[j]]<-solve(diag(q)+alpha[j]*K)
    }
  }
  
  W<-diag(alpha)
  Ghat<-solve(diag(q*s)+kronecker(W,K))%*%Y%*%A%*%solve(t(A)%*%A)
  
  #Estimation of the dimensions
  P<-vector("list",s)
  Gtilde<-vector("list",s)
  PY<-NULL
  MY<-NULL
  c<-NULL
  sigmahat2<-0
  for(j in 1:s)
  {
    M<-eigen(S[[j]])$vectors
    m1<-rep(1,q)/sqrt(q)
    St<-sqrt(sum((t-mean(t))^2))
    m2<-(t-mean(t)*rep(1,q))/St
    M[,1]<-m1
    M[,2]<-m2
    Yj<-Y[((j-1)*q+1):(j*q),]
    c[j]<-gcv2.dim(Yj,A,M)$c
    P[[j]]<-M[,1:c[j]]%*%t(M[,1:c[j]])
    Gtilde[[j]]<-P[[j]]%*%Yj%*%A%*%solve(t(A)%*%A)
    PY<-rbind(PY,P[[j]]%*%Yj)
    MY<-rbind(MY,t(M[,1:c[j]])%*%Yj)
    sigmahat2<-sigmahat2+sum(diag(t(Yj)%*%(diag(q)-P[[j]])%*%Yj))*(1/(n*(q-c[j])))
  }
  
  Omegahat<-MY%*%A%*%solve(t(A)%*%A)
  
  ctot<-sum(c)
  D<-rbind(1,-diag(m-1))
  C.Omega<-NULL
  cb<-0
  for(j in 1:s)
  {
    ca<-cb+1; cb<-cb+c[j]
    C.Omega<-rbind(C.Omega,cbind(0,diag(c[j]-1))%*%Omegahat[ca:cb,])
  }
  COD<-C.Omega%*%D
  
  Qstar<-sum(diag(COD%*%solve(t(D)%*%solve(t(A)%*%A)%*%D)%*%t(COD)))
  
  df1<-(ctot-s)*(m-1)
  df2<-n*(s*q-ctot)
  
  FF<-(Qstar/df1)/sigmahat2
  FF<-c(FF)

  p.value<-1-pf(FF,df1,df2)
  
  res<-list(F.value=FF,df1=df1,df2=df2,p.value=p.value,sigmahat2=sigmahat2,alpha=alpha,c=c,Gtilde=Gtilde)
  class(res) <- "mgcss"
  return(res)
}
