#' Testing of cubic smoothing splines and semi-parametric regression models
#'
#' \code{css} is used to test the cubic smoothing splines and semi-parametric regression models.
#'
#' @param y an \eqn{n}-dimensional vector of responses. 
#' @param x an \eqn{n}-dimensional vector of measuring points. 
#' @param U a full rank \eqn{n\times k} matrix of values of \eqn{k} explanatory variables. 
#' If the matrix \eqn{U} is given, the tests for cubic smoothing spline (test 1), semi-parametric model (test 2) and linear model (test 3)
#' are performed.  If the matrix \eqn{U} is not given, only the test for cubic 
#' smoothing spline (test 1) is performed.
#' @details 
#' Here are the details of the function...
#' @return A list containing the following components:
#' \describe{
#' \item{y}{\eqn{n}-dimensional vector of responses.}
#' \item{x}{\eqn{n}-dimensional vector of measuring points.}
#' \item{alpha}{the value of \eqn{\alpha}.}
#' \item{c}{estimated number of eigenvectors.}
#' \item{gcv}{values of the gcv criteria.}
#' \item{ghat}{estimate of cubic smoothing spline.}
#' \item{F.value1}{observed value of the F test statistic (Test 1).}
#' \item{df1.1}{degrees of freedom of the numerator of the F test statistic (Test 1).}
#' \item{df2.1}{degrees of freedom of the denominator of the F test statistic (Test 1).}
#' \item{p.value1}{p-value (Test 1).}
#' \item{F.value2}{observed value of the F test statistic (Test 2).}
#' \item{df1.2}{degrees of freedom of the numerator of the F test statistic (Test 2).}
#' \item{df2.2}{degrees of freedom of the denominator of the F test statistic (Test 2).}
#' \item{p.value2}{p-value (Test 2).}
#' \item{F.value3}{observed value of the F test statistic (Test 3).}
#' \item{df1.3}{degrees of freedom of the numerator of the F test statistic (Test 3).}
#' \item{df2.3}{degrees of freedom of the denominator of the F test statistic (Test 3).}
#' \item{p.value3}{p-value (Test 3).}
#' }
#' @importFrom stats pf
#' @export
css<-function(y, x, U=NULL)
{
  y<-c(y)
  x<-c(x)
  if(!is.null(U)&is.vector(U))
    U<-matrix(U,ncol=1)
  
  #The number of observations
  n <- length(y)
  
  #Incidence matrix N
  N<-1*outer(unique(x),x, "==")
  N<-t(N)
  
  #Roughness matrix K
  K<-roughness(unique(x))
  
  #Smoother matrix S (Definition (25))
  alpha<-1
  S<-N%*%solve(t(N)%*%N+alpha*K)%*%t(N)
  
  #Estimate the number of eigenvectors using gcv
  Tn<-eigen(S)$vectors
  m1<-rep(1,n)/sqrt(n)
  St<-sqrt(sum((x-mean(x))^2))
  m2<-(x-mean(x)*rep(1,n))/St
  Tn[,1]<-m1
  Tn[,2]<-m2
  res.gcv<-gcv1.dim(y, Tn)
  c<-res.gcv$c
  
  #Test 1
  T2<-Tn[,1:2]
  Tc<-Tn[,1:c]
  ghat<-(Tc%*%t(Tc))%*%y #Cubic spline curve estimate
  sslin<-y%*%(diag(n)-T2%*%t(T2))%*%y
  ssspl<-y%*%(diag(n)-Tc%*%t(Tc))%*%y
  df1.1<-c-2
  df2.1<-n-c
  F.value1<-(sslin-ssspl)/df1.1
  F.value1<-F.value1/(ssspl/df2.1)  
  p.value1<-1-pf(F.value1,df1.1,df2.1)
  
  #Test 2
  if(!is.null(U))
  {
    k<-ncol(U)
    Uh<-(diag(n)-Tc%*%t(Tc))%*%U
    U2<-cbind(1,x,U)
    P2<-U2%*%solve(t(U2)%*%U2)%*%t(U2)
    M3<-Uh%*%solve(t(Uh)%*%Uh)%*%t(Uh)
    ssmin<-y%*%(diag(n)-(Tc%*%t(Tc)+M3))%*%y
    ssreg2<-y%*%(diag(n)-P2)%*%y
    df1.2<-c-2
    df2.2<-n-k-c
    F.value2<-(ssreg2-ssmin)/df1.2
    F.value2<-F.value2/(ssmin/df2.2)
    p.value2<-1-pf(F.value2,df1.2,df2.2)
  }
  else
  {
    df1.2<-NA
    df2.2<-NA
    F.value2<-NA
    p.value2<-NA
  }

  #Test 3
  if(!is.null(U))
  {
    U1<-cbind(1,U)
    P1<-U1%*%solve(t(U1)%*%U1)%*%t(U1)
    ssreg1<-y%*%(diag(n)-P1)%*%y
    df1.3<-1
    df2.3<-n-k-2
    F.value3<-(ssreg1-ssreg2)/df1.3
    F.value3<-F.value3/(ssreg2/df2.3)
    p.value3<-1-pf(F.value2,df1.3,df2.3)
  }
  else
  {
    df1.3<-NA
    df2.3<-NA
    F.value3<-NA
    p.value3<-NA
  }

  res<-list(y=y,x=x,alpha=alpha,c=c,gcv=res.gcv$gcv,ghat=ghat,
            F.value1=F.value1,df1.1=df1.1,df2.1=df2.1,p.value1=p.value1,
            F.value2=F.value2,df1.2=df1.2,df2.2=df2.2,p.value2=p.value2,
            F.value3=F.value3,df1.3=df1.3,df2.3=df2.3,p.value3=p.value3)
  class(res) <- "css"
  return(res)
}
