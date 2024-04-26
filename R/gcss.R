#' Testing of spline growth curve model
#'
#' \code{gcss} is used to test the mean growth curves using cubic smoothing splines. 
#'
#' @param Y a \eqn{q\times n} matrix of response vectors. The \eqn{i}th column contains the 
#' \eqn{q}-variate response vector of the \eqn{i}th individual. The \eqn{j}th row contains the responses 
#' of the individuals in the \eqn{j}th time point.
#' @param A an \eqn{n\times m} between-individual design matrix. You should give either 
#' \code{A} or \code{grp}, not both.
#' @param grp an \eqn{n}-dimensional vector giving the \eqn{m} treatment groups of the observations. 
#' The \eqn{n\times m} between-individual design matrix \code{A} is constructed using the values of this vector:
#' The \eqn{i}th row of \code{A} is a unit vector with 1 in position \code{grp[i]} and zeros elsewhere.
#' @param t a \eqn{q}-dimensional vector of the time points.
#' @param alpha fixed value of \eqn{alpha}. If it is not given, \eqn{alpha} is estimated using gcv criteria. 
#' The minimum of gcv criteria is found using grid search.
#' @param alpha.min the lower bound of grid of \eqn{alpha} values when estimating \eqn{alpha} using gcv criteria.
#' @param alpha.max the upper bound of grid of \eqn{alpha} values when estimating \eqn{alpha} using gcv criteria.
#' @param c fixed value of the number of eigenvectors.
#' @param model covariance structure. The default "unif" assumes the uniform covariance structure
#' \eqn{R = d^2 1_q1_q'+I_q}. 
#' @details 
#' Here are the details of the function...
#' @return A list containing the following components:
#' \describe{
#' \item{F.value}{the observed value of the F test statistic.}
#' \item{df1}{the degrees of freedom of the numerator of the F test statistic.}
#' \item{df2}{the degrees of freedom of the denominator of the F test statistic.}
#' \item{p.value}{approximate p-value.}
#' \item{p.perm}{permutation test p-value.}
#' \item{alpha}{the value of \eqn{\alpha}.}
#' \item{c}{the value of c.}
#' \item{Ghat}{the spline fit.}
#' \item{Gtilde}{the approximated spline fit.}
#' }
#' @importFrom stats pf
#' @export
gcss<-function(Y, A=NULL, grp=NULL, t=1:dim(Y)[1], alpha=NULL, alpha.min=0.1, alpha.max=1000, c=NULL, model="unif")
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
  print(A)
  q<-nrow(Y)
  n<-ncol(Y)
  m<-ncol(A)
  K<-roughness(t)
  
  if(is.null(alpha))
  {
    #Estimation of alpha using gcv criteria
    gcv.res<-gcv.alpha(Y,t,A,alpha.min=alpha.min,alpha.max=alpha.max,len=100)
    if(gcv.res$minpoint==(-1))
      cat("The minimum gcv was found at the point alpha.min=",alpha.min,"\n")
    if(gcv.res$minpoint==1)
      cat("The minimum gcv was found at the point alpha.max=",alpha.max,"\n")
    alpha<-gcv.res$alpha.hat
  }

  S<-solve(diag(q)+alpha*K)
  
  Ghat<-S%*%Y%*%A%*%solve(t(A)%*%A)
  
  M<-eigen(S)$vectors
  m1<-rep(1,q)/sqrt(q)
  St<-sqrt(sum((t-mean(t))^2))
  m2<-(t-mean(t)*rep(1,q))/St
  M[,1]<-m1
  M[,2]<-m2
  
  if(is.null(c))
  {
    #Esimation of the number of eigenvectors c using gcv criteria
    c<-gcv2.dim(Y,A,M)$c
  }
  
  Mstar<-M[,1:c]
  Pm<-Mstar%*%t(Mstar)
  Gtilde<-Pm%*%Y%*%A%*%solve(t(A)%*%A)
  
  Omegahat<-t(Mstar)%*%Y%*%A%*%solve(t(A)%*%A)
  
  C<-cbind(rep(0,c-1),diag(c-1)) #drop the first eigenvector corresponding to the constant term
  D<-rbind(rep(1,m-1),-diag(m-1))  #test if the progression is the same in the m groups
  nu<-nrow(C)
  g<-ncol(D)
  COD<-C%*%Omegahat%*%D
  Qstar<-sum(diag(COD%*%solve(t(D)%*%solve(t(A)%*%A)%*%D)%*%t(COD)))
  sigmahat2<-sum(diag(t(Y)%*%(diag(q)-Pm)%*%Y))/(n*(q-c)); sigmahat2
  FF<-(Qstar/(nu*g))/sigmahat2
  FF<-c(FF)
  df1<-nu*g
  df2<-n*(q-c)
  p.value<-1-pf(FF,df1,df2)
  
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
  res<-list(F.value=FF,df1=df1,df2=df2,p.value=p.value,p.perm=p.perm,
            alpha=alpha,c=c,Ghat=Ghat,Gtilde=Gtilde)
  class(res) <- "gcss"
  return(res)
}
