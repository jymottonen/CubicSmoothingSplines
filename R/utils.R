#' Roughness matrix
#'
#' \code{roughness} is used to compute the roughness matrix. 
#'
#' @param x a \code{q}-dimensional vector of the time points.
#' @details 
#' Here are the details of the function...
#' @return A \eqn{q\times q} roughness matrix \eqn{K=\nabla\Delta^{-1}\nabla'}.
#' @keywords internal
roughness<-function(x)
{
  q<-length(x)
  h<-x[-1]-x[-q]
  
  Nabla<-matrix(0,q,q-2)
  for(i in 1:q){
    for(j in 1:(q-2)){
      if(i==j)
        Nabla[i,j]<-1/h[i]
      if(i==(j+1))
        Nabla[i,j]<- -(1/h[j]+1/h[i])
      if(i==(j+2))
        Nabla[i,j]<- 1/h[j+1]
    }
  }
  
  Delta<-matrix(0,q-2,q-2)
  for(i in 1:(q-2)){
    for(j in 1:(q-2)){
      if(i==j)
        Delta[i,j]<-(h[i]+h[i+1])/3
      if((i+1)==j)
        Delta[i,j]<-h[j]/6
      if(i==(j+1))
        Delta[i,j]<-h[i]/6
    }
  }

  roughness<-Nabla%*%solve(Delta)%*%t(Nabla)
  return(roughness)
}

#' Generalized cross-validation estimate of the smoothing parameter (growth data model)
#'
#' \code{gcv.alpha} estimates the smoothing parameter \eqn{\alpha} by using generalized cross-validation criteria
#'
#' @param Y a \eqn{q \times n} matrix of independent \eqn{q \times 1} response vectors.
#' @param x a \code{q}-dimensional vector of the time points.
#' @param A an \eqn{n \times m} between-individual design matrix.
#' @param alpha.min the minimum value of the grid of \eqn{\alpha}s.
#' @param alpha.max the maximum value of the grid of \eqn{\alpha}s.
#' @param len the number of values in the grid of \eqn{\alpha}s. 
#' @details 
#' Here are the details of the function...
#' @return A list containing the following components:
#' \describe{
#' \item{gcv}{a vector of gcv values.}
#' \item{alpha}{a vector of alpha values.}
#' \item{alpha.hat}{the alpha value that corresponds to the minimum of gcv.}
#' \item{ind.min}{the index of alpha that corresponds to the minimum of gcv.}
#' \item{ind.min}{the index of alpha that corresponds to the minimum of gcv.}
#' \item{minpoint}{Returns \eqn{0} if \eqn{1<ind.min<length(alpha)}, \eqn{-1} if \eqn{ind.min=1} and \eqn{+1} if \eqn{ind.min=length(alpha)}.}
#' }
#' @keywords internal
gcv.alpha<-function(Y,x,A,alpha.min,alpha.max,len=100)
{
  n<-ncol(Y)
  m<-ncol(A)
  K<-roughness(x)
  Pa<-A%*%solve(t(A)%*%A)%*%t(A)
  q<-length(x)
  alpha<-seq(alpha.min,alpha.max,length=len)
  gcv<-NULL
  for(i in 1:len)
  {
    S<-solve(diag(q)+alpha[i]*K)
    Yhat<-S%*%Y%*%Pa
    D<-Y-Yhat
    trDD<-sum(diag(D%*%t(D)))
    edf<-sum(diag(S))
    gcv[i]<-(1/(n*q))*trDD/(1-m*edf/(n*q))^2
  }
  ind.min<-which.min(gcv)
  minpoint<-0
  if(ind.min==1)
    minpoint <- -1
  if(ind.min==len)
    minpoint <- 1
  alpha.hat<-alpha[ind.min]
  res<-list(gcv=gcv,alpha=alpha,alpha.hat=alpha.hat,ind.min=ind.min,minpoint=minpoint)
  return(res)
}

#' Generalized cross-validation (cubic smoothing splines)
#'
#' \code{gcv1.dim} estimates the number of eigenvectors by using generalized cross-validation criteria (cubic smoothing splines)
#'
#' @param Y a \eqn{q \times n} matrix of independent \eqn{q \times 1} response vectors.
#' @param A an \eqn{n \times m} between-individual design matrix.
#' @param M matrix of eigenvectors of \code{S}.
#' @details 
#' Here are the details of the function...
#' @return A list containing the following components:
#' \describe{
#' \item{c}{the number of eigenvectors.}
#' \item{gcv}{the values of the gcv criteria.}
#' }
#' @keywords internal
gcv1.dim<-function(y,Tn)
{
  n<-length(y)
  gcv<-NULL
  for(i in 1:(n-1))
  {
    Mc<-Tn[,1:i]%*%t(Tn[,1:i])
    ss<-mean((y-Mc%*%y)^2)
    gcv[i]<-ss/(1-(i/n))^2
  }
  ind.min<-which.min(gcv)
  c<-(1:n)[ind.min]
  res<-list(c=c,gcv=gcv)
  return(res)
}

#' Generalized cross-validation estimate of the number of eigenvectors (growth data model)
#'
#' \code{gcv2.dim} estimates the number of eigenvectors by using generalized cross-validation criteria (growth data model)
#'
#' @param Y a \eqn{q \times n} matrix of independent \eqn{q \times 1} response vectors.
#' @param A an \eqn{n \times m} between-individual design matrix.
#' @param M matrix of eigenvectors of \code{S}.
#' @details 
#' Here are the details of the function...
#' @return A list containing the following components:
#' \describe{
#'   \item{c}{the estimated number of eigenvectors.}
#'   \item{c.grid}{the grid of the number of eigenvectors for the gcv criteria.}
#'   \item{gcv}{the gcv values at the points \eqn{c.grid}}
#' }
#' @keywords internal
gcv2.dim<-function(Y,A,M)
{
  n<-ncol(Y)
  q<-nrow(Y)
  m<-ncol(A)
  Pa<-A%*%solve(t(A)%*%A)%*%t(A)
  gcv<-NULL
  for(i in 3:(q-1))
  {
    Ti<-M[,1:i]%*%t(M[,1:i])
    Yhat<-Ti%*%Y%*%Pa
    D<-Y-Yhat
    trDD<-sum(diag(D%*%t(D)))
    gcv[i]<-(1/(n*q))*trDD/(1-m*i/(n*q))^2
  }
  ind.min<-which.min(gcv)
  c<-(1:(q-1))[ind.min]
  res<-list(c=c,c.grid=(3:(q-1)),gcv=gcv)
  return(res)
}

