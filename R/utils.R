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

#' Generalized cross-validation
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
#' \item{minpoint}{a logical value. If TRUE, the minimum of gcv was found at  
#' \eqn{alpha[ind.min]}, where \eqn{1<ind.min<length(alpha)}. If FALSE, the 
#' minimum of gcv was found at \eqn{alpha[1]} or \eqn{alpha[length(alpha)]}.}
#' }
#' @keywords internal
gcv.alpha<-function(Y,x,A,alpha.min,alpha.max,len=20)
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
  minpoint<-TRUE
  if((ind.min==1)|(ind.min==len))
    minpoint<-FALSE
  alpha.hat<-alpha[ind.min]
  plot(alpha,gcv,type="l")
  abline(v=alpha.hat,lty=2)
  res<-list(gcv=gcv,alpha=alpha,alpha.hat=alpha.hat,ind.min=ind.min,minpoint=minpoint)
  return(res)
}


#' Generalized cross-validation
#'
#' \code{gcv.dim} estimates the number of eigenvectors by using generalized cross-validation criteria
#'
#' @param Y a \eqn{q \times n} matrix of independent \eqn{q \times 1} response vectors.
#' @param A an \eqn{n \times m} between-individual design matrix.
#' @param M matrix of eigenvectors of \code{S}.
#' @details 
#' Here are the details of the function...
#' @return the number of eigenvectors.
#' @keywords internal
gcv.dim<-function(Y,A,M)
{
  n<-ncol(Y)
  q<-nrow(Y)
  m<-ncol(A)
  Pa<-A%*%solve(t(A)%*%A)%*%t(A)
  gcv<-NULL
  for(i in 1:q)
  {
    Ti<-M[,1:i]%*%t(M[,1:i])
    Yhat<-Ti%*%Y%*%Pa
    D<-Y-Yhat
    trDD<-sum(diag(D%*%t(D)))
    gcv[i]<-(1/(n*q))*trDD/(1-m*i/(n*q))^2
  }
  ind.min<-which.min(gcv)
  cc<-(1:q)[ind.min]
  res<-cc
  return(res)
}
