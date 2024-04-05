#' plot.gcss
#'
#' Plot method for objects of class "gcss". 
#'
#' @param x an object of class gcss.
#' @param t a \eqn{q}-dimensional vector of the time points.
#' @param xlab a title for the x axis.
#' @param ylab a title for the y axis.
#' @param ... further arguments passed to or from other methods.
#' @details 
#' Here are the details of the function...
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline axis lines points title
#' @export
plot.gcss<-function(x, t=1:dim(x$Gtilde)[1], xlab="Time points", ylab="Gtilde", ...)
{
  q<-nrow(x$Gtilde)
  m<-ncol(x$Gtilde)
  
  ymin<-min(x$Gtilde)
  ymax<-max(x$Gtilde)

  plot(t,x$Gtilde[,1],type="l", ylim=c(ymin,ymax), xlab=xlab, ylab=ylab)
  for(i in 2:m)
  {
    points(t,x$Gtilde[,i],type="l",lty=i)
  }
}

