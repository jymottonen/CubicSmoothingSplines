#' plot.css
#'
#' Plot method for objects of class "css". 
#'
#' @param x an object of class css.
#' @param xlab a title for the x axis.
#' @param ylab a title for the y axis.
#' @param main a main title.
#' @param xlim numeric vector of length 2, giving the x coordinate range.
#' @param ylim numeric vector of length 2, giving the y coordinate range.
#' @param both logical. If TRUE, both smoothing curve and observed points are plotted. 
#' If FALSE, only smoothing curve is plotted.
#' @param ... further arguments passed to or from other methods.
#' @details 
#' Here are the details of the function...
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline axis lines points title
#' @export
plot.css<-function(x, xlab="x", ylab="y", xlim=c(min(x$x),max(x$x)), ylim=c(min(x$y),max(x$y)), 
                   main="Approximate spline fit", both=FALSE, ...)
{
  if(both)
  {
    plot(x$y~x$x, ylab=ylab, xlab=xlab, xlim=xlim, ylim=ylim, main=main)
    lines(x$ghat~x$x, type="l")
  }
  else
    plot(x$ghat~x$x, type="l", ylab=ylab, xlab=xlab, xlim=xlim, ylim=ylim, main=main)
}

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


#' plot.mgcss
#'
#' Plot method for objects of class "mgcss". 
#'
#' @param x an object of class mgcss.
#' @param t a \eqn{q}-dimensional vector of the time points.
#' @param xlab a title for the x axis.
#' @param ylab a title for the y axis.
#' @param ... further arguments passed to or from other methods.
#' @details 
#' Here are the details of the function...
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline axis lines points title par
#' @export
plot.mgcss<-function(x, t=1:dim(x$Gtilde[[1]])[1], xlab="Time points", ylab="Gtilde", ...)
{
  s<-length(x$Gtilde)
  q<-nrow(x$Gtilde[[1]])
  m<-ncol(x$Gtilde[[1]])
  
  #op<-par(mar=c(1,1,1,1))
  par(ask=TRUE)
  for(j in 1:s)
  {
    plot(t,x$Gtilde[[j]][,1],type="l", xlab=xlab, ylab=ylab)
    for(i in 2:m)
    {
      points(t,x$Gtilde[[j]][,i],type="l",lty=i)
    }
  }
  par(ask=FALSE)
  #par(op)
}


