#' summary.css
#'
#' summary.css is used to print the summary of the
#' testing of cubic smoothing splines and semi-parametric regression models
#'
#' @param object an object of class gcss.
#' @param ... further arguments passed to or from other methods.
#' @param digits minimal number of significant digits.
#' @details 
#' Here are the details of the function...
#' @export
summary.css<-function(object, ..., digits=3)
{
  cat("Testing of cubic smoothing splines and semi-parametric regression models\n\n")
  cat("The smoothing parameter alpha:",object$alpha,"\n")
  cat("The number of eigenvectors c:",object$c,"\n\n")
  
  test.table<-cbind(object$F.value1,object$df1.1,object$df2.1,object$p.value1)
  colnames(test.table)<-c("F-value","df1","df2","p-value")
  rownames(test.table)<-c("")
  cat("Test 1\n")
  print(test.table)
  
  if(!is.na(object$F.value2))
  {
    test.table2<-cbind(object$F.value2,object$df1.2,object$df2.2,object$p.value2)
    colnames(test.table2)<-c("F-value","df1","df2","p-value")
    rownames(test.table2)<-c("")
    cat("\nTest 2\n")
    print(test.table2)
  }
  
  if(!is.na(object$F.value3))
  {
    test.table3<-cbind(object$F.value3,object$df1.3,object$df2.3,object$p.value3)
    colnames(test.table3)<-c("F-value","df1","df2","p-value")
    rownames(test.table3)<-c("")
    cat("\nTest 3\n")
    print(test.table3)
  }
}

#' summary.gcss
#'
#' summary.gcss is used to print the summary of the
#' testing of spline growth curve model.
#'
#' @param object an object of class gcss.
#' @param ... further arguments passed to or from other methods.
#' @param digits minimal number of significant digits.
#' @details 
#' Here are the details of the function...
#' @export
summary.gcss<-function(object, ..., digits=3)
{
  test.table<-cbind(object$F.value,object$df1,object$df2,object$p.value,object$p.perm)
  colnames(test.table)<-c("F-value","df1","df2","p-value","p-perm")
  rownames(test.table)<-c("")
  cat("Testing of spline growth curve model\n\n")
  cat("The smoothing parameter alpha:",object$alpha,"\n")
  cat("The number of eigenvectors c:",object$c,"\n\n")
  print(test.table)
}
