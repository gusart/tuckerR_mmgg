#' Concatenate data frame in array and matrix by cases, variables and
#' environments
#'
#' @param datos original data from data frames
#' @param I the numbers of cases
#' @param J the numbers of variables
#' @param K the numbers of environment or conditions
#'
#' @description Concatenate data frame in array and matrix by cases, variables
#'  and environments to performs three mode principal components with the
#'  function \code{\link{tucker2R}}.
#'
#' @return \code{matrizlista} return a list with: the array "m" with all the
#' data concatenate in array. X1 the data is concatenate by cases, X2  the data
#' concatenate by variables and X3 the data concatenate by environments.
#'
#' @details This process is also knowing as 'matricizing'  or 'unfolding'.
#'
#' @examples
#' data(maize_pop,package = "tuckerR.mmgg")
#' conc_matrix <- matrition(maize_pop,I=30,J=10,K=2)
#' conc_matrix$m  #get m array
#' conc_matrix$X1 #get matrix by cases
#' conc_matrix$X2 #get matrix by variables
#' conc_matrix$X3 #get matrix by environments
#' @export
matrition  <- function(datos,I,J,K){

  m <- array(as.matrix(datos),dim = c(I,J,K))


  X1 <- as.matrix(m)
  dim(X1) <- c(I,J*K)

  X2 <- array(0, dim=c(J,I*K))
  for (jj in 1:J)
  {
    c<-1;
    for (ii in 1:I)
    {
      for (kk in 1:K)
      {
        X2[jj,c]<-X2[jj,c]+m[ii,jj,kk];
        c<-c+1
      }
    }
  }


  x3m <- as.matrix(m)
  X3 <- matrix(x3m,K,I*J,byrow = TRUE)
  matrizlista <- list(X1=X1,X2=X2,X3=X3,m=m)
  return(matrizlista)
}
