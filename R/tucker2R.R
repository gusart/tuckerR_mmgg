#' Three-Mode Principal Components: Tucker 2 Model
#'
#' @description This function performs Three-Mode Principal Components  using
#' Tucker-2 Model.Compute all the output necessary to plot  interactive
#' Biplot.The Three-Mode Principal Component Analysis, provides both useful
#' analytic and graphic tools to study and characterize phytogenetic resources,
#' especially when the influence of environmental factors are possible.
#'
#' @usage tucker2R(datos, amb= 2, stand = TRUE, nc1 = 2, nc2 = 2, niter = 10000)
#'
#' @param datos a data frame with n rows for individuals and p variable for
#' columns. All the conditions must be the same variables names and
#' individuals.
#'
#' @param  amb The diferent conditions, in which the same variables and
#' individuals had been studied. By default is 2.
#' @param  stand a boolean value, if it is TRUE (value set by default) each
#' variable is centered and scale by variable.
#' @param nc1 number of components in the first mode, by default is 2
#' @param nc2 number of components in the second mode, by default is 2
#' @param niter the iteration number for the Tuckals algorithm, by default
#' 10000 iteration
#'
#' @details To determine the number of components that are going to be retained,
#'  we use previously to the algorithm applications,method called DifFit. The
#'  number of components in the third mode is obtained from the number of
#'  conditions.The labels of the variables must be the same for all conditions
#'  in the data frame.
#'
#' @return \code{Resultado} a list which stores the name of the individual and
#' the variables, the number of iterations, the variability explained by the
#' model, and the total variability.
#' \code{Proyeccion} It is a list which holds the projection of individuals and
#' variables to see if the biplot is difficult to understand because of  the
#' huge number of cases or plotted vectors.
#' \code{saltuck} is a list with the results of the algorithm to plot the biplot,
#'  where the names of the conditions are well kept.
#'
#' @references
#' \describe{
#'  \item{Marticorena, M.; Bramardi, S.; Defacio, R. 2010.}{Characterization of maize populations in different environmental
#' conditions by means of Three-Mode Principal Components Analysis. Revista Ciencia e Investigacion   Agraria. 37(3): 93-105.}
#'  \item{Timmerman, M.E., and H. Kiers. 2000.}{Three-mode principal components analysis. Choosing numbers of components and sensitivity to local
#' optima. The British Journal of the Mahematical and Statistical Psychology 53: 1-16.}
#' }
#'
#' @author Marta Marticorena, Gustavo Gimenez, Cecilia Gonzalez, Sergio Bramardi
#'
#' @seealso The function plot.marta for a complete analisis.
#'
#' @examples
#' data(maize_pop,package = "tuckerR.mmgg")
#' (output <- tucker2R(maize_pop,amb=2,stand=TRUE,nc1=3,nc2=3))
#'
#' @keywords kwd1
#'
#' @export

tucker2R <- function(datos,amb=2,stand=TRUE,nc1=2,nc2=2,niter=10000){
  niter <- niter
  nc3 <- amb

  if(!is.data.frame(datos)){stop("datos must be a data frame")}

  if(ncol(datos)%%amb!=0){stop("The variables must be the same in amb")}

  if (any(is.na(datos))){stop("There is at least one NA  'datos' must be complete")}

  if(amb>4){print("when amb > 4 The Tucker 3 is recommended!")}

  colnume <- ncol(datos)
  colnamb <- colnume/amb
  ambi <- amb-1
  for (j in 1:ambi){
    t <- -1+j
    t1 <- t*colnamb
    for (i in 1:colnamb){
      a <- i+t1
      b <- colnamb +i + t1
      fff <- colnames(datos)[a] == colnames(datos)[b]
      if (fff==FALSE){stop("datos must be the same variables")}
    }
  }

  if (stand == TRUE){datos.stan <- scale(datos)
  center <- attributes(datos.stan)$`scaled:center`
  Escala <- attributes(datos.stan)$`scaled:scale`
  datos <- as.data.frame(datos.stan)
  }else{
    datos <- datos
  }

  if(!(nc1<=nc2*nc3 && nc3<=nc1*nc2 && nc2 <= nc1*nc3)){stop("in the combination of components for a solution must use a Diffit criteria")}

  I <-dim(datos)[1]
  J <-dim(datos)[2]
  J <- J/amb
  K<-amb
  col <- J

  n1<-nc1
  n2<-nc2
  n3<-nc3


  etiq_var <- colnames(datos)[1:J]
  etiq_varXamb <- rep(etiq_var,K)
  eti <- rep("amb",K)
  neti <- seq(1:K)
  etiq_amb <- paste(eti,neti,sep=".")
  etiq_ambXJ <- rep(etiq_amb,each=J)
  etiq <- paste(etiq_var,etiq_ambXJ,sep="-")
  etiq_ind <- rownames(datos)


  matricex <- matrition(datos,I,J,K)
  m <- matricex$m
  X1 <- matricex$X1
  X2 <- matricex$X2
  X3 <- matricex$X3

  rownames(m) <- etiq_ind
  colnames(m) <- etiq_var
  dimnames(m)[3] <-list(etiq_amb)


  W1 <- 0
  x <- X1%*%t(X1)
  y <- X2%*%t(X2)
  z <- X3%*%t(X3)

  p <- 0
  for (u in 1:I)
  {
    p<-p+x[u,u]
  }

  tam1  <- nrow(X1)
  tam2  <- ncol(X1)


  descompx <- svd(x)
  a  <- descompx$u
  d1 <- descompx$d
  v  <- descompx$v

  descompy <- svd(y)
  b  <- descompy$u
  d2 <- descompy$d
  v  <- descompy$v

  corrida_tuckal <- tuckal(a,b,d2,n1,n2,n3,X1,X2,p,niter)
  iter <- corrida_tuckal$iter
  if (iter >= niter) print("Warnings: not converg")

  a <- corrida_tuckal$a
  A <- corrida_tuckal$A
  B <- corrida_tuckal$B
  C <- corrida_tuckal$C


  SCE  <- corrida_tuckal$SCE

  G1 <- corrida_tuckal$G1
  K <- kronecker(C,B)
  D <- G1%*%t(K)
  D <- t(D)

  sca  <- 0
  scb  <- 0
  sca  <- sum(sum(A^2))
  scb  <- sum(sum(D^2))
  sca  <- sca/tam1
  scb  <- scb/tam2
  scf  <- sqrt(sqrt(scb/sca));
  IND  <- a*scf
  INDi <- IND[,1:2]
  colnames(INDi) <- c("Dim-1","Dim-2")
  D <- D/scf
  rownames(D)    <- etiq
  rownames(IND)  <- etiq_ind
  Dvar           <- D[,1:2]
  colnames(Dvar) <- c("Dim-1","Dim-2")

  if (stand == TRUE){
    Resultado <- list(Scales=Escala,MEANS=center,individuos=etiq_ind,variables=etiq,
                      iteraciones=iter,SCExplicada=SCE,vartot=p)
  } else {
    Resultado <- list(individuos=etiq_ind,variables=etiq,
                      iteraciones=iter,SCExplicada=SCE,vartot=p)
  }

  Proyeccion <- list(variables=Dvar,Individuos=INDi)
  saltuck <- list(IND=IND,D=D,Ambientes=etiq_amb,
                  Resultados=Resultado,Proyecciones=Proyeccion,matrizG=G1)#modificado
  class(saltuck) <- "marta"
  return(saltuck)
}
