#' Aplication of the Diffit method.
#'
#' @aliases print.diff
#'
#' @description The diffit method is used to apply when we need to know the
#' axis number to be gathered in the P mode, and Q mode. The third mode, K
#' it is related to the environment  numbers.
#' The diffit method consist on fitting each value with the Tuckle algorithm.
#'
#' @param datos datos original data from data frames
#' @param amb numbers of environment
#' @param stand a boolean value, if it is TRUE (value set by default) each
#' variable is centered and scale by variable.
#' @param niter the iteration number for the Tuckals algorithm, by default
#' 10000 iteration.
#'
#' @details The final result is the model which has the most coefficient diffits
#' the greatest variability explained and the one which exceed the threshold.
#'
#' @return \code{saldiffit} a list with a combination numbers of axis,
#' percentage of variability explained and Diffit rate. The critic value or
#' threhold is also return.
#'
#' @references
#' \describe{
#'  \item{MARTICORENA, M.; BRAMARDI, S.; DEFACIO, R. 2010.}{Characterization of
#'  maize populations in different environmental conditions by means of
#'  Three-Mode Principal Components Analysis. Revista Ciencia e Investigacion
#'  Agraria. 37(3): 93-105.}
#'  \item{Timmerman, M.E., and H. Kiers. 2000.}{Three-mode principal components
#'   analysis. Choosing numbers of components and sensitivity to local optima.
#'   The British Journal of the Mahematical and Statistical Psychology 53: 1-16.}
#' }
#'
#' @author Marta Marticorena, Gustavo Gimenez, Cecilia Gonzalez, Sergio Bramardi
#'
#' @examples
#' #Copy and paste this example in your console without the comment
#' #data(maize_pop,package = "tuckerR.mmgg")
#' #dif_sal <- diffit(maize_pop,amb=2)
#' #print(dif_sal) the best combination is 3 3 2
#'
#' @export

diffit<- function(datos,amb=2,stand = TRUE,niter=10000){
  if (stand == TRUE){datos.stan <- scale(datos)
  center <- attributes(datos.stan)$`scaled:center`
  Escala <- attributes(datos.stan)$`scaled:scale`
  datos <- as.data.frame(datos.stan)
  }else{
    datos <- datos
  }
  K <- amb
  Q <- ncol(datos)/amb
  P <- nrow(datos)

  I <-dim(datos)[1]
  J <-dim(datos)[2]
  J <- J/amb
  K<-amb
  col <- J

  matricex <- matrition(datos,I,J,K)
  m <- matricex$m
  X1 <- matricex$X1
  X2 <- matricex$X2
  X3 <- matricex$X3

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

  all.c <- expand.grid(P=1:P,Q=1:Q,K=K,stringsAsFactors = FALSE)

  t <- 0
  retenidos <- data.frame(P=0,Q=0,K=0)
  for ( i in 1:nrow(all.c)) {
    if(all.c[i,1]<=all.c[i,2]*all.c[i,3] && all.c[i,3]<=all.c[i,1]*all.c[i,2] && all.c[i,2]<=all.c[i,1]*all.c[i,3]){
      t <- t+1
      retenidos[t,] <- as.vector(all.c[i,])
    }
  }

  retenidos$S <-  rowSums(retenidos)
  retenidos.ord <- retenidos[order(retenidos$S,retenidos$P),]

  li <- nrow(retenidos.ord)

  SCE_diffit <- NULL
  iter_diffit <- NULL
  for (ii in 1:li) {
    n1<-retenidos.ord[ii,"P"]
    n2<-retenidos.ord[ii,"Q"]
    n3<-retenidos.ord[ii,"K"]

    corrida_tuckal <- tuckal(a,b,d2,n1,n2,n3,X1,X2,p,niter)

    SCE_diffit <- c(SCE_diffit ,as.vector(round(corrida_tuckal$SCE,2)))
    iter_diffit <- c(iter_diffit ,as.vector(corrida_tuckal$iter))
  }
  datos_salida <- cbind(retenidos.ord[1:li,],SCE_diffit,iter_diffit)


  datos_salida2 <- datos_salida[datos_salida$iter_diffit != niter,]
  saltos <- unique(datos_salida2$S)
  modelos <- data.frame(P=0,Q=0,K=0,S=0,SCE_diffit=0,iter_diffit=0)
  it <- 0
  for (ij in saltos){
    dattt <- datos_salida2[datos_salida2$S == ij,]
    pos <- which.max(dattt$SCE_diffit)
    it <- it + 1
    modelos[it,] <- as.vector(dattt[pos,])
  }

  diffit_t <- c(modelos[1,"SCE_diffit"],diff(modelos[,"SCE_diffit"]))
  modelos$diffit <- diffit_t


  diffit_coc <- rep(NA,length(diffit_t))
  for (cd in 1:(length(diffit_t)-1)){
    if (diffit_t[cd] >  diffit_t[cd+1]) {diffit_coc[cd] <- diffit_t[cd]/diffit_t[cd+1]}
  }
  modelos$coc_diffit <- diffit_coc

  critic_value <- 100/(max(retenidos$S)-min(retenidos$S))

  saldiffit <- list(models=modelos,critic_value = critic_value)

  class(saldiffit) <- "diff"

  return(saldiffit)

}
