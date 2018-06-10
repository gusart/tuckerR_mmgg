#' Plot an interactive Biplot
#'
#' @aliases plot.marta
#'
#' @description The interactive Biplot consists of combining two of the modes,
#' obtaining markers for the individuals  and vectors for the variables that
#'  were concatenated with the conditions.
#'  To plot the interactive Biplot this function need the output for the
#'  tucker2R function.
#'
#' @param saltuck is a list with the results of the algorithm to plot the
#' biplot, where the names of the conditions are well kept.
#' @param  ... Arguments to be passed to plot.
#'
#' @details The interactive Biplot allows for the visualization of the
#' inter structure of the differents data tables.
#'
#' @references
#' \describe{
#'  \item{MARTICORENA, M.; BRAMARDI, S.; DEFACIO, R. 2010.}{Characterization of
#'   maize populations in different environmental conditions by means of Three
#'   Mode Principal Components Analysis. Revista Ciencia e Investigacion
#'   Agraria. 37(3): 93-105.}
#'  \item{Gabriel, K.R. 1971.}{The biplot graphic display of matrices with
#'  appications to principal components analysis. Biometrika. 58, 453-467.}
#' }
#' @author Marta Marticorena, Gustavo Gimenez, Cecilia Gonzalez, Sergio Bramardi
#'
#' @examples
#' data(maize_pop,package = "tuckerR.mmgg")
#' prueba1 <- tucker2R(maize_pop, amb=2, stand=TRUE, nc1=3, nc2=3)
#' plot(prueba1)
#'
#' @importFrom graphics par plot.default rect text points arrows
#' @importFrom graphics abline legend title
#' @rdname plot
#' @export
plot <- function(saltuck, ...){
  UseMethod("plot")
}
#'@return \code{NULL}
#'
#'@rdname plot
#'@export
plot.marta <- function(saltuck){
  if(class(saltuck)!="marta"){
    print("Object must be of class marta")
  }else{
    X1 <- saltuck$IND
    X <- as.data.frame(X1)
    Z1 <- saltuck$D
    Z <- as.data.frame(Z1)
    numcol <- dim(X)[2]
    numfil <- dim(X)[1]
    etiq  <- saltuck$Resultados$individuos
    etiq2 <- saltuck$Resultados$variables
    ambsa <- saltuck$Ambientes
    ambs  <- length(ambsa)

    plot.default(X[,1],
         X[,2],
         ylim=c(min(X,Z),max(X,Z)),
         xlim=c(min(X,Z),max(X,Z)), type="n",
         ylab="DIM 2", xlab="DIM 1"
    )

    title(main = list("Interactive Biplot Tucker-v2", cex = 1.5,
                      col = "red", font = 3))

    abline(h = 0, lty = 2); abline(v = 0, lty = 2)


    for (i in 1:(numfil)) {
      points(X[i,1],X[i,2],pch=20, col="blue")
      text(X[i,1],X[i,2],labels=etiq[i],
           cex=0.6,pos=2)
    }

    dim(Z)[2] -> numcolz
    dim(Z)[1] -> numfilz

    ga <- 0:ambs
    fa <- ambs:1
    filXamb <- numfilz/ambs
    colores <- c("red","green","cyan","brown","coral","darkmagenta","hotpink")
    for (tt in 1:ambs){
      f <- 0
      gh <- 0
      f <- fa[tt]
      gh <- (filXamb)*ga[tt]
      gi <- ga[tt]+1
      for (k in (1+gh):(filXamb*gi)){
        arrows(0,0,Z[k,1],Z[k,2],length=0.05,
               lwd=1.5, col=colores[tt])
        text(Z[k,1],Z[k,2],labels=etiq2[k],
             cex=0.6,pos=2)
      }
    }
    legend("topright",legend=c(ambsa),fill=colores[1:ambs])

  }
}
