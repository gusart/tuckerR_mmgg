plotear_tucker2 <-
function(saltuck){
    if(class(saltuck)!="marta") stop("Object must be of class 'marta'")
    X1 <- saltuck$IND
    X <- as.data.frame(X1)
    Z1 <- saltuck$D
    Z <- as.data.frame(Z1)
    numcol <- dim(X)[2]
    numfil <- dim(X)[1]
    etiq <- rownames(saltuck$IND)
    etiq2 <- rownames(saltuck$D)
    ambsa <- saltuck$Ambientes
    ambs <- length(ambsa)
                                        #Grafica el marco pero no los puntos
    plot(X[,1],
         X[,2],
         ylim=c(min(X,Z),max(X,Z)), 
         xlim=c(min(X,Z),max(X,Z)), type="n",
         ylab="DIM 2", xlab="DIM 1"
         )

    title(main = list("Interactive Biplot Tucker-v2", cex = 1.5,
              col = "red", font = 3))

    abline(h = 0, lty = 2); abline(v = 0, lty = 2) #ejes punteados
###grafica los puntos correspondientes a los individuos

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
                   lwd=1.7, col=colores[tt])
            text(Z[k,1],Z[k,2],labels=etiq2[k],
                 cex=0.6,pos=2)
        }
    }
    legend("topright",legend=c(ambsa),fill=colores[1:ambs])

}
