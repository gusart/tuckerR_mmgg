tucker2R <- function(datos,amb=2,stand=TRUE,nc1=2,nc2=2,niter=1000){

   {nc3 <- amb}
     ##La condición para que el marco de datos sea únicamente un marco de datos
    if(!is.data.frame(datos)){stop("datos must be a data frame")}
    ##Condición para que tenga el mismo número de variables en los ambientes
    if(ncol(datos)%%amb!=0){stop("The variables must be the same in amb")}
    ##Condición para agregar un dato, en caso que sea un dato faltante.
    if (any(is.na(datos))){stop("There is at least one NA  'datos' must be complete")}
    ##Recomendación para utilizar el tucker 3 en lugar del tucker 2
    if(amb>4){print("when amb > 4 The Tucker 3 is recommended!")}
    ##Condición para tener las mismas variables en todos los ambientes
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

    ##Argumento de estandarización
    if (stand == TRUE){datos.stan <- scale(datos)
                       center <- attributes(datos.stan)$`scaled:center`
                       Escala <- attributes(datos.stan)$`scaled:scale`
                       datos <- as.data.frame(datos.stan)
                   }else{
                       datos <- datos
                   }
      ##Condición para que se cumpla el criterio de Diffit
    if(!(nc1<=nc2*nc3 && nc3<=nc1*nc2 && nc2 <= nc1*nc3)){stop("in the combination of components for a solution must use a Diffit criteria")}
#######Define las categorías de cada modo
    I<-dim(datos)[1]   #número de categorías del primer modo
    J <-dim(datos)[2]   #número de categorías del segundo modo
    J <- J/amb 
    K<-amb   #número de categorías del tercer modo,general
    col <- J
#######Define el número de componentes de cada modo por default
    n1<-nc1   #número de componentes del primer modo
    n2<-nc2 #número de componentes del segundo modo
    n3<-nc3   #número de componentes del tercer modo

#######Etiquetas de los ambientes    
    etiq_var <- colnames(datos)[1:J]
    etiq_varXamb <- rep(etiq_var,K)
    eti <- rep("amb",K)
    neti <- seq(1:K)
    etiq_amb <- paste(eti,neti,sep=".")
    etiq_ambXJ <- rep(etiq_amb,each=J)
    etiq <- paste(etiq_var,etiq_ambXJ,sep="-")
    etiq_ind <- rownames(datos)

    m <- array(0, dim=c(I,J,K))
    for (k in 1:K){
        n <- k-1
        l <- n*col
        for (i in 1:I)  
            {
                for(j in 1:J)
                    {
                        m[i,j,k]<-m[i,j,k]+datos[i,j+l];
                    }
            }
    }

    rownames(m) <- etiq_ind
    colnames(m) <- etiq_var
    dimnames(m)[3] <-list(etiq_amb)

                                        #Concatenación de matrices: Obtenemos X1, X2 y X3
    ##Concatena por individuo X1
    X1 <- array(0, dim=c(I,J*K))
    for (ii in 1:I)
        {
            c<-1; 
            for (kk in 1:K)
                {
                    for (jj in 1:J)
                        {
                            X1[ii,c]<-X1[ii,c]+m[ii,jj,kk];
                            c<-c+1
                        }
                }
        }

                                        #Concatena por variable
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

                                        #Concatena por Ambiente
    X3 <- array(0, dim=c(K,I*J))
    for (kk in 1:K)
        {
            c<-1; 
            for (jj in 1:J)
                {
                    for (ii in 1:I)
                        {
                            X3[kk,c]<-X3[kk,c]+m[ii,jj,kk];
                            c<-c+1
                        }
                }
        }


                                        #cálculo de las matrices  a y b 
    W1=0 
    x=X1%*%t(X1)

                                        #inercia total
                                        #calculo la traza de x
    p<-0
    for (u in 1:I)
        {
            p<-p+x[u,u]
        }
    y=X2%*%t(X2);
    z=X3%*%t(X3);
    tam1 =nrow(X1);
    tam2 =ncol(X1);
    j1=1;

                                        #descomposición  en valores singulares de la matriz x
    a=svd(x)$u
    d1=svd(x)$d
    v=svd(x)$v

                                        #descomposición en valores singulares de la matriz y
    b=svd(y)$u
    d2=svd(y)$d
    v=svd(y)$v

    a=a[,1:n1]
    b=b[,1:n2] 
    c=diag(n3) 

    iter <-  0 #iter contador de iteraciones
    while (abs(j1)>=0.05){
        k1=kronecker(c,b) 
        G1=t(a)%*%X1%*%k1 
        r1=kronecker(t(c),t(b))
        S1=a%*%G1%*%r1
        t1=(X1-S1)%*%t(X1-S1)
        l1=sum(diag(t1))
        m=l1
        j1=l1-W1
        k1=kronecker(c,b)
        x=X1%*%k1
        x=x%*%t(x) 
        a1=svd(x)$u
        d2=svd(x)$d
        v2= svd(x)$v
        a=a1[ ,1:n1]
        k1=kronecker(c,b)
        G1=t(a)%*%X1%*%k1
        r1=kronecker(t(c),t(b))
        S1=a%*%G1%*%r1
        t1=(X1-S1)%*%t(X1-S1)
        l2=sum(diag(t1))
        m=l2
        j1=l2-l1
        k2=kronecker(c,a)
        y=X2%*%k2
        y=y%*%t(y)
        b1=svd(y)$u
        d3=svd(y)$d
        v3= svd(y)$v
        b=b1[ ,1:n2]
        k1=kronecker(c,b)
        G1=t(a)%*%X1%*%k1
        r1=kronecker(t(c),t(b))
        S1=a%*%G1%*%r1
        t1=(X1-S1)%*%t(X1-S1)
        l3=sum(diag(t1))
        m=l3
        j1=l3-l2
        iter=iter + 1
        if (iter >= niter) print("Warnings: not converg")
        if (iter >= niter) break
        W1=l3}
    A=a
    B=b
    C=c

                                        #BIPLOT
                                        # Obtención de los marcadores para el Biplot Interactivo: Djk
    K=kronecker(C,B);
    D=G1%*%t(K);
    D=t(D);
                                        # disp(D)

                                        #  suma de cuadrados explicada, en porcentaje.
    SCE=(p-m)/p * 100

                                        #Reescalamiento óptimo
    sca=0; 
    scb=0;
    sca=sum(sum(A^2));
    scb=sum(sum(D^2));
    sca=sca/tam1;
    scb=scb/tam2;
    scf=sqrt(sqrt(scb/sca));
    IND=a*scf
    IND
    INDi <- IND[,1:2]
    colnames(INDi) <- c("Dim-1","Dim-2")
    D=D/scf
    rownames(D) <- etiq
    rownames(IND) <- etiq_ind
    Dvar <- D[,1:2]
    colnames(Dvar) <- c("Dim-1","Dim-2")

    if (stand == TRUE){
        Resultado <- list(Scales=Escala,MEANS=center,individuos=etiq_ind,variables=etiq,
                          iteraciones=iter,SCExplicada=SCE,vartot=p)
    }
    if (stand == FALSE){
        Resultado <- list(individuos=etiq_ind,variables=etiq,
                          iteraciones=iter,SCExplicada=SCE,vartot=p)
    }

    Proyeccion <- list(variables=Dvar,Individuos=INDi)
    saltuck <- list(IND=IND,D=D,Ambientes=etiq_amb,
                    Resultados=Resultado,Proyecciones=Proyeccion)
    class(saltuck) <- "marta"
    return(saltuck)
}
