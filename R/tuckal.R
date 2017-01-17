#########################Tuckal v2
#This function fits tuckals algorithm
tuckal <- function(a,b,d2,n1,n2,n3,X1,X2,p,niter){
  W1 <- 0
  a <- a[,1:n1]
  b <- b[,1:n2]
  c <- diag(n3)

  j1     <- 1
  iter   <-  0
  while (abs(j1)>=0.05 && iter < niter ){
    k1 <- kronecker(c,b)
    G1 <- crossprod(a,X1)%*%k1
    r1 <- kronecker(t(c),t(b))
    S1 <- a%*%G1%*%r1
    t1 <- tcrossprod(X1-S1)
    l1 <- sum(diag(t1))
    m  <- l1
    j1 <- l1-W1
    x  <- X1%*%k1
    x  <- x%*%t(x)
    descompx1 <- svd(x)
    a1 <- descompx1$u
    d2 <- descompx1$d
    v2 <- descompx1$v
    a  <-  a1[ ,1:n1]
    G1 <- crossprod(a,X1)%*%k1
    S1 <- a%*%G1%*%r1
    t1 <- tcrossprod(X1-S1)
    l2 <- sum(diag(t1))
    m  <-  l2
    j1 <- l2-l1
    k2 <- kronecker(c,a)
    y  <- X2%*%k2
    y  <- y%*%t(y)
    descompy1 <- svd(y)
    b1 <- descompy1$u
    d3 <- descompy1$d
    v3 <- descompy1$v
    b  <- b1[ ,1:n2]
    k1 <- kronecker(c,b)
    G1 <- crossprod(a,X1)%*%k1
    r1 <- kronecker(t(c),t(b))
    S1 <- a%*%G1%*%r1
    t1 <- tcrossprod(X1-S1)
    l3 <- sum(diag(t1))
    m  <- l3
    j1 <- l3-l2
    iter <- iter + 1
    W1 <- l3
  }

  A <- a
  B <- b
  C <- c

  SCE <- (p-m)/p * 100

  sal_tuckal <- list(a=a,b=b,c=c,A=A,B=B,C=C,m=m,G1=G1,SCE=SCE,iter=iter)
  return(sal_tuckal)
}
