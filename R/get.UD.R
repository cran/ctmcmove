get.UD <- function(R, method="lu", maxiter, start, tol){
  
  ## subset if needed
  zero.idx = which(colSums(R)==0)
  n.0=length(zero.idx)
  if(length(zero.idx)>0){
    R=R[-zero.idx,-zero.idx]
  }
  
  ## get infinitessimal generator
  
  n=nrow(R)
  one=Matrix(1,n,1)
  
  if(method=="lu"){
    G=-R
    diag(G)=as.numeric(R%*%one)
    
    ## get stationary distribution of CTMC
    ## G' %*% pi = 0
    
    ## See pg 455 of Harrod and Plemmons (1984) "Comparison of some direct methods for computing stationary distributions of Markov chains"
    
    LU.decomp=lu(t(G))
    LU=expand(LU.decomp)
    P=LU$P
    Q=LU$Q
    L=LU$L
    U=LU$U
    ## so t(G)=P'LUQ
    ## max(abs(t(G)-t(P)%*%L%*%U%*%Q)/max(G))
    
    ## unit vector
    e.n=Matrix(0,n,1,sparse=T)
    e.n[n,1]=1
    
    ## solve (UQ)pi=e.n
    U[n,n]=1
    Qpi=solve(U,e.n)
    pi=t(Q)%*%Qpi
  } else if(method=="limit"){
    if(missing(tol)) tol = 1.0e-7
    P=R
    P = Diagonal(x=1/as.numeric(P%*%one)) %*% P
    if(missing(start)){
      pi1 = one/sum(one)
    } else pi1 = start
    if(missing(maxiter)) maxiter=100
    conv_test=1
    iter=0
    while(conv_test > tol & iter<maxiter){
      pi = crossprod(P, pi1)
      conv_test = sqrt(sum(abs(pi-pi1)^2))
      # cat( conv_test, "\n")
      iter=iter+1
      pi1=pi
    }
    if(iter>=maxiter) warning("Max iter reached, tolerance achieved: ", conv_test)
    pi = pi / as.numeric(R%*%one)
  } else stop("method must be either 'lu' or 'limit'")
  ## normalize
  pi=pi/sum(pi)
  
  ## fix pi when there are "zero" cells
  
  if(n.0>0){
    pi.full=rep(0,n.0+n)
    pi.full[-zero.idx]=pi
    pi=pi.full
  }
  
  ## return stationary distribution
  
  as.numeric(pi)
}

