get.UD <- function(R,zero.idx=integer()){

    ## subset if needed

    n.0=length(zero.idx)
    if(length(zero.idx)>0){
        R=R[-zero.idx,-zero.idx]
    }

    ## get infinitessimal generator

    n=nrow(R)
    one=Matrix(1,n,1)
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
    ## normalize
    pi=pi/sum(pi)

    ## fix pi when there are "zero" cells

    if(n.0>0){
        pi.full=rep(0,n.0+n)
        pi.full[-zero.idx]=pi
        pi=pi.full
    }

    ## return stationary distribution
    
    pi
}

