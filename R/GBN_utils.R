# Computes the total effects L=(I-W)^-1
# from a direct effect matrix W
Lfct=function(W) {
  d=dim(W)[1];
  I=diag(rep(1,d));
  L=I; tmp=I; i=0;
  while (i<d-1) {
    tmp=tmp%*%W; L=L+tmp; i=i+1;
  }
  return(L);
}

# Computes the covariance matrix of a GBN
# from a vector of residual standard deviation
# and a matrix of total effects Cov = L^T diag(sd^2) L
Sfct=function(sd,L) {
  return(t(L)%*%diag(sd^2)%*%L);
}


#' Generates a reference GBN from a DAG
#'
#' @param dag a DAG object
#' @param w.range range for the direct effect (uniform distribution)
#' @param m.range range for the residual means (uniform distribution)
#' @param s.mean mean of the log of residual
#' @param s.sd sd of the log of residual
#' @return a GBN
#'
#' @importFrom stats rnorm runif
#' @importFrom parental as.adjacency
#'
#' @examples
#' dag=str2dag("[][1][1,2]")
#' dag2ref(dag)
#'
#' @export
dag2ref=function(dag,w.range=c(-3,3),m.range=c(-2,2),
                 s.mean=log(1),s.sd=0.3) {
  p=length(dag)
  W=as.adjacency(dag)
  W[W==1]=runif(sum(W==1),min=w.range[1],max=w.range[2])
  s=rnorm(p,mean=s.mean,sd=s.sd)
  L=Lfct(W)
  Sigma=Sfct(exp(s),L)
  m=runif(p,min=m.range[1],max=m.range[2])
  mu=(m%*%L)[1,]
  return(list(dag=dag,m=m,sigma=exp(s),s=s,W=W,
              L=L,mu=mu,Sigma=Sigma,p=length(m)))
}

#' Simulate a GBN with a specific intervention design
#'
#' @importFrom MASS mvrnorm
#'
#' @param m vector of residual means
#' @param sigma vector of residual standard deviations
#' @param W matrix of direct effects
#' @param N either an integer or a vector of integers
#' @param I a list of intervention sets
#' @param fill_int a function of \code{n} returning a vector of \code{n} data to fill the intervention genes
#'
#' @return the simulated data matrix with \code{sum(N)} rows
#'
#' @examples
#' dag=bn(c(),1,2)
#' set.seed(42)
#' ref=dag2ref(dag)
#' N=c(50,5,5)
#' I=list(numeric(0),c(1,3),2)
#' x=simul_gbn(ref$m,ref$sigma,ref$W,N,I)
#' x=simul_gbn(ref$m,ref$sigma,ref$W,N,I,fill_in=function(n)rnorm(n))
#'
#' @export
simul_gbn=function(m,sigma,W,N,I=list(c()),fill_int=function(n){rep(0,n)}) {
  if (length(I)!=length(N)) stop("N and I mismatch.")
  p=length(m)
  if (length(sigma)!=p) stop("m and sigma mismatch")
  if (nrow(W)!=p) stop("m and W mismatch")
  if (ncol(W)!=p) stop("m and W mismatch")
  res=matrix(NA,sum(N),p)
  for (j in 1:length(N)) {
    localW=W
    localsigma=sigma
    int=I[[j]]
    localsigma[int]=0
    localW[,int]=0
    L=Lfct(localW)
    Sigma=Sfct(localsigma,L)
    nu=matrix(rep(m,N[j]),nrow=N[j],byrow=TRUE)
    if (length(int)>0) nu[,int]=fill_int(N[j]*length(int))
    mu=nu%*%L
    x=mvrnorm(N[j],mu=rep(0,p),Sigma=Sigma)+mu
    res[c(0,cumsum(N))[j]+1:N[j],]=x
  }
  return(res)
}

