# input: a real vector l. output: log(sum(exp(l)))
logsumexp=function(l) {
  i=which.max(l);
  res=l[i]+log1p(sum(exp(l[-i]-l[i])));
  if (is.nan(res)) res=-Inf;
  return(res);
}

#' function creating the workspace for fitting our models
#' the output must be named "global" and used as a global variable
#' @importFrom hashmap hashmap
#'
#' @param x the (n,p) real value (n experiments, p genes)
#' @param N vector of length k, number of replication per design
#' @param I list of length k, intervention pattern per design (\code{NULL} means "no intervention")
#' e.g. \code{N=c(20,5,5)} and \code{I=list(NULL,c(1,3),c(2,3))} => 20 WT + 5 KO on 1 and 3 + 5 KO on 2 and 3
#' @return a list containing
#' - x the data, p number of genes, n number of experiments
#' - idx a list of length p such that \code{idx[[j]]} contains all experiments with no intervention on j
#' - cache is a empty hashmap with numeric values, further used for speeding up computations
#' @export
compute_global=function(x,N,I) {
  p=ncol(x)
  tmp=matrix(FALSE,nrow=sum(N),ncol=p)
  pos=0
  for (i in 1:length(N)) {
    tmp[pos+1:N[i],I[[i]]]=TRUE
    pos=pos+N[i]
  }
  n=apply(!tmp,2,sum);
  idx=lapply(as.list(data.frame(!tmp)),which)
  names(idx)=NULL
  cache=hashmap(character(0),numeric(0))
  return(list(x=x,p=p,n=n,idx=idx,cache=cache))
}

#' simple function estimating the model for a given dag
#' this function use the global variable "global" which must be created by the compute_global function
#'
#' @importFrom stats lm
#'
#' @param dag a DAG
#' @return an estimated GBN
#' @export
estim_gbn=function(dag) {
  global=NULL; rm(list="global") # a fix to avoid global warnings
  if (length(dag)!=global$p) stop("Incompatible dag length !")
  reg=vector(global$p,mode="list")
  for (j in 1:global$p) {
    pa=dag[[j]]
    idx=global$idx[[j]]
    data=data.frame(global$x[idx,c(j,pa),drop=FALSE])
    colnames(data)=paste0("X",c(j,pa))
    reg[[j]]=lm(formula = as.formula(paste0("X", j, "~.")), data = data)
  }
  sigma=sapply(reg,function(z)summary(z)$sigma)
  m=sapply(reg,function(z)as.numeric(coefficients(z)[1]))
  W=matrix(0,global$p,global$p)
  for (j in 1:global$p) {
    W[dag[[j]],j]=coefficients(reg[[j]])[-1]
  }
  loglik=sapply(reg,logLik)
  dim=sapply(dag,length)+2
  return(list(loglik=loglik,dim=dim,m=m,sigma=sigma,W=W))
}


#' main function computing the local score i.e. the regression on a specific gene
#' this function use the global variable "global" which must be created by the compute_global function
#' the function check if the result is stored in global$cache and return or compute it
#'
#' @importFrom stats as.formula coefficients logLik
#'
#' @param j one gene index
#' @param dag a dag object (we only need the parents of j, but it is simpler to pass the whole dag)
#' @param crit either the extended BIC (default) or the BIC. The extended BIC should select more parsimonious models

#' @return the local score which is the likelihood of the j component of the model plus the penalty

#' @export
cached_local_score=function(j,dag,crit=c("eBIC","BIC")) {
  global=NULL; rm(list="global") # a fix to avoid global warnings
  crit=match.arg(crit)
  pa=dag[[j]]
  idx=global$idx[[j]]
  n=global$n[j]
  str=paste0(j,"|",paste0(sort(pa),collapse=","))
  res=global$cache[[str]]
  if (is.na(res)) {
    data=data.frame(global$x[idx,c(j,pa),drop=FALSE])
    colnames(data)=paste0("X",c(j,pa))
    reg=lm(formula = as.formula(paste0("X", j, "~.")), data = data)
    res=as.numeric(logLik(reg))
    global$cache[[str]]<<-res
  }
  if (crit=="BIC") {
    res=res-log(n)*(length(pa)+2)/2
  } else if (crit=="eBIC") {
    res=res-log(n)*(length(pa)+2)/2-lchoose(global$p-1,length(pa))
  }
  return(res)
}

# main function computing the score of a dag
# this function use the global variable "global" which must be created by the compute_global function
# uses cached_local_score (see above)
cached_dag_score=function(dag,crit=c("eBIC","BIC")) {
  global=NULL; rm(list="global") # a fix to avoid global warnings
  crit=match.arg(crit)
  res=sapply(1:global$p,cached_local_score,dag,crit)
  return(res)
}

# fast updating with only pa_j changed only for j in J
# this function use the global variable "global" which must be created by the compute_global function
# even faster update when the dag is slightly changing.
update_cached_dag_score=function(dag,res,J,crit=c("eBIC","BIC")) {
  crit=match.arg(crit)
  res[J]=sapply(J,cached_local_score,dag,crit)
  return(res)
}


#' The MC3 implementation using the cache mecanism
#' this function use the global variable "global" which must be created by the compute_global function
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom parental sampleBN
#'
#' @param itermax max number of iterations
#' @param burnin the number of iterations to discard at the end
#' @param initial a starting DAG or NULL for random start
#' @param maxParents the maximum number of parents for each node
#' @param constraintT a matrix of additional constraints (0: no constraint, 1: edge must be present, -1: must be absent)
#' @param verbose boolean for verbose output
#' @param exact boolean for exact posterior on the approximate support
#' @param crit criterion for the score
#'
#' @return the MC3 GBN estimate
#'
#' @export
mc3_gbn=function(itermax=5000,burnin=1000,initial=NULL,maxParents=global$p-1,
                        constraintT=matrix(0,global$p,global$p),verbose=TRUE,exact=TRUE,
                        crit=c("eBIC","BIC")) {
  global=NULL; rm(list="global") # a fix to avoid global warnings
  crit=match.arg(crit)
  dag=vector(itermax,mode="list")
  score=rep(NA,itermax)
  accepted=rep(NA,itermax)

  start=proc.time()
  if (verbose) cat("Running",itermax,"MC3 iterations ...\n")
  if (verbose) {
    pb=txtProgressBar(max=itermax,style=3)
    setTxtProgressBar(pb,0)
  }
  if (is.null(initial)) initial=sampleBN(n=global$p,maxNumberParents=maxParents)
  dag0=initial
  fit0=cached_dag_score(dag0,crit)
  routes0=routes(dag0)
  adj0=as.adjacency(dag0)
  logMoves0=logNumMHNeighbours(routes0,adj0,constraintT,maxParents)
  for (iter in 1:itermax) {
    # proposal move
    prop=proposal(routes0,adj0,constraintT,maxParents)
    # update proposal dag and aux data
    dag1=prop$dag
    routes1=prop$routes
    adj1=prop$adj
    logMoves1=logNumMHNeighbours(routes1,adj1,constraintT,maxParents)
    J=prop$changed
    # update fitting using the fast update function
    fit1=update_cached_dag_score(dag1,fit0,J,crit)
    # acceptance rate
    logAR=sum(fit1[J])-sum(fit0[J])-logMoves1+logMoves0

    #accept=runif(1)<min(1,exp(beta*logAR))
    accept=runif(1)<min(1,exp(logAR))

    if (accept) {
      dag0=dag1
      fit0=fit1
      routes0=routes1
      adj0=adj1
      logMoves0=logMoves1
    }
    accepted[iter]=accept
    dag[[iter]]=dag0
    score[iter]=sum(fit0)

    if (verbose) setTxtProgressBar(pb,value=iter)
  }
  elapsed=(proc.time()-start)[[3]]
  if (verbose) cat("\nRunning time",elapsed, "seconds\n")

  ratio=mean(accepted)

  post=sort(table(sapply(dag[burnin:itermax],dag2str))/(itermax-burnin+1),decreasing=TRUE)
  support=names(post)
  if (exact) {
    tmp=sapply(support,function(z) {sum(cached_dag_score(str2dag(z),crit))})
    post=exp(tmp-logsumexp(tmp))
  }
  return(list(dag=dag,score=score,accepted=accepted,ratio=ratio,
              support=lapply(support,str2dag),post=post))
}


