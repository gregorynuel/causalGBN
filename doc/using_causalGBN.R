## -----------------------------------------------------------------------------
# load the library
# note that the rjbgoudie/parental package from github is loaded
# and produces the associated output
library(causalGBN)

## -----------------------------------------------------------------------------
# create a basic dag structure using the parental:bn function
dag1=bn(c(),c(1),c(1,2))

# or through the str2dag function
dag2=str2dag("[][1][1,2]")
identical(dag1,dag2)

# turn a dag into a string
dag2str(dag1)

# or into a dot description
dag2dot(dag1)

# which can be used to create a dag image
# if graphviz is installed on your system
suppressWarnings({
  check.dot=system("which dat",intern=TRUE)
})
if (length(check.dot)>0) {
  dag2dot(dag1,file="mydag.dot")
  system("dot -Grankdir=LR -Tpdf mydag.dot > mydag.pdf")
}

# we can also check if two DAGs belong to the same I-MEC
is_mec(bn(c(),1,2,2),bn(2,3,c(),2))
is_mec(bn(c(),1,2,2),bn(2,3,c(),2),I=list(2))
is_mec(bn(c(),1,2,2),bn(2,3,c(),2),I=list(4))

## -----------------------------------------------------------------------------
# and create a random GBN from a dag
dag=bn(c(),c(1),c(1,2))
set.seed(42)
ref=dag2ref(dag)

# a simulate some data with intervention sets
N=c(5,2,2)
I=list(numeric(0),c(1,3),2)
x=simul_gbn(ref$m,ref$sigma,ref$W,N,I)
x

# by default the intervention genes are filled with 0
# but we can provide something else
x=simul_gbn(ref$m,ref$sigma,ref$W,N,I,
            fill_int=function(n) runif(n,min=0,max=0.001))
x

## -----------------------------------------------------------------------------
# the first thing to do is to set up the workspace for efficient repeated estimations
global=compute_global(x,N,I)

fit=estim_gbn(dag)

plot(unlist(list(ref$m,ref$sigma,ref$W)),
     unlist(list(fit$m,fit$sigma,fit$W)),
     xlab="ref",ylab="fit",main=paste0(sum(N)," samples"))
abline(0,1,lty=2)

x=simul_gbn(ref$m,ref$sigma,ref$W,100*N,I)
global=compute_global(x,100*N,I)
fit=estim_gbn(dag)
plot(unlist(list(ref$m,ref$sigma,ref$W)),
     unlist(list(fit$m,fit$sigma,fit$W)),
     xlab="ref",ylab="fit",main=paste0(sum(100*N)," samples"))
abline(0,1,lty=2)

## -----------------------------------------------------------------------------
dag=bn(c(),1,2,2,3,4,5)
set.seed(42)
ref=dag2ref(dag)
N=c(200,5)
I=list(numeric(0),c(4))
x=simul_gbn(ref$m,ref$sigma,ref$W,N,I)
global=compute_global(x,N,I)
fit=mc3_gbn(verbose=FALSE)

post_cpdag=list(fit$post[1])
post_rep=list(str2dag(names(fit$post[1])))
for (i in 2:length(fit$post)) {
  dag=str2dag(names(fit$post[i]))
  post_prob=sapply(post_cpdag,function(z) z[[1]])
  idx=rep(FALSE,length(post_rep))
  # test MEC only for dag with very similar probability
  idx1=abs(fit$post[i]-post_prob)<1e-10
  if (sum(idx1)>0) {
    idx2=sapply(post_rep[idx1],function(z)is_mec(z,dag,I))
    idx[idx1]=idx2
  }
  if (sum(idx)==0) {
    post_cpdag=c(post_cpdag,list(fit$post[i]))
    post_rep=c(post_rep,list(str2dag(names(fit$post[i]))))
  } else {
    post_cpdag[[which(idx)]]=c(post_cpdag[[which(idx)]],fit$post[i])
  }
}

post_prob=sapply(post_cpdag,sum)
names(post_prob)=sapply(post_rep,dag2str)
idx=sort(post_prob,decreasing=TRUE,index.return=TRUE)$ix
post_prob=post_prob[idx]
post_cpdag=post_cpdag[idx]

post_prob[1:3]
post_cpdag[1:3]

