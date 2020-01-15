#' convert a dag object into a string
#'
#' @param dag a dag object
#' @return a string
#' @examples
#' dag2str(bn(c(),c(1),c(1,2)))
#'
#' @export
dag2str=function(dag) {
  res=""
  for (j in 1:length(dag)) res=paste0(res,"[",paste(dag[[j]],collapse=","),"]")
  return(res)
}

#' convert a string into a dag object
#'
#' @param str a string (ex: \code{"[][2][3]"})
#' @return a dag object
#' @examples
#' str2dag("[][1][1,2]")
#'
#'
#' @export
str2dag=function(str) {
  str=gsub("\\[","",str)
  str=strsplit(str,"\\]")[[1]]
  str=lapply(str,function(z) as.numeric(strsplit(z,",")[[1]]))
  str <- lapply(str, as.integer)
  str <- lapply(str, sort.int)
  class(str) <- c("bn","parental")
  return(str)
}

#' convert a dag object into dot file
#'
#' @param dag a dag object
#' @param file a file for the output
#' @param w optional weights for the edges
#' @return a dot description of the dag
#' @examples
#' dag2dot(bn(c(),c(1),c(1,2)))
#'
#' @export
dag2dot=function(dag,file="",w=NULL) {
  content="digraph {\n"
  for (j in 1:length(dag)) {
    content=paste0(content,j,";\n")
    if (is.null(w)) {
      for (i in dag[[j]]) content=paste0(content,i,"->",j,";\n")
    } else {
      for (i in dag[[j]]) content=paste0(content,i,"->",j,"[label=\"",signif(w[i,j],3),"\"];\n")
    }
  }
  content=paste0(content,"}\n")
  cat(content,file=file)
}

# returns a list of the v-struct of a dag
# (c,a,b) such that a -> c <- b with no edge between a and b
v_struct=function(dag) {
  res=list()
  for (i in which(sapply(dag,length)>1)) {
    for (j in 1:(length(dag[[i]])-1)) for (k in (j+1):length(dag[[i]])) {
      a=dag[[i]][j]; b=dag[[i]][k]
      if (sum(is.element(dag[c(a,b)[c(a,b)>0]],c(a,b)))==0) res=c(res,list(c(i,a,b)))
    }
  }
  return(res)
}

#' tells if two DAGs are I-Markov equivalent or not
#'
#' @param dag1 a DAG object
#' @param dag2 a DAG object
#' @param I a list of intervention (in addition to emptyset which is always added)
#' @return a boolean
#'
#' @examples
#' is_mec(bn(c(),1,2,2),bn(2,3,c(),2))
#' is_mec(bn(c(),1,2,2),bn(2,3,c(),2),I=list(2))
#' is_mec(bn(c(),1,2,2),bn(2,3,c(),2),I=list(4))
#' dag1=bn(c(),c(1,3),c(1))
#' dag2=bn(c(2,3),c(3),c())
#' @export
is_mec=function(dag1,dag2,I=list()) {
  idx=sapply(I,length)>0
  I=I[idx]
  if (length(dag1)!=length(dag2)) stop("dag1 and dag2 must be defined on the same vertices")
  to1=to2=vector(length(dag1),mode="list")
  for (i in 1:length(dag1)) {
    to1[[i]]=which(sapply(dag1,function(z) is.element(i,z)))
    to2[[i]]=which(sapply(dag2,function(z) is.element(i,z)))
  }
  # check skeletton
  vois1=lapply(Map('union',dag1,to1),sort)
  vois2=lapply(Map('union',dag2,to2),sort)
  if (!identical(vois1,vois2)) return(FALSE)
  # add intervention edges
  if (length(I)>0) {
    for (i in 1:length(I)) {
      for (j in I[[i]]) {
        dag1[[j]]=c(-i,dag1[[j]])
        dag2[[j]]=c(-i,dag2[[j]])
      }
    }
    dag1=lapply(dag1,sort)
    dag2=lapply(dag2,sort)
  }
  return(identical(v_struct(dag1),v_struct(dag2)))
}

