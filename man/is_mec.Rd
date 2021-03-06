% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dag_utils.R
\name{is_mec}
\alias{is_mec}
\title{tells if two DAGs are I-Markov equivalent or not}
\usage{
is_mec(dag1, dag2, I = list())
}
\arguments{
\item{dag1}{a DAG object}

\item{dag2}{a DAG object}

\item{I}{a list of intervention (in addition to emptyset which is always added)}
}
\value{
a boolean
}
\description{
tells if two DAGs are I-Markov equivalent or not
}
\examples{
is_mec(bn(c(),1,2,2),bn(2,3,c(),2))
is_mec(bn(c(),1,2,2),bn(2,3,c(),2),I=list(2))
is_mec(bn(c(),1,2,2),bn(2,3,c(),2),I=list(4))
dag1=bn(c(),c(1,3),c(1))
dag2=bn(c(2,3),c(3),c())
}
