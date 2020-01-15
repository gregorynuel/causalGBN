# functions extracted from structmcmc

# @param constraint A matrix of dimension ncol(data) x ncol(data) giving
#   constraints to the sample space. The (i, j) element is:
#   \describe{
#     \item{1}{if the edge i -> j is required}
#     \item{-1}{if the edge i -> is excluded.}
#     \item{0}{if the edge i -> j is not constrained.}
#   }

transposeEdgeIsTogglable=function (routes, adjacency, constraintT, maxNumberParents = Inf)
{
  removable <- transposeEdgeIsRemovable(routes = routes, adjacency = adjacency,
                                        constraintT = constraintT)
  addable <- transposeEdgeIsAddable(routes = routes, adjacency = adjacency,
                                    constraintT = constraintT, maxNumberParents = maxNumberParents)
  removable + addable == T
}
transposeEdgeIsRemovable=function (routes, adjacency, constraintT)
{
  routes == 0 & t(adjacency) == 1 & constraintT == 0
}
transposeEdgeIsAddable=function (routes, adjacency, constraintT, maxNumberParents)
{
  addable <- routes == 0 & t(adjacency) == 0 & constraintT ==
    0
  canAddMoreParents <- colSums(adjacency) < maxNumberParents
  addable[!canAddMoreParents, ] <- F
  addable
}
edgeIsFlippable=function (routes, adjacency, constraintT, maxNumberParents)
{
  flippable <- routes == 1 & adjacency == 1 & constraintT ==
    0
  canAddMoreParents <- colSums(adjacency) < maxNumberParents
  flippable[!canAddMoreParents, ] <- F
  flippable
}
routes=function (x)
{
  stopifnot("bn" %in% class(x))
  nNodes <- nNodes(x)
  nodesSeq <- seq.int(nNodes)
  routes <- matrix(0, nNodes, nNodes)
  diag(routes) <- 1
  for (head in nodesSeq) {
    for (tail in x[[head]]) {
      routes <- routes + outer(routes[, tail], routes[head,
                                                      ])
    }
  }
  routes
}
routesRemoveEdge=function (x, i, j)
{
  #x - x[, i] %*% .Internal(t.default((x[j, ])))
  x-tcrossprod(x[,i],x[j,])
}
routesAddEdge=function (x, i, j)
{
  #x + x[, i] %*% .Internal(t.default((x[j, ])))
  x+tcrossprod(x[,i],x[j,])
}
logNumMHNeighbours=function (routes, adjacency, constraintT, maxNumberParents = Inf)
{
  removable <- transposeEdgeIsRemovable(routes = routes, adjacency = adjacency,
                                        constraintT = constraintT)
  addable <- transposeEdgeIsAddable(routes = routes, adjacency = adjacency,
                                    constraintT = constraintT, maxNumberParents = maxNumberParents)
  flippable <- edgeIsFlippable(routes = routes, adjacency = adjacency,
                               constraintT = constraintT, maxNumberParents = maxNumberParents)
  log(length(adjacency[flippable | addable | removable]))
}

proposal=function(routes0,adj0,constraintT,maxParents) {
  canAddOrRemove=transposeEdgeIsTogglable(routes0,adj0,constraintT,maxParents)
  canFlip=edgeIsFlippable(routes0,adj0,constraintT,maxParents)

  nonCycleInducing=which(canAddOrRemove, arr.ind = T)
  nonCycleInducingFlips=which(canFlip, arr.ind = T)

  nNonCycleInducing=nrow(nonCycleInducing)
  nNonCycleInducingFlips=nrow(nonCycleInducingFlips)
  numberOfMoves=nNonCycleInducing + nNonCycleInducingFlips

  select=sample.int(numberOfMoves,size=1)
  adj1=adj0
  routes1=routes0
  if (select <= nNonCycleInducing){
    # changing edge status (careful with the transposition)
    j=nonCycleInducing[[select, 1]]
    i=nonCycleInducing[[select, 2]]
    # only pa[[j]] are changed
    changed=j
    if (adj1[i,j]==1) {
      adj1[i,j]=0
      routes1=routesRemoveEdge(routes1,i,j)
    } else {
      adj1[i,j]=1
      routes1=routesAddEdge(routes1,i,j)
    }
  } else {
    # flipping edge
    i=nonCycleInducingFlips[[select-nNonCycleInducing, 1]]
    j=nonCycleInducingFlips[[select-nNonCycleInducing, 2]]
    # both pa[[i]] and pa[[j]] are changed
    changed=c(i,j)
    adj1[i,j]=0
    adj1[j,i]=1
    routes1=routesRemoveEdge(routes1,i,j)
    routes1=routesAddEdge(routes1,j,i)
  }
  return(list(dag=parental::as.bn(adj1),adj=adj1,routes=routes1,changed=changed))
}
