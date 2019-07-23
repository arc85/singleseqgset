#' Calculate metric for GSEA
#'
#' This function calculates the metric to be used for the gene set enrichment analysis. The most basic metric is the log fold-change in gene expression between a selected cluster and all other clusters, and is implemented here. Other possibile metrics could be used.
#'
#' @param cluster.ids Named vector of the cluster identities assigned for each cell
#'
#' @param expr.mat Expression matrix with rows as genes and columns as cells
#'
#' @return Returns the log fold change in gene expression for each gene for cells in a cluster vs cells not in the cluster
#'
#' @export

logFC <- function(cluster.ids,expr.mat) {

if (is.factor(cluster.ids)==F) {
  cluster.ids <- as.factor(cluster.ids)
}

if (!class(expr.mat)=="dgCMatrix") {
	expr.mat <- Matrix(expr.mat,sparse=TRUE)
}

cluster.cells <- vector(mode="list",length=length(levels(cluster.ids)))
names(cluster.cells) <- levels(cluster.ids)

for (i in names(cluster.cells)) {
	cluster.cells[[i]] <- colnames(expr.mat)[which(cluster.ids==i)]
}

mean.cluster.exp <- vector(mode="list",length=length(cluster.cells))
cluster.size <- vector(length=length(mean.cluster.exp))
names(mean.cluster.exp) <- names(cluster.size) <- levels(cluster.ids)

for (i in 1:length(cluster.cells)) {
	mean.cluster.exp[[i]] <- apply(expr.mat[,cluster.cells[[i]]],1,mean)
  cluster.size[i] <- ncol(expr.mat[,cluster.cells[[i]]])
}

sym_diff <- function(a,b) unique(setdiff(a,b),setdiff(b,a))

log.fc.cluster <- vector(mode="list",length=length(cluster.cells))
names(log.fc.cluster) <- names(cluster.cells)

for (i in 1:length(cluster.cells)) {
  mean.clust <- mean.cluster.exp[[i]]
	non.clust.names <- sym_diff(names(cluster.cells),names(cluster.cells)[i])
  non.clust.means <- vector("list",length=length(non.clust.names))

  for (z in 1:length(non.clust.names)) {
    non.clust.means[[z]] <- mean.cluster.exp[[non.clust.names[z]]]/cluster.size[non.clust.names[z]]
  }

	log.fc.cluster[[i]] <- mean.clust-Reduce("+",non.clust.means)/length(non.clust.means)
}

return(list("cluster.cells"=cluster.cells,"log.fc.cluster"=log.fc.cluster))

}
