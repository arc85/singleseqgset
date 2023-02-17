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

if (!inherits(expr.mat, "dgCMatrix")) {
#if (!class(expr.mat)=="dgCMatrix") {
	expr.mat <- Matrix(expr.mat,sparse=TRUE)
}

cluster.cells <- vector(mode="list",length=length(levels(cluster.ids)))
names(cluster.cells) <- levels(cluster.ids)

for (i in names(cluster.cells)) {
	cluster.cells[[i]] <- colnames(expr.mat)[which(cluster.ids==i)]
}

log.fc.cluster <- vector(mode="list",length=length(cluster.cells))
names(log.fc.cluster) <- names(cluster.cells)

sym_diff <- function(a,b) { unique(setdiff(a,b),setdiff(b,a)) }

for (i in 1:length(cluster.cells)) {

	other.clusters <- sym_diff(names(cluster.cells),names(cluster.cells)[i])

	mean.cluster.exp <- apply(
		expr.mat[,cluster.cells[[i]]],
		1, function(x) {
			log(mean(expm1(x))+1)
	})

	mean.exp.other <- apply(
		expr.mat[,c(unlist(cluster.cells[other.clusters]))],
		1, function(x) {
			log(mean(expm1(x))+1)
	})


	log.fc.cluster[[i]] <- mean.cluster.exp-mean.exp.other

}

return(list("cluster.cells"=cluster.cells,"log.fc.cluster"=log.fc.cluster))

}
