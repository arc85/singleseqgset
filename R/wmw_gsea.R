#' Wilcoxon Mann Whitney Correlation Corrected GSEA
#'
#' The function performs a correlation corrected GSEA for each cluster and returns a list of correlation corrected statistics and p values
#'
#' @param expr.mat Expression matrix for the dataset
#'
#' @param cluster.cells List of clusters containing the assigned cells
#'
#' @param log.fc.cluster Log fold change in gene expression between clusters
#'
#'
#' @param gene.sets Named list of gene sets containing the genes defining each set
#'
#' @return Returns the log fold change in gene expression for each gene for cells in a cluster vs cells not in the cluster
#'
#' @export

wmw.gsea <- function(expr.mat,cluster.cells,log.fc.cluster,gene.sets)

{

results <- vector(mode="list",length=length(cluster.cells))

for (z in 1:length(cluster.cells)) {

gene.statistics = log.fc.cluster[[z]]

data <- t(expr.mat[,cluster.cells[[z]]])
n = nrow(data)
p = ncol(data)
gene.set.indexes = gene.sets

num.gene.sets = length(gene.set.indexes)
n = nrow(data)
p.values = matrix(0,nrow=num.gene.sets)
rownames(p.values) = names(gene.set.indexes)
gene.set.statistics = matrix(T,nrow=num.gene.sets)
rownames(gene.set.statistics) = names(gene.set.indexes)

for (i in 1:num.gene.sets) {

	x.test <- gene.statistics[!is.na(match(names(gene.statistics),gene.set.indexes[[i]]))]
	y.test <- gene.statistics[is.na(match(gene.statistics,gene.set.indexes[[i]]))]

	if (length(x.test)<5) {
		gene.set.statistics[i] <- NA
		p.values[i] <- 1
	} else {

	wilcox.results <- wilcox.test(x=x.test,y=y.test,alternative="two.sided",exact=F,correct=F)

	cor.mat <- cor(as.matrix(data[,names(x.test)]))
	cor.mat <- cor.mat[upper.tri(cor.mat,diag=F)]
	cor.mat <- cor.mat[which(!is.na(cor.mat))]
	mean.cor <- mean(cor.mat)

	m1 <- length(x.test)
	m2 <- length(y.test)

	rank.sum <- wilcox.results$statistic
	var.rank.sum <- ((m1*m2)/(2*pi))*(asin(1) + (m2 - 1)*asin(.5) + (m1-1)*(m2-1)*asin(mean.cor/2) +(m1-1)*asin((mean.cor+1)/2))
	z.stat = (rank.sum - (m1*m2)/2)/sqrt(var.rank.sum)
	gene.set.statistics[i] = z.stat
	lower.p = pnorm(z.stat,lower.tail=T)
	upper.p = pnorm(z.stat,lower.tail=F)
	p.values[i] = 2*min(lower.p,upper.p)
	}
}

results[[z]]$p.values <- p.values
results[[z]]$statistics <- gene.set.statistics

}

p.val.res <- data.frame(lapply(results,function(x) x[1]))
statistics.res <- data.frame(lapply(results,function(x) x[2]))
colnames(p.val.res) <- colnames(statistics.res) <- paste("Cluster_",seq(1,length(cluster.cells),1),sep="")

return(list("GSEA_statistics"=statistics.res,"GSEA_p_values"=p.val.res))

}
