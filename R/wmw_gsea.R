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
#' @param gene.sets Named list of gene sets containing the genes defining each set
#'
#' @return Returns the z scores and p values for gene set enrichment
#'
#' @export

wmw_gsea <- function(expr.mat,cluster.cells,log.fc.cluster,gene.sets)

{

if (!class(expr.mat)=="matrix") {
	expr.mat <- as.matrix(expr.mat)
}

results <- vector(mode="list",length=length(cluster.cells))

for (z in 1:length(cluster.cells)) {

gene.statistics <- log.fc.cluster[[z]]

data <- as.matrix(t(expr.mat[,cluster.cells[[z]]]))
n <- nrow(data)
p <- ncol(data)
gene.set.indexes = gene.sets

num.gene.sets <- length(gene.set.indexes)
p.values <- matrix(1,nrow=num.gene.sets)
rownames(p.values) <- names(gene.set.indexes)
gene.set.statistics <- matrix(0,nrow=num.gene.sets)
rownames(gene.set.statistics) <- names(gene.set.indexes)

for (i in 1:num.gene.sets) {

	x.test <- gene.statistics[!is.na(match(names(gene.statistics),gene.set.indexes[[i]]))]
	y.test <- gene.statistics[is.na(match(gene.statistics,gene.set.indexes[[i]]))]

	if (length(x.test)<5) {
		gene.set.statistics[i] <- 0
		p.values[i] <- 1
	} else {

		data.sub <- as.matrix(data[,names(x.test)])
		data.values <- apply(data.sub,2,sum)
		data.val.pos <- length(which(data.values>0))

		if (data.val.pos<5) {
			gene.set.statistics[i] <- 0
			p.values[i] <- 1
		} else {

			data.val.index <- which(data.values>0)

			wilcox.results <- wilcox.test(x=x.test,y=y.test,alternative="two.sided",exact=F,correct=F)

			cor.mat <- cor(as.matrix(data.sub[,data.val.index]))
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
	}

results[[z]]$p.values <- p.values
results[[z]]$statistics <- gene.set.statistics

}

p.val.res <- data.frame(lapply(results,function(x) x[1]))
statistics.res <- data.frame(lapply(results,function(x) x[2]))
colnames(p.val.res) <- colnames(statistics.res) <- names(cluster.cells)

return(list("GSEA_statistics"=statistics.res,"GSEA_p_values"=p.val.res))

}
