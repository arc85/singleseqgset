% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wmw_gsea.R
\name{wmw_gsea}
\alias{wmw_gsea}
\title{Wilcoxon Mann Whitney Correlation Corrected GSEA}
\usage{
wmw_gsea(expr.mat, cluster.cells, log.fc.cluster, gene.sets)
}
\arguments{
\item{expr.mat}{Expression matrix for the dataset}

\item{cluster.cells}{List of clusters containing the assigned cells}

\item{log.fc.cluster}{Log fold change in gene expression between clusters}

\item{gene.sets}{Named list of gene sets containing the genes defining each set}
}
\value{
Returns the z scores and p values for gene set enrichment
}
\description{
The function performs a correlation corrected GSEA for each cluster and returns a list of correlation corrected statistics and p values
}
