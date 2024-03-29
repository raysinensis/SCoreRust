% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_module_scores.R
\name{get_background}
\alias{get_background}
\title{Construct a Module/Pathway background matrix from gene list}
\usage{
get_background(background, features, bins, nrep)
}
\arguments{
\item{background}{calculated background expression matrix}

\item{features}{vector of genes in pathway}

\item{bins}{gene bin vector to use specific order}

\item{nrep}{number of repeats to calculate sampling of bins}
}
\value{
matrix of sampled background # of feature genes x cells
}
\description{
Construct a Module/Pathway background matrix from gene list
}
