% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_signal_activity.R
\name{getSiganlLinks}
\alias{getSiganlLinks}
\title{getSiganlLinks}
\usage{
getSiganlLinks(mulNet)
}
\arguments{
\item{mulNet}{List, the multilayer signal network between of specific celltype pair.}
}
\value{
List. The sublist contains the upstream signal pairs of specific target genes.
}
\description{
get the upstream signal pairs of downstream target genes in the multilayer signal network
}
