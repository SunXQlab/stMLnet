% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_signal_importance.R
\name{getSiganlImport}
\alias{getSiganlImport}
\title{getSiganlImport}
\usage{
getSiganlImport(
  SiganlActivity,
  Lable,
  OutputDir = NULL,
  NCores = NULL,
  AutoPara = TRUE,
  NTrees = 500,
  NTrys = 10,
  TreeMethod = "variance",
  NodeSize = 5,
  NPert = 10
)
}
\arguments{
\item{SiganlActivity}{List, The sublist contains the activity of upstream signal pairs and expression of downstream target genes in specific celltype pair.}

\item{Lable}{Character, Denotes which celltype pair to study.}

\item{OutputDir}{Character, The output path of the currently running job where temporary and final results will be saved.}

\item{NCores}{Numercial, set the cores for the parallel process.}

\item{AutoPara}{Logical, Whether to do automatic optimization parameter.}

\item{NTrees}{Numercial, number of trees in random forests model, see Seurat package for more details.}

\item{NTrys}{Numercial, Only relevant if AutoPara=F. Number of variables to possibly split at in each node, see Seurat package for more details.}

\item{TreeMethod}{Character, Only relevant if AutoPara=F. Splitting rule, see Seurat package for more details.}

\item{NodeSize}{Numercial, Only relevant if AutoPara=F. Minimal node size, see Seurat package for more details.}

\item{NPert}{Numercial, Number of shuffle.}
}
\value{
List, The first sublist (df_im) is the importance of upstream signal pairs (Ligand/Receptor pairs), the second sublist (df_pim) is the importance of upstream signals (ligands/receptors).
}
\description{
calculate the upstream signal pairs or signals importance in the multilay signal network of cell communication
}
