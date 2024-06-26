% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_cell_communication.R
\name{runstMLnet}
\alias{runstMLnet}
\title{Title}
\usage{
runstMLnet(
  ExprMat,
  AnnoMat,
  DistMat,
  LigClus = NULL,
  RecClus = NULL,
  Normalize = F,
  NormMethod = "LogNormalze",
  logfc.ct = 0.1,
  pct.ct = 0.05,
  pval.ct = 0.05,
  expr.ct = 0.1,
  OutputDir = NULL,
  Databases = NULL,
  RecTF.method = "Fisher",
  TFTG.method = "Fisher",
  TGList = NULL,
  LigList = NULL,
  RecList = NULL,
  NCores = 6,
  AutoPara = TRUE,
  NTrees = 500,
  NTrys = 10,
  TreeMethod = "variance",
  NodeSize = 5,
  NPert = 10
)
}
\arguments{
\item{ExprMat}{Matrix, The expression matrix with rows as genes (gene symbols) and columns as cells.}

\item{AnnoMat}{Matrix, The annotation matrix for cell type annotation. The first column is barcode and the second is cell type.}

\item{DistMat}{Matrix, The distance matrix represents the distance between two cells (e.g., Euclidean distance).}

\item{LigClus}{Vector, cell types that sending signals.}

\item{RecClus}{Vector, cell types that receiving signals.}

\item{Normalize}{Logical, Whether to do normalization.}

\item{NormMethod}{Character, Denotes which normalization method to use. Available options are: LogNormalze, SCTransform. See Seurat package for more details.}

\item{logfc.ct}{Numercial, Screening threshold for FindMarkers in Seurat. The default setting is 0.1.}

\item{pct.ct}{Numercial, Screening threshold for FindMarkers in Seurat. The default setting is 0.05.}

\item{pval.ct}{Numercial, Screening threshold for FindMarkers in Seurat. The default setting is 0.05.}

\item{expr.ct}{Numercial, Screening threshold for high expressed gene in groups of cells. Default is 0.05.}

\item{OutputDir}{Character, The output path of the currently running job where temporary and final results will be saved.}

\item{Databases}{List, The prior database used by running jobs for now. Databases includes Ligand-Receptor interactions (LigRec.DB), Receptor-TF interactions (RecTF.DB) and TF-Target interactions (TFTG.DB).}

\item{RecTF.method}{Character, Denotes the strategy for filtering downstream pairing signals (Receptor-TF).  Available options are: Fisher(default, meaning Fisher exact test) and Search (meaning searching in database).}

\item{TFTG.method}{Character, Denotes the strategy for filtering downstream pairing signals (TF-Target).  Available options are: Fisher(default, meaning Fisher exact test) and Search (meaning searching in database).}

\item{TGList}{List, The target genes of interest in groups of cells (RecClus).}

\item{LigList}{List, The potential ligands in groups of cells (LigClus).}

\item{RecList}{List, The potential receptors in groups of cells (RecClus).}

\item{NCores}{Numercial, set the cores for the parallel process.}

\item{AutoPara}{Logical, Whether to do automatic optimization parameter.}

\item{NTrees}{Numercial, number of trees in random forests model, see Seurat package for more details.}

\item{NTrys}{Numercial, Only relevant if AutoPara=F. Number of variables to possibly split at in each node, see Seurat package for more details.}

\item{TreeMethod}{Character, Only relevant if AutoPara=F. Splitting rule, see Seurat package for more details.}

\item{NodeSize}{Numercial, Only relevant if AutoPara=F. Minimal node size, see Seurat package for more details.}

\item{NPert}{Numercial, Number of shuffle.}
}
\value{
List, containning the output folder of the multilayer signal network (MLnetDir), the signal activity (ActivityDir) and the signal importance (ImportDir).
}
\description{
Title
}
