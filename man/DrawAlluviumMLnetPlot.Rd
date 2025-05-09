% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualize_cell_communication.R
\name{DrawAlluviumMLnetPlot}
\alias{DrawAlluviumMLnetPlot}
\title{DrawAlluviumMLnetPlot}
\usage{
DrawAlluviumMLnetPlot(
  MLnetDir = MLnetDir,
  ImportDir = ImportDir,
  Sender,
  Receiver,
  ColorDB = colordb,
  Check = TRUE,
  top.n = 30,
  p_height = 7.5,
  p_width = 7
)
}
\arguments{
\item{MLnetDir}{Character, the path where the multilayer network is stored.}

\item{ImportDir}{Character, the path where the signals importance is stored.}

\item{Sender}{Character, sender cell type are shown.}

\item{Receiver}{Character, receiver cell type are shown.}

\item{ColorDB}{Vector, containing the color of cell types, key siganls and others.}

\item{Check}{Logical, Whether to do check visualization.}

\item{top.n}{Numercial, Only relevant if Check=T. Number of top regulation between LR pairs and target genes to show according to the signals importance.}

\item{p_height}{Numercial, the height of plot.}

\item{p_width}{Numercial, the width of plot.}
}
\description{
Draw Alluvium MLnet Plot
}
