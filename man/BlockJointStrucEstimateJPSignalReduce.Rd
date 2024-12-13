% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BlockJointStrucEstimateJPSignalReduce.R
\name{BlockJointStrucEstimateJPSignalReduce}
\alias{BlockJointStrucEstimateJPSignalReduce}
\title{Estimate a Specific Joint Block Basis}
\usage{
BlockJointStrucEstimateJPSignalReduce(
  blockIn,
  datablock,
  dataname,
  VBars,
  UBars,
  phiBars,
  psiBars,
  rBars,
  outMap,
  theta0 = 45,
  optArgin = NULL,
  iprint = FALSE,
  figdir = ""
)
}
\arguments{
\item{blockIn}{A logical vector indicating shared blocks.}

\item{datablock}{List of data blocks.}

\item{dataname}{Vector of names of data blocks.}

\item{VBars}{A list of adjusted signal row spaces.}

\item{UBars}{List of adjusted signal column spaces for each data block.}

\item{phiBars}{A numeric vector of perturbation angles for each data matrix.}

\item{psiBars}{Numeric vector of loadings perturbation angles for each data matrix.}

\item{rBars}{A numeric vector of adjusted signal ranks for each data matrix.}

\item{outMap}{A mapping between block index sets and estimated partially shared structures.}

\item{theta0}{A numeric value representing the angle between estimated spaces and the optimized vector.}

\item{optArgin}{A list of optimization tuning parameters (optional). Includes:
\itemize{
  \item \code{tau0}: Initial tuning parameter.
  \item \code{tau_max}: Maximum tuning parameter.
  \item \code{mu}: Step size.
  \item \code{t_max}: Maximum number of iterations.
  \item \code{tol}: Tolerance level for optimization.
  \item \code{delta}: Perturbation parameter.
}}

\item{iprint}{Logical. Whether to save figures.}

\item{figdir}{Directory for saving figures. Used only when `iprint = TRUE`.}
}
\value{
A list containing:
\itemize{
  \item \code{Vi}: Estimated basis matrix.
  \item \code{curRanks}: Updated cumulative ranks for each data block.
}
}
\description{
This function is called internally within `DJIVEJointStrucEstimateJPLoadInfo`.
It estimates a specific joint block basis based on provided inputs.
It calculates an estimated basis matrix and updates cumulative ranks for data blocks.
}
