#' Estimate Joint and Partially Joint Structure
#'
#' This function establishes a Difference of Convex (DC) programming problem to estimate each partially
#' joint structure. It uses the Penalty Convex-Concave Procedure (CCP) algorithm to solve the DC programming problem.
#'
#' @param VBars List of adjusted signal row spaces.
#' @param UBars List of adjusted signal column spaces for each data block.
#' @param phiBars Numeric vector of perturbation angles for each data matrix.
#' @param psiBars Numeric vector of loadings perturbation angles for each data matrix.
#' @param rBars A numeric vector of adjusted signal ranks for each data matrix.
#' @param datablock List of data blocks.
#' @param iprint Logical. Whether to print results.
#' @param figdir Directory for saving figures. Set to `NULL` if do not want to save.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{outMap}: A mapping between joint block index sets and estimated partially shared structures.
#'   \item \code{keyIdxMap}: A mapping between joint block indices and data blocks.
#'   \item \code{jointBlockOrder}: A character vector recording the order of the joint blocks.
#' }
#'
#' @details
#' The function solves the DC programming problem iteratively, leveraging the Penalty CCP algorithm
#' to optimize the estimation of joint and partially joint structures. This approach is particularly
#' useful for complex multi-block data structures.
#'
#'
#' @export
DJIVEJointStrucEstimateJPLoadInfo <- function(
    VBars, UBars, phiBars, psiBars, rBars, datablock, iprint = FALSE, figdir = ""
) {

  # Check input dimensions
  if(!methods::is(datablock, "list")) stop("Input datablock has to be a list with length >= 2.")
  if(length(datablock) < 2) stop("Input datablock has to be a list with length >= 2.")
  if(max(sapply(datablock, ncol)) != min(sapply(datablock, ncol))) stop("All data blocks have to have the same number of columns (samples).")

  # Check data block names
  nb <- length(datablock)
  dataname <- names(datablock)
  if(is.null(dataname)){
    warning("Input datablock is unnamed, generic names for data blocks generated.")
    dataname <- paste0("Datablock", 1:nb)
  }

  theta0 = 45
  optArgin = list()


  allIdx <- 1:nb

  outMap <- list()
  keyIdxMap <- list()
  anglesMap <- list()
  jointBlockOrder <- character(0)
  # flag <- FALSE


  # TODO: can limit searching via flag (see also MATLAB code for curRanks)
  for (len in nb:1) {
    lenIdces <- utils::combn(allIdx, len) #, simplify = FALSE
    nlen <- ncol(lenIdces)

    for (i in 1:nlen) {
      blockIdx <- lenIdces[, i]
      blockIn <- allIdx %in% blockIdx

      Vi_angles <- BlockJointStrucEstimateJPSignalReduce(
        blockIn, datablock, dataname, VBars, UBars, phiBars, psiBars, rBars,
        outMap, theta0, optArgin, iprint, figdir
      )

      Vi <- Vi_angles[[1]]
      angles <- Vi_angles[[2]]


      if (ncol(Vi) > 0) {
        t <- Idx2numMJ(blockIn)
        outMap[[as.character(t)]] <- Vi
        keyIdxMap[[as.character(t)]] <- blockIdx
        anglesMap[[as.character(t)]] <- angles
        jointBlockOrder <- c(jointBlockOrder, as.character(t))
      }
    }


    # if (flag) {
    #   break
    # }
  }

  return(list(outMap = outMap, keyIdxMap = keyIdxMap, anglesMap = anglesMap, jointBlockOrder = jointBlockOrder))
}
