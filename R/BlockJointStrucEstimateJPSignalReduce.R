#' Estimate a Specific Joint Block Basis
#'
#' This function is called internally within `DJIVEJointStrucEstimateJPLoadInfo`.
#' It estimates a specific joint block basis based on provided inputs.
#' It calculates an estimated basis matrix and updates cumulative ranks for data blocks.
#'
#'
#' @param blockIn A logical vector indicating shared blocks.
#' @param datablock List of data blocks.
#' @param dataname Vector of names of data blocks.
#' @param VBars A list of adjusted signal row spaces.
#' @param UBars List of adjusted signal column spaces for each data block.
#' @param phiBars A numeric vector of perturbation angles for each data matrix.
#' @param psiBars Numeric vector of loadings perturbation angles for each data matrix.
#' @param rBars A numeric vector of adjusted signal ranks for each data matrix.
#' @param outMap A mapping between block index sets and estimated partially shared structures.
#' @param theta0 A numeric value representing the angle between estimated spaces and the optimized vector.
#' @param optArgin A list of optimization tuning parameters (optional). Includes:
#' \itemize{
#'   \item \code{tau0}: Initial tuning parameter.
#'   \item \code{tau_max}: Maximum tuning parameter.
#'   \item \code{mu}: Step size.
#'   \item \code{t_max}: Maximum number of iterations.
#'   \item \code{tol}: Tolerance level for optimization.
#'   \item \code{delta}: Perturbation parameter.
#' }
#' @param iprint Logical. Whether to save figures.
#' @param figdir Directory for saving figures. Used only when `iprint = TRUE`.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{Vi}: Estimated basis matrix.
#'   \item \code{curRanks}: Updated cumulative ranks for each data block.
#' }
#'
BlockJointStrucEstimateJPSignalReduce <- function(
    blockIn, datablock, dataname,
    VBars, UBars, phiBars, psiBars,
    rBars, outMap, theta0 = 45,
    optArgin = NULL, iprint = FALSE, figdir = ""
) {
  # if (is.logical(blockIn)) {
  #   blockIn <- which(blockIn)
  # }

  nb <- length(blockIn)
  allIdx <- 1:nb
  blockIdx <- allIdx[blockIn]
  blockName <- dataname[blockIn]
  cat(paste("Find joint structure shared only among", paste(blockName, collapse = ", "), ".\n"))

  n <- nrow(VBars[[1]])
  blockLen <- sum(blockIn)
  curBasisSizes <- rBars

  Vorth <- NULL
  Vnorth <- NULL
  if (length(Vorth) == 0 || is.null(Vorth)) {
    Vorth <- matrix(0, n, 1)
  }

  for (len in nb:(blockLen + 1)) {
    if (len > length(allIdx)) next
    lenIdces <- utils::combn(allIdx, len)
    for (i in 1:ncol(lenIdces)) {
      bkIdx <- lenIdces[, i]
      bkIn <- allIdx %in% bkIdx
      t <- Idx2numMJ(bkIn)

      if (!is.element(as.character(t), names(outMap))) next

      if (all(blockIdx %in% bkIdx)) {
        Vorth <- cbind(Vorth, outMap[[as.character(t)]])
      } else {
        Vnorth <- cbind(Vnorth, outMap[[as.character(t)]])
      }
    }
  }

  #projOrth <- qr.Q(qr(Vorth))
  #projOrth <- svd(Vorth)$u

  if (all(Vorth == 0)) {
    projOrth <- matrix(0, n, 0)
  } else {
    projOrth <- svd(Vorth)$u
  }
  VorthDim <- ncol(projOrth)

  Mo1 <- NULL
  Mo2 <- NULL
  Qc1 <- vector("list", nb)
  Qc2 <- vector("list", nb)
  Qc1Load <- vector("list", nb)
  Qc2Load <- vector("list", nb)
  precompute_load_matrix <- vector("list", nb)

  for (ib in seq_len(nb)) {
    d <- nrow(UBars[[ib]])
    precompute_load_matrix[[ib]] <- t(datablock[[ib]]) %*% UBars[[ib]] %*% t(UBars[[ib]]) %*% datablock[[ib]]

    if (blockIn[ib]) {
      Mo2 <- cbind(Mo2, VBars[[ib]])
      Qc1[[ib]] <- diag(n)

      curVBars <- VBars[[ib]] %*% pracma::nullspace(t(Vorth) %*% VBars[[ib]])
      # if (is.null(curVBars) || ncol(curVBars) == 0) {
      #   warning(paste("No valid curVBars found for datablock", ib))
      #   next
      # }
      # print(curVBars[1:min(5, nrow(curVBars)), 1:min(5, ncol(curVBars))])

      curBasisSizes[ib] <- ncol(curVBars)

      Qc2[[ib]] <- curVBars %*% t(curVBars) / cos(phiBars[ib] * pi / 180)^2
      Qc1Load[[ib]] <- t(datablock[[ib]]) %*% datablock[[ib]]
      Qc2Load[[ib]] <- precompute_load_matrix[[ib]] / cos(psiBars[ib] * pi / 180)^2
    } else {
      Qc1[[ib]] <- VBars[[ib]] %*% t(VBars[[ib]]) / cos(phiBars[ib] * pi / 180)^2
      Qc2[[ib]] <- diag(n)
      Qc1Load[[ib]] <- matrix(0, n, n)
      Qc2Load[[ib]] <- matrix(0, n, n)
    }
  }

  # print('Qc1')
  # print(Qc2[[3]])
  #
  Qo1 <- if (is.null(Mo1)) 1e-6 * matrix(1, n, n) else Mo1 %*% t(Mo1)
  Qo2 <- if (is.null(Mo2)) 1e-6 * matrix(1, n, n) else Mo2 %*% t(Mo2)

  if (!is.null(Vnorth)) {
    Qc1[[length(Qc1) + 1]] <- Vnorth %*% t(Vnorth) / cos(theta0)^2
    Qc2[[length(Qc2) + 1]] <- diag(n)
  }

  # print('Qc2')
  # print(Qc3)
  #
  searchNext <- TRUE

  if (blockLen == 1) {
    if (VorthDim > 0) {
      V0 <- svd(datablock[[blockIn]] %*% (diag(nrow(Vorth)) - projOrth %*% t(projOrth)))$v[, 1:max(rBars[blockIn])]
    } else {
      V0 <- svd(datablock[[blockIn]])$v[, 1:max(rBars[blockIn])]
    }
  } else {
    V0 <- svd(Qo2)$v[, 1:max(rBars[blockIn])]
  }

  # print(V0)

  Vi <- NULL
  angles <- matrix(0, nrow = nb, ncol = ncol(V0))# matrix(90, nb, ncol(V0))
  j <- 0
  #print(angles)

  while (searchNext) {
    j <- j + 1
    if (j > max(rBars[blockIn])) {
      break
    }
    cat(sprintf("Search Direction %d:\n", j))
    Vorth <- cbind(Vorth, Vi)

#################### 之前全对
    for (ib in seq_len(nb)) {
      if (blockIn[ib]) {
        #curVBars <- VBars[[ib]] %*% qr.Q(qr(t(Vorth) %*% VBars[[ib]]))
        curVBars <- VBars[[ib]] %*% pracma::nullspace(t(Vorth) %*% VBars[[ib]])
        curBasisSizes[ib] <- ncol(curVBars)
        Qc2[[ib]] <- curVBars %*% t(curVBars) / cos(phiBars[ib] * pi / 180)^2
      }
    }

    # print(Qc2)

    output <- penaltyCCPJPEarlyStopLoadInfo(V0[, j], Qo1, Qo2, Qc1, Qc2, Qc1Load, Qc2Load, Vorth, optArgin)
    opt_v <- output[[1]]
    cache_v <- output[[2]]
    converge <- output[[5]]

    if (converge <= 1) {
      angleHats <- ccpOutAnalysisMJ(cache_v, VBars)

      figname <- paste(paste(blockName, collapse = "-"), "-joint-optV", j, sep = "")
      ccpOutVisualMJ(angleHats, phiBars, dataname, iprint, figdir, figname)
      for (ib in seq_len(nb)) {
        angles[ib, j] <- angleHats[[ib]][length(angleHats[[ib]])]
      }
    }

    if (converge != 1) {
      cat(sprintf("Direction %d does not converge. Stop searching current joint block.\n", j))
      break
    }

    cat(sprintf("Direction %d converges.\n", j))
    Vi <- cbind(Vi, opt_v)
  }

  if(!is.null(Vi)){
    angles <- angles[, 1:ncol(Vi), drop=F]
  } else {
    angles <- angles[,1,drop=F]
  }


  return(list(Vi = Vi, angles = angles))
}
