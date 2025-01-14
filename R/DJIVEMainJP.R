DJIVEMainJP <- function(datablock, paramstruct = list(), truth = NULL) {

  # Initialize parameters
  nb <- length(datablock)
  dataname <- vector("list", nb)
  for (ib in seq_len(nb)) {
    dataname[[ib]] <- paste0("datablock", ib)
  }
  nsim <- 400
  theta0 <- 45
  optArgin <- list(0.5, 1000, 1.05, 50, 1e-3, 1e-3)
  iprint <- TRUE
  colCent <- FALSE
  rowCent <- FALSE
  figdir <- ""
  filterPerc <- 1 - (2 / (1 + sqrt(5))) # "Golden Ratio"
  noisepercentile <- rep(0.5, nb)

  if (!is.null(paramstruct)) {
    if (!is.null(paramstruct$dataname)) {
      dataname <- paramstruct$dataname
    }
    if (!is.null(paramstruct$nsim)) {
      nsim <- paramstruct$nsim
    }
    if (!is.null(paramstruct$theta0)) {
      theta0 <- paramstruct$theta0
    }
    if (!is.null(paramstruct$optArgin)) {
      optArgin <- paramstruct$optArgin
    }
    if (!is.null(paramstruct$iprint)) {
      iprint <- paramstruct$iprint
    }
    if (!is.null(paramstruct$figdir)) {
      figdir <- paramstruct$figdir
    }
    if (!is.null(paramstruct$colCent)) {
      colCent <- paramstruct$colCent
    }
    if (!is.null(paramstruct$rowCent)) {
      rowCent <- paramstruct$rowCent
    }
    if (!is.null(paramstruct$filterPerc)) {
      filterPerc <- paramstruct$filterPerc
    }
    if (!is.null(paramstruct$noisepercentile)) {
      noisepercentile <- paramstruct$noisepercentile
    }
  }

  rowSpaces <- vector("list", nb)
  datablockc <- vector("list", nb)
  for (ib in seq_len(nb)) {
    rowSpaces[[ib]] <- 0
    datablockc[[ib]] <- MatCenterJP(datablock[[ib]], colCent, rowCent)
  }

  if (!is.null(truth)) {
    for (ib in seq_len(nb)) {
      rowSpaces[[ib]] <- qr.Q(qr(t(truth[[ib]])))
    }
  }

  # Step 1: Estimate signal space and perturbation angle
  Phase1 <- DJIVESignalExtractJP(
    datablock = datablockc, nsim = nsim,
    iplot = FALSE, colCent = colCent, rowCent = rowCent, cull = filterPerc, noisepercentile = noisepercentile
  )
  VBars <- Phase1[[1]]
  UBars <- Phase1[[2]]
  phiBars <- Phase1[[3]]
  psiBars <- Phase1[[4]]
  rBars <- Phase1[[6]]
  VVHatCacheBars <- Phase1[[10]]
  UUHatCacheBars <- Phase1[[11]]

  # Step 2: Estimate joint (and partially joint) structure
  Phase2 <- DJIVEJointStrucEstimateJP(
    VBars, UBars, phiBars, psiBars, rBars, dataname, theta0, optArgin, iprint, figdir
  )

  outMap <- Phase2[[1]]
  keyIdxMap <- Phase2[[2]]
  jointBlockOrder <- Phase2[[4]]

  # Step 3: Reconstruct DJIVE decomposition
  outstruct <- DJIVEReconstructMJ(datablockc, dataname, outMap, keyIdxMap, jointBlockOrder, 0)

  outstruct$rBars <- rBars
  outstruct$phiBars <- phiBars
  outstruct$psiBars <- psiBars
  outstruct$VBars <- VBars
  outstruct$UBars <- UBars
  outstruct$VVHatCacheBars <- VVHatCacheBars
  outstruct$UUHatCacheBars <- UUHatCacheBars
  outstruct$jointBasisMapRaw <- outMap

  return(outstruct)
}

