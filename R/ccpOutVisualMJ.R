ccpOutVisualMJ <- function(angleHats, phiBars, dataname, iprint = NULL, figdir = NULL, figname = NULL) {
  # Set default figure name if not provided
  if (is.null(figname) || figname == "") {
    figname <- "opt_progress"
  }

  nb <- length(phiBars)
  T <- length(angleHats[[1]])
  idx <- 1:T


  #par(mfrow = c(1, nb), oma = c(0, 0, 2, 0), mar = c(1, 1, 1, 1))
  par(mfrow = c(1, nb), oma = c(0, 0, 1, 0), mar = c(2, 2, 1, 1))
  for (ib in 1:nb) {
    # for each data block
    plot(idx, angleHats[[ib]], type = "l", lwd = 2,
         xlab = "Index", ylab = "Projected Angle",
         main = paste0(dataname[[ib]], "\n", figname),
         xlim = c(1, T), ylim = c(0, max(angleHats[[ib]], phiBars[ib], na.rm = TRUE)))

    # Add a horizontal line for phiBars
    abline(h = phiBars[ib], col = "green", lty = 2, lwd = 2)

    legend("topright", legend = c(paste0("Estimated Angle ", ib), paste0("Perturbation Angle ", ib)),
           col = c("black", "green"), lty = c(1, 2), lwd = c(2, 2), bty = "n")
  }


  if (!is.null(iprint) && iprint == 1) {
    # If figdir is not provided or doesn't exist, set it to the current working directory
    if (is.null(figdir) || !dir.exists(figdir)) {
      message("No valid figure directory found! Saving to the current folder.")
      figdir <- getwd()  # Use the current working directory as a fallback
    }

    savestr <- file.path(figdir, paste0(figname, ".png"))

    tryCatch({
      png(savestr, width = 1500, height = 500)

      par(mfrow = c(1, nb), oma = c(0, 0, 2, 0), mar = c(1, 1, 1, 1))
      for (ib in 1:nb) {
        plot(idx, angleHats[[ib]], type = "l", lwd = 2,
             xlab = "Index", ylab = "Projected Angle",
             main = paste0(dataname[[ib]], "\n", figname),
             xlim = c(1, T), ylim = c(0, max(angleHats[[ib]], phiBars[ib], na.rm = TRUE)))
        abline(h = phiBars[ib], col = "green", lty = 2, lwd = 2)
      }
      dev.off()
      message("Figure saved successfully in: ", savestr)
    }, error = function(e) {
      message("Failed to save figure! Error: ", e)
    })
  }
}
