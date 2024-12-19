#' Convert Logical Index Array to Number (Binary to Decimal)
#'
#' This function takes a logical vector and converts it to a decimal number,
#' interpreting the logical values as a binary representation.
#'
#' @param blockIn A logical vector indicating which block index is in (TRUE/FALSE).
#' @return A numeric value representing the decimal conversion of the binary input.
Idx2numMJ <- function(blockIn) {
  # Ensure input is a logical vector
  if (!is.logical(blockIn)) {
    stop("Input must be a logical vector.")
  }

  nb <- length(blockIn)
  t <- 0
  for (i in 1:nb) {
    if (blockIn[i]) {
      t <- t + 2^(i - 1)
    }
  }
  return(t)
}

ccpSubOptJPLoadInfo <- function(v0, Qo1, Qo2, Qc1, Qc2, Qc1Load, Qc2Load, Vo, tau) {
  # Inputs:
  #   v0 - Initial direction
  #   Qo1 - a matrix of the 1st quadratic form in the objective function
  #   Qo2 - a matrix of the 2nd quadratic form in the objective function
  #   Qc1, Qc2 - list of matrices for quadratic constraints
  #   Qc1Load, Qc2Load - list of matrices for quadratic constraints on load
  #   Vo - a matrix of orthogonal constraints
  #   tau - multiplier of slack variables

  nc <- length(Qc1)
  ncload <- length(Qc1Load)
  n <- length(v0)
  loadslackscale <- rep(1, ncload)
  for (ipre in seq_len(ncload)) {
    if (sum(Qc1Load[[ipre]] != 0) > 0) {
      #loadslackscale[ipre] <- svd(Qc1Load[[ipre]])$d[1]
      loadslackscale[ipre] <- RSpectra::svds(Qc1Load[[ipre]], k = 1, nu = 0, nv = 0)$d
    }
  }
  #print(loadslackscale)
  ro <- ncol(Vo)

  v <- CVXR::Variable(n)
  slack <- CVXR::Variable(nc + ncload + 2)
  objective <- CVXR::quad_form(v, Qo1) - 2 * t(v0) %*% Qo2 %*% v + CVXR::quad_form(v0, Qo2) + tau * sum(slack)

  constraints <- list()
  # Add constraints for Qc1 and Qc2
  for (ic in 1:nc) {
    constraints <- c(constraints,
                     CVXR::quad_form(v, Qc1[[ic]]) - 2 * t(v0) %*% Qc2[[ic]] %*% v + CVXR::quad_form(v0, Qc2[[ic]]) <= slack[ic])
  }

  for (ic in 1:ncload) {
    constraints <- c(constraints,
                     CVXR::quad_form(v, Qc1Load[[ic]]) - 2 * t(v0) %*% Qc2Load[[ic]] %*% v + CVXR::quad_form(v0, Qc2Load[[ic]]) <= (slack[nc + ic]/loadslackscale[ic]))
  }

  constraints <- c(constraints,
                   CVXR::sum_squares(v) - 1 <= slack[nc + ncload + 1],
                   1 - 2 * t(v0) %*% v + CVXR::sum_squares(v) <= slack[nc + ncload + 2],
                   slack >= 0,
                   t(Vo) %*% v == matrix(0, nrow = ro, ncol = 1))

  problem <- CVXR::Problem(CVXR::Minimize(objective), constraints)
  result <- CVXR::solve(problem)  # solver = "ECOS"
  #result <- solve(problem, solver = "SCS")
  #result <- solve(problem, solver = "ECOS", reltol = 1e-8, abstol = 1e-8)

  opt_v <- result$getValue(v)
  opt_slack <- result$getValue(slack)
  cvx_objval <- t(opt_v) %*% Qo1 %*% opt_v - 2 * t(v0) %*% Qo2 %*% opt_v + t(v0) %*% Qo2 %*% v0
  cvx_optval <- result$value
  return(list(v = opt_v, slack = opt_slack, objval = cvx_optval, status = result$status))
}



# v0 <- matrix(c(1, 2, 3), nrow = 3, ncol = 1)
# Qo1 <- diag(3)
# Qo2 <- 0.5 * diag(3)
# Qc1 <- list(diag(3), 2 * diag(3))
# Qc2 <- list(diag(3), 0.5 * diag(3))
# Qc1Load <- list(3 * diag(3), 3 * diag(3))
# Qc2Load <- list(2 * diag(3), diag(3))
# Vo <- matrix(c(1, 0, 0, 0, 1, 0), nrow = 3, byrow = F)
# tau <- 1
#
# result <- ccpSubOptJPLoadInfo(v0, Qo1, Qo2, Qc1, Qc2, Qc1Load, Qc2Load, Vo, tau)
#
# # Print the results
# cat("Optimized v:\n")
# print(result$v)
# cat("Slack variables:\n")
# print(result$slack)
# cat("Objective value:\n")
# print(result$objval)
# cat("Status:\n")
# print(result$status)



takeNormOfEachColumnJP <- function(inmat) {

  d <- nrow(inmat)
  sumsqs <- sqrt(colSums(inmat^2))
  outmat <- inmat / matrix(rep(sumsqs, each = d), nrow = d)

  return(outmat)
}


# test_matrix <- matrix(c(3, 4, 0, 5, 1, 2), nrow = 3, ncol = 2)
#
# normalized_matrix <- takeNormOfEachColumnJP(test_matrix)
# print(normalized_matrix)



penaltyCCPJPEarlyStopLoadInfo <- function(v0, Qo1, Qo2, Qc1, Qc2, Qc1Load, Qc2Load, Vorth, optArgin = list(0.5, 1000, 1.05, 200, 1e-3, 1e-3)) {
  # Inputs:
  #   v0 - Initial direction (vector)
  #   Qo1, Qo2 - Matrices for the quadratic form in the objective function
  #   Qc1, Qc2 - Lists of matrices for the quadratic form in the constraints
  #   Qc1Load, Qc2Load - Lists of matrices for the load constraints
  #   Vorth - Matrix of orthogonal constraints
  #   optArgin - List of optimization tuning parameters
  #
  # Outputs:
  #   A list with optimized direction vector, cached vectors, objective values, slack variables, and convergence flag

  # Set default optimization parameters
  if(is.null(optArgin)){
    optargs <- list(0.5, 1000, 1.05, 200, 1e-3, 1e-3)
  } else {
    numvarargs <- length(optArgin)
    if(numvarargs != 6) stop("Have to specify 6 tuning parameters, ie optArgin has to be a vector of length 6.")
    optargs[1:numvarargs] <- optArgin
  }
  tau0 <- optargs[[1]]
  tau_max <- optargs[[2]]
  mu <- optargs[[3]]
  t_max <- optargs[[4]]
  tol <- optargs[[5]]
  delta <- optargs[[6]]

  # Orthogonalize initial direction if needed
  if (nrow(Vorth) == 0) {
    Vorth <- matrix(0, nrow(Qo1), 1)
  } else {
    orthBasis <- qr.Q(qr(Vorth))
    #print(Vorth)
    #print(orthBasis)
    v0 <- takeNormOfEachColumnJP((diag(nrow(Vorth)) - orthBasis %*% t(orthBasis)) %*% v0)
  }

  #print(v0)

  # Initialize caches
  cache_v <- vector("list", t_max)
  cache_cvx_objval <- rep(Inf, t_max)
  cache_slack <- vector("list", t_max)

  cache_v[[1]] <- v0
  nc <- length(Qc1)
  cache_slack[[1]] <- rep(Inf, nc + length(Qc1Load) + 2)

  converge <- 0 # flag to ensure two steps without change for convergence
  tau <- tau0

  updatePrintFormat <- paste0("%d ", paste(rep("%f ", length(Qc1) + length(Qc1Load) + 3), collapse = ""), "\n")

  for (t in 2:t_max) {
    if (t %% 10 == 0) {
      cat(sprintf("Iteration %d\n", t))
    }

    result <- ccpSubOptJPLoadInfo(cache_v[[t - 1]], Qo1, Qo2, Qc1, Qc2, Qc1Load, Qc2Load, Vorth, tau)
    cache_v[[t]] <- result[[1]]
    cache_slack[[t]] <- result[[2]]
    cache_cvx_objval[t] <- result[[3]]
    cvx_status <- result[[4]]

    if (cvx_status != "optimal") {            #Solved
      cat(sprintf("Iteration %d Status: %s\n", t, cvx_status))
      loadslackscale <- rep(1, length(Qc1Load))
      for (ipre in seq_along(Qc1Load)) {
        if (sum(Qc1Load[[ipre]] != 0)) {
          loadslackscale[ipre] <- svd(Qc1Load[[ipre]])$d[1]
        }
      }
      print(loadslackscale)


      if (!grepl("inaccurate", cvx_status)) { #Inaccurate/Solved, Inaccurate/Unbounded, Inaccurate/Infeasible
        cat(sprintf(updatePrintFormat, t, cache_slack[[t]], cache_cvx_objval[t]))
        t <- t - 1
        converge <- 2
        break
      }
    }

    # Normalize the current solution
    cache_v[[t]] <- cache_v[[t]] / norm(cache_v[[t]], type = "2")
    # Compute the objective value difference
    curr_objval <- cache_cvx_objval[t] + tau * sum(cache_slack[[t]])
    if (is.infinite(cache_cvx_objval[t - 1])) {
      pre_objval <- 1e6  # Assign a large finite value if the first value is Inf
    } else {
      pre_objval <- cache_cvx_objval[t - 1] + (tau / mu) * sum(cache_slack[[t - 1]])
    }

    # Print the iteration details
    cat(do.call(sprintf, c(list(updatePrintFormat, t), as.list(cache_slack[[t]]), list(curr_objval - pre_objval))))

    # Check convergence based on objective value change and slack variables
    curr_slack <- cache_slack[[t]]
    if ((abs(curr_objval - pre_objval) < tol && sum(cache_slack[[t]]) <= delta) ||
        (sum(curr_slack[1:(length(Qc1) + length(Qc1Load))]) <= delta^2)) {
      if (converge == 1) {
        break
      } else {
        converge <- 1
      }
    } else {
      converge <- 0
    }

    tau <- min(mu * tau, tau_max)
  }

  # Trim cache to actual iterations
  cache_v <- cache_v[1:t]
  cache_cvx_objval <- cache_cvx_objval[1:t]
  cache_slack <- cache_slack[1:t]
  opt_v <- cache_v[[length(cache_v)]]

  # Return the results
  output <- list(opt_v, cache_v, cache_cvx_objval, cache_slack, converge)
  return(output)
}




# v0 <- matrix(c(1, 2, 3, 4), nrow = 4, ncol = 1)
# Qo1 <- diag(4)
# Qo2 <- 0.5 * diag(4)
# Qc1 <- list(diag(4), 2 * diag(4))
# Qc2 <- list(diag(4), 0.5 * diag(4))
# Qc1Load <- list(3 * diag(4), 4 * diag(4))
# Qc2Load <- list(2 * diag(4), diag(4))
# Vorth <- matrix(c(1, 0, 0, 1, 0, 1, 0, 0), nrow = 4, byrow = TRUE)  # Fixed orthogonal matrix
#
# result <- penaltyCCPJPEarlyStopLoadInfo(v0, Qo1, Qo2, Qc1, Qc2, Qc1Load, Qc2Load, Vorth)

# Print the results
# cat("Optimized v:\n")
# print(result[[1]])
# cat("Slack variables:\n")
# print(result[[4]])
# cat("Objective values:\n")
# print(result[[3]])
# cat("Converged after iterations:\n")
# print(length(result[[2]]))



#' Calculate the Projection Angle between a Vector and a Basis Matrix
#'
#' The function `projAngleMJ` calculates the angle between a given vector and its projection
#' onto the subspace spanned by a basis matrix. The resulting angle provides a measure of
#' how closely the vector aligns with the subspace.
#'
#' @param x A numeric vector of length \code{n}, representing an \eqn{n \times 1} vector.
#' @param V A numeric matrix of dimension \eqn{n \times r}, representing a basis matrix of the subspace.
#'
#' @return A numeric value representing the projected angle of \code{x} on \code{V}, in degrees.
#' The value is in the range \eqn{[0, 180]}.
#'
#' @export
projAngleMJ <- function(x, V) {
  # Ensure inputs are compatible
  if (!is.numeric(x) || !is.matrix(V)) {
    stop("Input x must be a numeric vector, and V must be a numeric matrix.")
  }

  if (length(x) != nrow(V)) {
    stop("The length of vector x must match the number of rows in matrix V.")
  }

  # Project x onto the subspace spanned by V
  px <- V %*% t(V) %*% x

  # Calculate the angle between the vector x and its projection px
  norm_px <- norm(px, type = "2")
  norm_x <- norm(x, type = "2")

  if (norm_px == 0 || norm_x == 0) {
    warning("The norm of the projection or input vector is zero. Returning angle of 90 degrees.")
    return(90)
  }

  cos_theta <- sum(px * x) / (norm_px * norm_x)
  # Clip cos_theta to avoid numerical issues with acos
  cos_theta <- max(min(cos_theta, 1), -1)

  angle <- acos(cos_theta) * (180 / pi)

  return(angle)
}



#' Analyze Optimization Output by Calculating Projection Angles
#'
#' The function `ccpOutAnalysisMJ` computes the projected angles between each cached optimization vector
#' and each adjusted signal row space. The projection angle serves as a metric for evaluating alignment
#' between the vectors obtained during optimization and the reference row spaces.
#'
#' @param cache_v A list of optimization vectors obtained from iterations of the optimization algorithm.
#' @param VBars A list of adjusted signal row spaces (matrices).
#'
#' @return A list of length equal to the number of signal row spaces. Each element is a numeric vector of
#' angles for each iteration of the optimization process.
#'
#' @export
ccpOutAnalysisMJ <- function(cache_v, VBars) {
  # Number of adjusted signal row spaces
  nb <- length(VBars)
  # Number of cached optimization vectors
  T <- length(cache_v)
  # Initialize angleHats as a list to store results
  angleHats <- vector("list", nb)

  # Loop over each adjusted signal row space
  for (ib in 1:nb) {
    # Initialize angles vector with -1 for each iteration (default initialization)
    angles <- rep(-1, T)

    # Loop over each cached optimization vector to compute the projection angle
    for (t in 1:T) {
      # Calculate the projected angle of each vector on the row space
      angles[t] <- projAngleMJ(cache_v[[t]], VBars[[ib]])
    }

    # Store the calculated angles in the angleHats list
    angleHats[[ib]] <- angles
  }

  return(angleHats)
}


#' Visualize the Optimization Progress of Projected Angles
#'
#' The `ccpOutVisualMJ` function visualizes the projected angles for each data block
#' over optimization iterations. The function allows for saving the plot as a PNG file
#' to the specified directory.
#'
#' @param angleHats A list of numeric vectors representing the projected angles of each data block at different iterations.
#' @param phiBars A numeric vector representing the perturbation angles for each data matrix.
#' @param dataname A character vector containing the names of the data blocks.
#' @param iprint An integer value (0 or 1). If \code{1}, the plot will be saved to the disk. Default is \code{NULL}.
#' @param figdir A character string specifying the directory to save the figure. Default is \code{NULL}, which saves to the current directory if \code{iprint = 1}.
#' @param figname A character string specifying the name of the figure file. Default is \code{"opt_progress"}.
#'
#' @return None. The function generates and optionally saves a plot.
#'
#' @importFrom grDevices png dev.off
#' @importFrom graphics par plot abline legend
#' @export
ccpOutVisualMJ <- function(angleHats, phiBars, dataname, iprint = FALSE, figdir = NULL, figname = NULL) {
  if (is.null(figname) || figname == "") {
    figname <- "opt_progress"
  }

  nb <- length(phiBars)
  T <- length(angleHats[[1]])
  idx <- 1:T

  par(mfrow = c(1, nb), oma = c(0, 0, 2, 0), mar = c(3, 3, 2, 1))

  for (ib in 1:nb) {
    plot(idx, angleHats[[ib]], type = "l", lwd = 2,
         xlab = "Iteration Index", ylab = "Projected Angle",
         main = paste0(dataname[[ib]], "\n", figname),
         xlim = c(1, T), ylim = c(0, max(angleHats[[ib]], phiBars[ib], na.rm = TRUE)))

    abline(h = phiBars[ib], col = "green", lty = 2, lwd = 2)
    legend("topright", legend = c("Estimated Angle", "Perturbation Angle"),
           col = c("black", "green"), lty = c(1, 2), lwd = c(2, 2), bty = "n")
  }

  if (iprint) {
    if (is.null(figdir) || figdir == "") {
      figdir <- paste0(getwd(), "/DIVAS_figures")

      message("No figure directory Given. Will save to ", figdir)

      # figdir <- path.expand("~/DIVAS_figures")  # default
    }

    figdir <- gsub("//", "/", figdir)

    if (!dir.exists(figdir)) {
      figdir <- paste0(getwd(), "/DIVAS_figures")
      message("Figure directory does not exist. Creating directory: ", figdir)
      dir.create(figdir, recursive = TRUE)
    }

    if (file.access(figdir, mode = 2) != 0) {
      stop("No write permission for the specified directory: ", figdir)
    }

    savestr <- file.path(figdir, paste0(figname, ".png"))
    #cat("Figure directory:", figdir, "\n")
    #cat("Save path:", savestr, "\n")


    tryCatch({
      png(savestr, width = 1500, height = 500)
      par(mfrow = c(1, nb), oma = c(0, 0, 2, 0), mar = c(3, 3, 2, 1))
      for (ib in 1:nb) {
        plot(idx, angleHats[[ib]], type = "l", lwd = 2,
             xlab = "Iteration Index", ylab = "Projected Angle",
             main = paste0(dataname[[ib]], "\n", figname),
             xlim = c(1, T), ylim = c(0, max(angleHats[[ib]], phiBars[ib], na.rm = TRUE)))
        abline(h = phiBars[ib], col = "green", lty = 2, lwd = 2)
      }

      dev.off()
      message("Figure saved successfully in: ", savestr)
    }, error = function(e) {
      message("Failed to save figure! Error: ", e$message)
    })
  }
}



# ccpOutVisualMJ <- function(angleHats, phiBars, dataname, iprint = NULL, figdir = NULL, figname = NULL) {
#   if (is.null(figname) || figname == "") {
#     figname <- "opt_progress"
#   }
#
#   nb <- length(phiBars)
#   # Number of iterations or time points in angleHats
#   T <- length(angleHats[[1]])
#   idx <- 1:T
#
#   # Create a plot window with specified layout for the plots
#   par(mfrow = c(1, nb), oma = c(0, 0, 2, 0), mar = c(3, 3, 2, 1))
#
#   # Loop through each data block and create a plot for the corresponding angles
#   for (ib in 1:nb) {
#     plot(idx, angleHats[[ib]], type = "l", lwd = 2,
#          xlab = "Iteration Index", ylab = "Projected Angle",
#          main = paste0(dataname[[ib]], "\n", figname),
#          xlim = c(1, T), ylim = c(0, max(angleHats[[ib]], phiBars[ib], na.rm = TRUE)))
#
#     # Add a horizontal line indicating the perturbation angle
#     abline(h = phiBars[ib], col = "green", lty = 2, lwd = 2)
#
#     # Add a legend to indicate the lines plotted
#     legend("topright", legend = c("Estimated Angle", "Perturbation Angle"),
#            col = c("black", "green"), lty = c(1, 2), lwd = c(2, 2), bty = "n")
#   }
#
#
#
#
#
#   if (!is.null(iprint) && iprint == 1) {
#     # If figdir is not provided or doesn't exist, set it to the current working directory
#     # if (is.null(figdir) || !dir.exists(figdir) || figdir == "") {
#     #   message("No valid figure directory found! Saving to the current folder.")
#     #   figdir <- getwd()  # Use the current working directory as a fallback
#     # }
#
#     if (is.null(figdir) || figdir == "") {
#       warning("Figure directory is empty. Using current working directory.")
#       figdir <- getwd()
#     }
#     if (!dir.exists(figdir)) {
#       warning("Figure directory does not exist. Creating directory.")
#       dir.create(figdir, recursive = TRUE)
#     }
#
#     cat("Figure directory:", figdir, "\n")
#
#
#     figdir <- gsub("/+$", "", figdir)
#     savestr <- file.path(figdir, paste0(figname, ".png"))
#
#
#     print(paste("Figure directory:", figdir))
#     print(paste("Save path:", savestr))
#
#     # 检查写权限
#     if (file.access(figdir, mode = 2) != 0) {
#       stop("No write permission for the specified directory: ", figdir)
#     }
#
#     tryCatch({
#       png(savestr, width = 1500, height = 500)
#       par(mfrow = c(1, nb), oma = c(0, 0, 2, 0), mar = c(3, 3, 2, 1))
#       for (ib in 1:nb) {
#         plot(idx, angleHats[[ib]], type = "l", lwd = 2,
#              xlab = "Iteration Index", ylab = "Projected Angle",
#              main = paste0(dataname[[ib]], "\n", figname),
#              xlim = c(1, T), ylim = c(0, max(angleHats[[ib]], phiBars[ib], na.rm = TRUE)))
#         abline(h = phiBars[ib], col = "green", lty = 2, lwd = 2)
#       }
#
#       dev.off()
#       message("Figure saved successfully in: ", savestr)
#     }, error = function(e) {
#       message("Failed to save figure! Error: ", e)
#     })
#   }
# }
