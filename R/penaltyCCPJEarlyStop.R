penaltyCCPJPEarlyStop <- function(v0, Qo1, Qo2, Qc1, Qc2, Vorth, optArgin = list()) {
  # Inputs:
  #   v0 - Initial direction (vector)
  #   Qo1 - A matrix of the 1st quadratic form in the objective function
  #   Qo2 - A matrix of the 2nd quadratic form in the objective function
  #   Qc1 - A list of matrices for the 1st quadratic form in the constraints
  #   Qc2 - A list of matrices for the 2nd quadratic form in the constraints
  #   Vorth - A matrix of orthogonal constraints
  #   optArgin - A list of optimization tuning parameters: tau0, tau_max, mu, t_max, tol, delta
  #
  # Outputs:
  #   A list with the following elements:
  #   opt_v - The optimized direction vector
  #   cache_v - Cache of all vectors from each iteration
  #   cache_cvx_objval - Cache of objective values from each iteration
  #   cache_slack - Cache of slack variable values from each iteration
  #   converge - Flag indicating if the optimization converged

  optargs <- list(0.5, 1000, 1.05, 200, 1e-3, 1e-3)

  numvarargs <- length(optArgin)
  if (numvarargs > 0) {
    optargs[1:numvarargs] <- optArgin
  }

  tau0 <- optargs[[1]]
  tau_max <- optargs[[2]]
  mu <- optargs[[3]]
  t_max <- optargs[[4]]
  tol <- optargs[[5]]
  delta <- optargs[[6]]

  if (nrow(Vorth) == 0) {
    Vorth <- matrix(0, nrow(Qo1), 1)
  }

  cache_v <- vector("list", t_max)
  cache_cvx_objval <- rep(Inf, t_max)
  cache_slack <- vector("list", t_max)

  cache_v[[1]] <- v0
  nc <- length(Qc1)
  cache_slack[[1]] <- rep(Inf, nc + 2)

  converge <- 0  # flag to ensure two steps without change for convergence
  tau <- tau0

  updatePrintFormat <- paste0("%d ", paste(rep("%f ", length(Qc1) + 3), collapse = ""), "\n")

  for (t in 2:t_max) {
    if (t %% 10 == 0) {
      cat(sprintf("Iteration %d\n", t))
    }


    result_penalty <- ccpSubOptJPEarlyStop(cache_v[[t - 1]], Qo1, Qo2, Qc1, Qc2, Vorth, tau)

    ########
    # Before accessing the result$getValue(v), add a safety check
    if (!is.null(result_penalty$v) && length(result_penalty$v) > 0) {
      cache_v[[t]] <- result_penalty$v
    } else {
      cat(sprintf("Iteration %d - Solver returned an invalid value. Skipping.\n", t))
      break
    }
    #######
    #cache_v[[t]] <- result_penalty$v
    cache_slack[[t]] <- result_penalty$slack
    cache_cvx_objval[t] <- result_penalty$objval

    if (any(is.na(cache_v[[t]]))) {
      cat(sprintf("NaN solution appears in iteration %d\n", t))
      t <- t - 1
      break
    }

    cache_v[[t]] <- cache_v[[t]] / norm(cache_v[[t]], "2")

    curr_objval <- cache_cvx_objval[t] + tau * sum(cache_slack[[t]])
    if (is.infinite(cache_cvx_objval[t - 1])) {
      pre_objval <- 1e6  #  Inf
    } else {
      pre_objval <- cache_cvx_objval[t - 1] + (tau / mu) * sum(cache_slack[[t - 1]])
    }
    curr_slack <- cache_slack[[t]][1:(length(Qc1) + 2)]
    cat(do.call(sprintf, c(list(updatePrintFormat, t), as.list(curr_slack), list(curr_objval - pre_objval))))

    # Check convergence based on objective value change and slack variables
    if (norm(curr_objval - pre_objval, type = "2") < tol && sum(cache_slack[[t]]) <= delta) {
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

  cache_v <- cache_v[1:t]
  cache_cvx_objval <- cache_cvx_objval[1:t]
  cache_slack <- cache_slack[1:t]

  opt_v <- cache_v[[t]]
  return(list(opt_v = opt_v, cache_v = cache_v, cache_cvx_objval = cache_cvx_objval, cache_slack = cache_slack, converge = converge))
}

