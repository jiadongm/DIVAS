suppressWarnings(library(CVXR, warn.conflicts=FALSE))

ccpSubOptJPEarlyStop <- function(v0, Qo1, Qo2, Qc1, Qc2, Vo, tau) {
  # Inputs:
  #   v0 - Initial direction
  #   Qo1 - a matrix of the 1st quadratic form in the objective function
  #   Qo2 - a matrix of the 2nd quadratic form in the objective function
  #   Qc1 - a list of matrices of the 1st quadratic form in the constraints
  #   Qc2 - a list of matrices of the 2nd quadratic form in the constraints
  #   Vo  - a matrix of orthogonal constraints
  #   tau - multiplier of slack variables
  #
  # Outputs:
  #   result - list containing optimized direction, slack variable values, and objective value

  nc <- length(Qc1)
  n <- length(v0)
  ro <- ncol(Vo)

  # cvx_begin quiet;
  v <- CVXR::Variable(n)
  slack <- CVXR::Variable(nc + 2)

  objective <- CVXR::quad_form(v, Qo1) - 2 * t(v0) %*% Qo2 %*% v + CVXR::quad_form(v0, Qo2) + tau * sum(slack)
  constraints <- list()

  # constraints for each Qc1 and Qc2
  for (ic in 1:nc) {
    constraints <- c(constraints,
                     CVXR::quad_form(v, Qc1[[ic]]) - 2 * t(v0) %*% Qc2[[ic]] %*% v + CVXR::quad_form(v0, Qc2[[ic]]) <= slack[ic])
  }
  # Addition
  constraints <- c(constraints,
                   CVXR::sum_squares(v) - 1 <= slack[nc + 1],
                   1 - 2 * t(v0) %*% v + CVXR::sum_squares(v) <= slack[nc + 2],
                   slack >= 0,
                   t(Vo) %*% v == matrix(0, nrow = ro, ncol = 1))
  # cvx_end

  problem <- CVXR::Problem(Minimize(objective), constraints)
  result <- CVXR::solve(problem)
  #result <- solve(problem, solver = "ECOS")
  v_opt <- result$getValue(v)
  slack_opt <- result$getValue(slack)
  cvx_objval <- t(v_opt) %*% Qo1 %*% v_opt - 2 * t(v0) %*% Qo2 %*% v_opt + t(v0) %*% Qo2 %*% v0

  return(list(v = v_opt, slack = slack_opt, objval = cvx_objval,  status=result$status))
}
