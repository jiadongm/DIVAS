projAngleMJ <- function(x, V) {
  # Project x onto the subspace spanned by V
  px <- V %*% t(V) %*% x

  # Calculate the angle between the vector x and its projection px
  angle <- acos(sum(px * x) / (norm(px, type = "2") * norm(x, type = "2"))) * 180 / pi

  return(angle)
}



ccpOutAnalysisMJ <- function(cache_v, VBars) {
  # Number of adjusted signal row spaces
  nb <- length(VBars)
  # Number of cached optimization vectors
  T <- length(cache_v)
  angleHats <- vector("list", nb)


  for (ib in 1:nb) {
    # Initialize angles with -1 for each iteration
    angles <- rep(-1, T)
    # Loop over each cached optimization vector
    for (t in 1:T) {
      angles[t] <- projAngleMJ(cache_v[[t]], VBars[[ib]])
    }
    angleHats[[ib]] <- angles
  }

  return(angleHats)
}
