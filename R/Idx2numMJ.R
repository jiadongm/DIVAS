# Convert Logical Index Array to Number (Binary to Decimal)
# 
# This function takes a logical vector and converts it to a decimal number,
# interpreting the logical values as a binary representation.
# 
# Inputs:
#   blockIn - A logical vector indicating which block index is in (TRUE/FALSE).
# 
# Outputs:
#   t - A numeric value representing the decimal conversion of the binary input.

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


# # Test cases
# cat("Test 1: Simple binary conversion\n")
# blockIn <- c(TRUE, FALSE, TRUE)  # Binary 101
# expected <- 5
# result <- Idx2numMJ(blockIn)
# cat("Input:", blockIn, "\nExpected:", expected, "\nResult:", result, "\n\n")
# 
# cat("Test 2: All FALSE should return 0\n")
# blockIn <- c(FALSE, FALSE, FALSE)
# expected <- 0
# result <- Idx2numMJ(blockIn)
# cat("Input:", blockIn, "\nExpected:", expected, "\nResult:", result, "\n\n")
# 
# cat("Test 3: All TRUE should return maximum value\n")
# blockIn <- c(TRUE, TRUE, TRUE)  # Binary 111
# expected <- 7
# result <- Idx2numMJ(blockIn)
# cat("Input:", blockIn, "\nExpected:", expected, "\nResult:", result, "\n\n")
# 
# cat("Test 4: Mixed TRUE and FALSE\n")
# blockIn <- c(FALSE, TRUE, TRUE)  # Binary 110
# expected <- 6
# result <- Idx2numMJ(blockIn)
# cat("Input:", blockIn, "\nExpected:", expected, "\nResult:", result, "\n")
