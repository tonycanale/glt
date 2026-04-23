identify_lambda <- function(Lambda, p, nsample = NULL, 
                            return_delta = TRUE,
                            return_pivot_info = FALSE) {
  if (is.null(nsample)) {
    nsample <- length(fit$lambda)
  }
  
  # Initialize lists
  lambda_ordered <- list()
  delta_ordered <- list()
  pivot_location <- list()
  c <- 0
  
  # Step 1: Filter and order samples based on 3579 rule
  for (i in 1:nsample) {
    delta <- matrix(fit$delta[[i]][, c(fit$activeFactors[[i]] == TRUE)],
                    nrow = p)
    lambda <- matrix(fit$lambda[[i]][, c(fit$activeFactors[[i]] == TRUE)],
                     nrow = p)
    
    # Check if zero columns are present (not identifiable)
    if (sum(colSums(delta) == 0) == 0) {
      
      # If zero rows are present, check counting rule only on non-zero rows
      delta_tocheck <- as.matrix(delta[rowSums(delta) > 0, ])
      
      # Check 3579 rule
      if (counting_rule_holds(delta_tocheck)) {
        c <- c + 1
        
        # Order pivots decreasing and sign switch when needed
        neword <- reorder_and_sign_swap(delta, lambda)
        lambda_ordered[[c]] <- neword$lambda
        delta_ordered[[c]] <- neword$delta
        pivot_location[[c]] <- neword$pivot_location
      }
    }
  }
  
  # Step 2: Find modal pivot configuration
  modal_config <- find_modal_vector(pivot_location)
  modal_rank <- length(modal_config$mode)
  
  # Step 3: Keep only lambdas with modal configuration
  c <- 0
  lambda_identified <- list()
  delta_identified <- list()
  
  for (i in 1:length(lambda_ordered)) {
    if (length(pivot_location[[i]]) == modal_rank) {
      if (all(pivot_location[[i]] == modal_config$mode)) {
        c <- c + 1
        lambda_identified[[c]] <- lambda_ordered[[i]]
        delta_identified[[c]] <- delta_ordered[[i]]
      }
    }
  }
  
  # Prepare output
  result <- list(
    lambda = lambda_identified,
    modal_config = modal_config$mode,
    modal_rank = modal_rank,
    n_identified = length(lambda_identified)
  )
  
  if (return_delta) {
    result$delta <- delta_identified
  }
  
  if (return_pivot_info) {
    result$all_pivot_locations <- pivot_location
    result$lambda_ordered <- lambda_ordered
    result$delta_ordered <- delta_ordered
  }
  
  return(result)
}


# function to reorder and sign swap based on delta and lambda
reorder_and_sign_swap <- function(delta, lambda) {
  # Ensure both have same dimensions
  if (!all(dim(delta) == dim(lambda))) {
    stop("delta and lambda must have the same dimensions")
  }
  
  p <- nrow(lambda)
  k <- ncol(lambda)
  
  # Identify first nonzero position and sign for each column
  pivot_pos <- rep(NA_integer_, k)
  pivot_sign <- rep(1, k)
  
  for (h in seq_len(k)) {
    nz <- which(lambda[, h] != 0)
    if (length(nz) > 0) {
      pivot_pos[h] <- nz[1]
      pivot_sign[h] <- sign(lambda[nz[1], h])
    } else {
      pivot_pos[h] <- p + 1  # if column is all zeros, put it last
    }
  }
  
  # Order columns by first nonzero position
  order_idx <- order(pivot_pos)
  pivot_pos_sorted = sort(pivot_pos)
  
  # Apply the order
  delta_new <- delta[, order_idx, drop = FALSE]
  lambda_new <- lambda[, order_idx, drop = FALSE]
  pivot_sign <- pivot_sign[order_idx]
  
  # Flip sign where needed
  for (h in seq_len(k)) {
    if (pivot_sign[h] < 0) {
      lambda_new[, h] <- -lambda_new[, h]
    }
  }
  
  list(delta = delta_new, lambda = lambda_new, order = order_idx,
       pivot_location = pivot_pos_sorted)
}

# find modal vector
find_modal_vector <- function(pivot_list) {
  # convert each vector to a string key
  keys <- vapply(pivot_list, function(v) paste(v, collapse = ","), character(1))
  
  # count frequency
  freq <- table(keys)
  
  # find the most frequent key
  modal_key <- names(freq)[which.max(freq)]
  
  # convert back to numeric vector
  if (nzchar(modal_key)) {
    modal_vec <- as.integer(strsplit(modal_key, ",", fixed = TRUE)[[1]])
  } else {
    modal_vec <- integer(0)  # handle empty vector case
  }
  
  list(mode = modal_vec, count = max(freq), freq_table = freq)
}