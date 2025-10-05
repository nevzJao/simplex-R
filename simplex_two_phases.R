#' Solves a linear programming problem using the Two-Phase Simplex method.
#' This function is designed to maximize an objective function.
#'
#' @param A A numeric matrix of the constraint coefficients. Each row represents a constraint,
#' and each column represents a decision variable.
#' @param obj A numeric vector of the objective function coefficients.
#' @param rhs A numeric vector of the right-hand side values of the constraints.
#' @param relactions A character vector of the relational operators for each constraint.
#' Each element must be one of "<=", ">=", or "=".
#'
#' @return A list containing the optimal solution vector (`solution`) and the
#' optimal value of the objective function (`optimal_value`). Returns a string
#' message for infeasible or unbounded problems.

simplex <- function(A, obj, rhs, relactions) {
  
  # --- INITIALIZATION AND TABLEAU SETUP ---
  
  num_decision_vars <- ncol(A)
  num_constraints <- nrow(A)
  
  # Calculate the number of slack/excess and artificial variables needed
  num_slack <- sum(relactions == "<=") # This variable is informational
  num_artificial <- sum(relactions %in% c(">=", "="))
  
  # Define column names for all variables (decision, slack/surplus, artificial)
  var_names <- c(
    paste0("x", 1:num_decision_vars),
    if (any(relactions %in% c("<=", ">="))) paste0("s", 1:(num_constraints)), # 's' can be slack or surplus
    if (num_artificial > 0) paste0("a", 1:num_artificial)
  )
  
  total_cols <- num_decision_vars + num_constraints + num_artificial
  
  # Create the initial tableau matrix
  tableau <- matrix(0, nrow = num_constraints, ncol = total_cols)
  tableau[, 1:num_decision_vars] <- A
  
  # Vector to store the basic variable for each row
  base <- character(num_constraints)
  
  # Counters for slack/surplus and artificial variable columns
  col_s_count <- 0
  col_a_count <- 0
  original_artificial_indices <- c() # Stores the indices of artificial columns for later removal
  
  # Populate the tableau with slack, surplus, and artificial variables
  for (i in 1:num_constraints) {
    col_s_count <- col_s_count + 1
    
    if (relactions[i] == "<=") {
      # Add a slack variable
      tableau[i, num_decision_vars + col_s_count] <- 1
      base[i] <- paste0("s", col_s_count)
      
    } else if (relactions[i] == ">=") {
      # Add a surplus variable (-1) and an artificial variable (+1)
      col_a_count <- col_a_count + 1
      tableau[i, num_decision_vars + col_s_count] <- -1
      tableau[i, num_decision_vars + num_constraints + col_a_count] <- 1
      base[i] <- paste0("a", col_a_count)
      original_artificial_indices <- c(original_artificial_indices, num_decision_vars + num_constraints + col_a_count)
      
    } else if (relactions[i] == "=") {
      # Add an artificial variable (+1)
      col_a_count <- col_a_count + 1
      tableau[i, num_decision_vars + num_constraints + col_a_count] <- 1
      base[i] <- paste0("a", col_a_count)
      original_artificial_indices <- c(original_artificial_indices, num_decision_vars + num_constraints + col_a_count)
    }
  }
  
  colnames(tableau) <- var_names[1:total_cols]
  rhs1 <- rhs
  
  # --- PHASE 1: FIND A FEASIBLE BASIC SOLUTION ---
  
  if (num_artificial > 0) {
    
    # Phase 1 objective function: Minimize the sum of artificial variables.
    # In a maximization problem, this is equivalent to Maximizing W = -Sum(artificials).
    zobj_phase1 <- numeric(total_cols)
    zobj_phase1[original_artificial_indices] <- 1  
    
    z_rhs_phase1 <- 0 # Stores the right-hand side (RHS) value of the Phase 1 objective function
    
    # Ensure the objective function is consistent with the initial basis.
    # Zero out the coefficients of the basic variables in the objective function row.
    for (i in 1:num_constraints) {
      if (grepl("^a", base[i])) {
        z_rhs_phase1 <- z_rhs_phase1 - rhs1[i]
        zobj_phase1 <- zobj_phase1 - tableau[i, ]
      }
    }
    
    # Iterate until there are no more negative coefficients in the objective function
    negative_phase1 <- TRUE
    while (negative_phase1) {
      
      # If there are no negative coefficients, Phase 1 optimality is reached
      if (min(zobj_phase1) >= -1e-9) {
        negative_phase1 <- FALSE
        next
      }
      
      # Find the pivot column (entering variable)
      pivot_col <- which.min(zobj_phase1)
      
      # Minimum ratio test to find the pivot row (leaving variable)
      ratio_values <- rep(NA, num_constraints)
      for (i in 1:num_constraints) {
        if (tableau[i, pivot_col] > 1e-9) { # The divisor must be positive
          ratio_values[i] <- rhs1[i] / tableau[i, pivot_col]
        }
      }
      
      # If all coefficients in the pivot column are non-positive, the problem is unbounded
      if (all(is.na(ratio_values))) {
        return("UNBOUNDED PROBLEM (in Phase 1)")
      }
      
      pivot_row <- which.min(ratio_values)
      
      # Pivoting
      pivot_val <- tableau[pivot_row, pivot_col]
      z_factor <- zobj_phase1[pivot_col]
      
      # Update the objective function value
      z_rhs_phase1 <- z_rhs_phase1 - z_factor * (rhs1[pivot_row] / pivot_val)
      
      # Normalize the pivot row
      tableau[pivot_row, ] <- tableau[pivot_row, ] / pivot_val
      rhs1[pivot_row] <- rhs1[pivot_row] / pivot_val
      
      # Zero out the other elements of the pivot column
      for (i in 1:num_constraints) {
        if (i != pivot_row) {
          factor <- tableau[i, pivot_col]
          rhs1[i] <- rhs1[i] - factor * rhs1[pivot_row]
          tableau[i, ] <- tableau[i, ] - factor * tableau[pivot_row, ]
        }
      }
      
      # Update the objective function row
      zobj_phase1 <- zobj_phase1 - z_factor * tableau[pivot_row, ]
      base[pivot_row] <- colnames(tableau)[pivot_col]
    }
    
    # If the optimal value of Phase 1 is negative, it means the artificial vars could not be zeroed.
    # Therefore, there is no feasible solution for the original problem.
    if (z_rhs_phase1 < -1e-9) {
      return("ERROR: INFEASIBLE PROBLEM")
    }
    
    # Remove the artificial variable columns to start Phase 2
    tableau <- tableau[, -original_artificial_indices, drop = FALSE]
  }
  
  # --- PHASE 2: FIND THE OPTIMAL SOLUTION ---
  
  # Set up the original objective function of the problem (Max Z -> Z - C'X = 0)
  zobj <- numeric(ncol(tableau))
  zobj[1:num_decision_vars] <- -obj
  
  z_rhs <- 0 # The initial value of the original objective function is 0
  
  # Ensure the objective function is consistent with the current basis (post-Phase 1)
  for (i in 1:length(base)) {
    basic_var <- base[i]
    col_idx <- which(colnames(tableau) == basic_var)
    
    if (length(col_idx) > 0 && zobj[col_idx] != 0) {
      factor <- zobj[col_idx]
      z_rhs <- z_rhs - factor * rhs1[i]
      zobj <- zobj - factor * tableau[i, ]
    }
  }
  
  # Start the iterative process of the Simplex method
  negative <- TRUE
  while (negative) {
    # If there are no negative coefficients in the Z row, the optimal solution has been found
    if (min(zobj) >= -1e-9) {
      negative <- FALSE
      next
    }
    
    # Find the pivot column
    pivot_col <- which.min(zobj)
    
    # Minimum ratio test to find the pivot row
    ratio_values <- rep(NA, num_constraints)
    for (i in 1:num_constraints) {
      if (tableau[i, pivot_col] > 1e-9) {
        ratio_values[i] <- rhs1[i] / tableau[i, pivot_col]
      }
    }
    
    # If there is no candidate to leave the basis, the problem is unbounded
    if (all(is.na(ratio_values))) {
      return("UNBOUNDED PROBLEM: The objective function can be maximized infinitely.")
    }
    
    pivot_row <- which.min(ratio_values)
    
    # Pivoting
    pivot_val <- tableau[pivot_row, pivot_col]
    z_factor <- zobj[pivot_col]
    
    # Update the objective function value
    z_rhs <- z_rhs - z_factor * (rhs1[pivot_row] / pivot_val)
    
    # Normalize the pivot row
    tableau[pivot_row, ] <- tableau[pivot_row, ] / pivot_val
    rhs1[pivot_row] <- rhs1[pivot_row] / pivot_val
    
    # Update the other rows of the tableau
    for (i in 1:num_constraints) {
      if (i != pivot_row) {
        factor <- tableau[i, pivot_col]
        rhs1[i] <- rhs1[i] - factor * rhs1[pivot_row]
        tableau[i, ] <- tableau[i, ] - factor * tableau[pivot_row, ]
      }
    }
    
    # Update the objective function row
    zobj <- zobj - z_factor * tableau[pivot_row, ]
    base[pivot_row] <- colnames(tableau)[pivot_col]
  }
  
  # --- RESULT ---
  
  # Extract the final solution from the tableau
  solution <- rep(0, num_decision_vars)
  for (i in 1:length(base)) {
    # If a decision variable is in the basis, its value is the corresponding RHS
    if (grepl("^x", base[i])) {
      var_num <- as.numeric(sub("^x", "", base[i]))
      solution[var_num] <- rhs1[i]
    }
  }
  
  # The optimal value of the objective function is the final value of the Z row's RHS
  maxval <- z_rhs
  
  return(list(solution = solution, optimal_value = maxval))
}