######EXAMPLE###########

A<-matrix(nrow = 4,ncol = 5)  #Matrix of constraints (first rows) and objective (last row)
A[1,] = c(2,-1,-1,0,0)    #If >=, then slack variable is -1
A[2,] = c(3,1,0,-1,0)
A[3,] = c(2,-3,0,0,1)     #If <= then slack variable is +1
A[4,] = c(50,-20,0,0,0)

a<-c(-5,3,12,0)   # Vector of constants

#########################

simplex <- function(A, a) {
  obj <- A[nrow(A), ]  # Objective row
  while (min(obj) < 0) {
    col_pivo <- which.min(obj)  # Column with smallest value in objective
    razoes <- a / A[, col_pivo]
    razoes[A[, col_pivo] <= 0] <- Inf  # Ignore division by zero/negatives
    lin_pivo <- which.min(razoes)       # Pivot row
    
    # Update pivot row
    pivo <- A[lin_pivo, col_pivo]
    A[lin_pivo, ] <- A[lin_pivo, ] / pivo
    a[lin_pivo] <- a[lin_pivo] / pivo
    
    # Update other rows
    for (i in 1:nrow(A)) {
      if (i == lin_pivo) next
      fator <- -A[i, col_pivo]
      A[i, ] <- A[i, ] + fator * A[lin_pivo, ]
      a[i] <- a[i] + fator * a[lin_pivo]
    }
    obj <- A[nrow(A), ]  # Update objective
  }
  
  # Extract solution for original variables
  solucao <- numeric(ncol(A) - 1)  # Vector for x1, x2, ...
  for (j in 1:(ncol(A) - 1)) {
    if (sum(A[, j] == 1) == 1 && sum(A[, j]) == 1) {  # If column is pivot
      lin <- which(A[, j] == 1)
      solucao[j] <- a[lin]
    }
  }
  
  return(list(x = solucao, Z = a[nrow(A)]))  #Z is the last element
}

resultado <- simplex(A, a)
print(resultado)