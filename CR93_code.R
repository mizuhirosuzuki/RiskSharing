rm(list = ls())
library(knitr)

# Common settings ------------------
y_vals <- c(1, 2, 3)
theta_hat_vals <- c((y_vals[2] - y_vals[1]) / 2, 
                    (y_vals[3] - y_vals[1]) / 2, 
                    (y_vals[3] - y_vals[2]) / 2)

# v(theta)
v_theta <- function(theta, y_vals, pi, rho){
  output <- 0
  for (i in 1:3){
    if (i > 1){
      for (j in 1:(i-1)){
        output <- output + pi[i,j] * ((y_vals[i] - theta[i+j-2]) ^ (1 - rho) / (1 - rho) + 
                                      (y_vals[j] + theta[i+j-2]) ^ (1 - rho) / (1 - rho))
      } 
    } else {
      output <- output
    }
    output <- output + pi[i, i] * (y_vals[i]) ^ (1 - rho) / (1 - rho) 
  }
  return(output)
}

# Iteration for theta_star
iteration <- function(theta_star, theta_hat, y_val, y_vals, pi, rho, r, v_bar){
  new_theta_star_ij <- min(theta_hat, 
                           y_val - (y_val ^ (1 - rho) - 
                                    (v_theta(theta_star, y_vals, pi, rho) - v_bar) 
                                    * (1 - rho) / r) ^ (1 / (1 - rho)))
  return(new_theta_star_ij)
}

# Iteration loop for theta_star
theta_star_loop <- function(init_theta_star, theta_hat_vals, y_vals, pi, rho, r, v_bar){
  iter <- 0
  tol <- 1e-16
  error <- tol + 1
  theta_star <- init_theta_star
  new_theta_star <- rep(0, 3)
  while (iter < 500 & error > tol){
    new_theta_star[1] <- iteration(theta_star, theta_hat_vals[1], 
                                   y_vals[2], y_vals, pi, rho, r, v_bar)
    new_theta_star[2] <- iteration(theta_star, theta_hat_vals[2], 
                                   y_vals[3], y_vals, pi, rho, r, v_bar)
    new_theta_star[3] <- iteration(theta_star, theta_hat_vals[3], 
                                   y_vals[3], y_vals, pi, rho, r, v_bar)
    error = max(abs(new_theta_star - theta_star))
    theta_star <- new_theta_star
    iter <- iter + 1
  }
  return(theta_star)
}

# Table 1 -----------------------
# The 'highly covariate' income stream
table_1 <- matrix(rep(0, (5 + 11 + 14) * 6), ncol = 6)

pi <- matrix(c(0.2, 0.05, 0.05, 0.05, 0.3, 0.05, 0.05, 0.05, 0.2), nrow = 3) 
r_vals <- c(0.05, 0.15, 0.25)
rho_vals <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7)

# r = 0.05
r <- r_vals[1]
for (i in 1:5){
  # r
  table_1[i, 1] <- r
  # rho
  table_1[i, 2] <- rho_vals[i]
  # v_bar
  v_bar <- 0
  for (j in 1:3){
    for (k in 1:3){
      v_bar <- v_bar + pi[j,k] * (y_vals[j] ^ (1 - rho_vals[i]) / (1 - rho_vals[i]))
    }
  }
  # theta_star iteration
  init_theta_star <- theta_hat_vals
  theta_star <- theta_star_loop(init_theta_star, theta_hat_vals, 
                                y_vals, pi, rho_vals[i], r, v_bar)
  # theta_star_21
  table_1[i, 3] <- theta_star[1]
  # theta_star_31
  table_1[i, 4] <- theta_star[2]
  # theta_star_32
  table_1[i, 5] <- theta_star[3]
  # gamma
  table_1[i, 6] <- (v_theta(theta_star, y_vals, pi, rho_vals[i]) - v_bar) / 
                   (v_theta(theta_hat_vals, y_vals, pi, rho_vals[i]) - v_bar)
}

# r = 0.15
r <- r_vals[2]
for (i in 1:11){
  # r
  table_1[i+5, 1] <- r
  # rho
  table_1[i+5, 2] <- rho_vals[i]
  # v_bar
  v_bar <- 0
  for (j in 1:3){
    for (k in 1:3){
      v_bar <- v_bar + pi[j,k] * (y_vals[j] ^ (1 - rho_vals[i]) / (1 - rho_vals[i]))
    }
  }
  # theta_star iteration
  init_theta_star <- theta_hat_vals
  theta_star <- theta_star_loop(init_theta_star, theta_hat_vals, 
                                y_vals, pi, rho_vals[i], r, v_bar)
  # theta_star_21
  table_1[i+5, 3] <- theta_star[1]
  # theta_star_31
  table_1[i+5, 4] <- theta_star[2]
  # theta_star_32
  table_1[i+5, 5] <- theta_star[3]
  # gamma
  table_1[i+5, 6] <- (v_theta(theta_star, y_vals, pi, rho_vals[i]) - v_bar) / 
                     (v_theta(theta_hat_vals, y_vals, pi, rho_vals[i]) - v_bar)
}

# r = 0.25
r <- r_vals[3]
for (i in 1:14){
  # r
  table_1[i+16, 1] <- r
  # rho
  table_1[i+16, 2] <- rho_vals[i]
  # v_bar
  v_bar <- 0
  for (j in 1:3){
    for (k in 1:3){
      v_bar <- v_bar + pi[j,k] * (y_vals[j] ^ (1 - rho_vals[i]) / (1 - rho_vals[i]))
    }
  }
  # theta_star iteration
  init_theta_star <- theta_hat_vals
  theta_star <- theta_star_loop(init_theta_star, theta_hat_vals, 
                                y_vals, pi, rho_vals[i], r, v_bar)
  # theta_star_21
  table_1[i+16, 3] <- theta_star[1]
  # theta_star_31
  table_1[i+16, 4] <- theta_star[2]
  # theta_star_32
  table_1[i+16, 5] <- theta_star[3]
  # gamma
  table_1[i+16, 6] <- (v_theta(theta_star, y_vals, pi, rho_vals[i]) - v_bar) / 
                      (v_theta(theta_hat_vals, y_vals, pi, rho_vals[i]) - v_bar)
}

# Table 2 -----------------------
# The 'moderately covariate' income stream
table_2 <- matrix(rep(0, (5 + 7 + 9) * 6), ncol = 6)

pi <- matrix(c(0.1, 0.1, 0.1, 0.1, 0.2, 0.1, 0.1, 0.1, 0.1), nrow = 3) 
r_vals <- c(0.05, 0.15, 0.25)
rho_vals <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7)

# r = 0.05
r <- r_vals[1]
for (i in 1:5){
  # r
  table_2[i, 1] <- r
  # rho
  table_2[i, 2] <- rho_vals[i]
  # v_bar
  v_bar <- 0
  for (j in 1:3){
    for (k in 1:3){
      v_bar <- v_bar + pi[j,k] * (y_vals[j] ^ (1 - rho_vals[i]) / (1 - rho_vals[i]))
    }
  }
  # theta_star iteration
  init_theta_star <- theta_hat_vals
  theta_star <- theta_star_loop(init_theta_star, theta_hat_vals, 
                                y_vals, pi, rho_vals[i], r, v_bar)
  # theta_star_21
  table_2[i, 3] <- theta_star[1]
  # theta_star_31
  table_2[i, 4] <- theta_star[2]
  # theta_star_32
  table_2[i, 5] <- theta_star[3]
  # gamma
  table_2[i, 6] <- (v_theta(theta_star, y_vals, pi, rho_vals[i]) - v_bar) / 
                   (v_theta(theta_hat_vals, y_vals, pi, rho_vals[i]) - v_bar)
}

# r = 0.15
r <- r_vals[2]
for (i in 1:7){
  # r
  table_2[i+5, 1] <- r
  # rho
  table_2[i+5, 2] <- rho_vals[i]
  # v_bar
  v_bar <- 0
  for (j in 1:3){
    for (k in 1:3){
      v_bar <- v_bar + pi[j,k] * (y_vals[j] ^ (1 - rho_vals[i]) / (1 - rho_vals[i]))
    }
  }
  # theta_star iteration
  init_theta_star <- theta_hat_vals
  theta_star <- theta_star_loop(init_theta_star, theta_hat_vals, 
                                y_vals, pi, rho_vals[i], r, v_bar)
  # theta_star_21
  table_2[i+5, 3] <- theta_star[1]
  # theta_star_31
  table_2[i+5, 4] <- theta_star[2]
  # theta_star_32
  table_2[i+5, 5] <- theta_star[3]
  # gamma
  table_2[i+5, 6] <- (v_theta(theta_star, y_vals, pi, rho_vals[i]) - v_bar) / 
                     (v_theta(theta_hat_vals, y_vals, pi, rho_vals[i]) - v_bar)
}

# r = 0.25
r <- r_vals[3]
for (i in 1:9){
  # r
  table_2[i+12, 1] <- r
  # rho
  table_2[i+12, 2] <- rho_vals[i]
  # v_bar
  v_bar <- 0
  for (j in 1:3){
    for (k in 1:3){
      v_bar <- v_bar + pi[j,k] * (y_vals[j] ^ (1 - rho_vals[i]) / (1 - rho_vals[i]))
    }
  }
  # theta_star iteration
  init_theta_star <- theta_hat_vals
  theta_star <- theta_star_loop(init_theta_star, theta_hat_vals, 
                                y_vals, pi, rho_vals[i], r, v_bar)
  # theta_star_21
  table_2[i+12, 3] <- theta_star[1]
  # theta_star_31
  table_2[i+12, 4] <- theta_star[2]
  # theta_star_32
  table_2[i+12, 5] <- theta_star[3]
  # gamma
  table_2[i+12, 6] <- (v_theta(theta_star, y_vals, pi, rho_vals[i]) - v_bar) / 
                      (v_theta(theta_hat_vals, y_vals, pi, rho_vals[i]) - v_bar)
}

# Table 3 -----------------------
# The 'highly non-covariate' income stream
table_3 <- matrix(rep(0, (5 + 6 + 8) * 6), ncol = 6)

pi <- matrix(c(0.1, 0.1, 0.1, 0.1, 0.2, 0.1, 0.1, 0.1, 0.1), nrow = 3) 
r_vals <- c(0.05, 0.15, 0.25)
rho_vals <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7)

# r = 0.05
r <- r_vals[1]
for (i in 1:5){
  # r
  table_3[i, 1] <- r
  # rho
  table_3[i, 2] <- rho_vals[i]
  # v_bar
  v_bar <- 0
  for (j in 1:3){
    for (k in 1:3){
      v_bar <- v_bar + pi[j,k] * (y_vals[j] ^ (1 - rho_vals[i]) / (1 - rho_vals[i]))
    }
  }
  # theta_star iteration
  init_theta_star <- theta_hat_vals
  theta_star <- theta_star_loop(init_theta_star, theta_hat_vals, 
                                y_vals, pi, rho_vals[i], r, v_bar)
  # theta_star_21
  table_3[i, 3] <- theta_star[1]
  # theta_star_31
  table_3[i, 4] <- theta_star[2]
  # theta_star_32
  table_3[i, 5] <- theta_star[3]
  # gamma
  table_3[i, 6] <- (v_theta(theta_star, y_vals, pi, rho_vals[i]) - v_bar) / 
                   (v_theta(theta_hat_vals, y_vals, pi, rho_vals[i]) - v_bar)
}

# r = 0.15
r <- r_vals[2]
for (i in 1:6){
  # r
  table_3[i+5, 1] <- r
  # rho
  table_3[i+5, 2] <- rho_vals[i]
  # v_bar
  v_bar <- 0
  for (j in 1:3){
    for (k in 1:3){
      v_bar <- v_bar + pi[j,k] * (y_vals[j] ^ (1 - rho_vals[i]) / (1 - rho_vals[i]))
    }
  }
  # theta_star iteration
  init_theta_star <- theta_hat_vals
  theta_star <- theta_star_loop(init_theta_star, theta_hat_vals, 
                                y_vals, pi, rho_vals[i], r, v_bar)
  # theta_star_21
  table_3[i+5, 3] <- theta_star[1]
  # theta_star_31
  table_3[i+5, 4] <- theta_star[2]
  # theta_star_32
  table_3[i+5, 5] <- theta_star[3]
  # gamma
  table_3[i+5, 6] <- (v_theta(theta_star, y_vals, pi, rho_vals[i]) - v_bar) / 
                     (v_theta(theta_hat_vals, y_vals, pi, rho_vals[i]) - v_bar)
}

# r = 0.25
r <- r_vals[3]
for (i in 1:8){
  # r
  table_3[i+11, 1] <- r
  # rho
  table_3[i+11, 2] <- rho_vals[i]
  # v_bar
  v_bar <- 0
  for (j in 1:3){
    for (k in 1:3){
      v_bar <- v_bar + pi[j,k] * (y_vals[j] ^ (1 - rho_vals[i]) / (1 - rho_vals[i]))
    }
  }
  # theta_star iteration
  init_theta_star <- theta_hat_vals
  theta_star <- theta_star_loop(init_theta_star, theta_hat_vals, 
                                y_vals, pi, rho_vals[i], r, v_bar)
  # theta_star_21
  table_3[i+11, 3] <- theta_star[1]
  # theta_star_31
  table_3[i+11, 4] <- theta_star[2]
  # theta_star_32
  table_3[i+11, 5] <- theta_star[3]
  # gamma
  table_3[i+11, 6] <- (v_theta(theta_star, y_vals, pi, rho_vals[i]) - v_bar) / 
                      (v_theta(theta_hat_vals, y_vals, pi, rho_vals[i]) - v_bar)
}

# Generate output tables
colnames(table_1) <- c("$r$", "$\\rho$", "$\\theta_{21}^*$", 
                       "$\\theta_{31}^*$", "$\\theta_{32}^*$", "$\\gamma$")
colnames(table_2) <- c("$r$", "$\\rho$", "$\\theta_{21}^*$", 
                       "$\\theta_{31}^*$", "$\\theta_{32}^*$", "$\\gamma$")
colnames(table_3) <- c("$r$", "$\\rho$", "$\\theta_{21}^*$", 
                       "$\\theta_{31}^*$", "$\\theta_{32}^*$", "$\\gamma$")

kable(table_1, digits = 3, caption = "Highly covariate income streams")
kable(table_2, digits = 3, caption = "Moderately covariate income streams")
kable(table_3, digits = 3, caption = "Highly non-covariate income streams")
