library(deepgp)
library(lhs)
library(tictoc)
source("R/piston.R")

stop("Warning: This script takes a very long time to run. Comment or remove line 6 to run.")

N_vec <- c(5, 10, 20, 40, 100)
T_vec <- rep(NA, length(N_vec))
for(i in seq_along(N_vec)){
  N <- N_vec[i]
  X <- maximinLHS(N, 7)
  y <- apply(X, 1, piston, scale01=TRUE)
  z <- (y - mean(y))/sd(y)
  tic()
  fit <- fit_two_layer(X, z)
  toc <- toc()
  T_vec[i] <- toc$toc - toc$tic
}

lin_fit <- lm(T_vec ~ I(N_vec) + I(N_vec^2) + I(N_vec^3))
plot(N_vec, T_vec, pch=16, cex=2)
beta <- lin_fit$coefficients
curve(beta[1] + beta[2]*x + beta[3]*x^2 + beta[4]*x^3, add=TRUE)

# Check prediction for N=150
N <- 150
beta[1] + beta[2]*N + beta[3]*N^2 + beta[4]*N^3
X <- maximinLHS(N, 7)
y <- apply(X, 1, piston, scale01=TRUE)
z <- (y - mean(y))/sd(y)
tic()
fit <- fit_two_layer(X, z)
toc <- toc()

#Pretty close, but we'll adjust and calculate again
T_vec <- c(T_vec, toc$toc - toc$tic)
N_vec <- c(N_vec, N)
lin_fit <- lm(T_vec ~ I(N_vec) + I(N_vec^2) + I(N_vec^3))
plot(N_vec, T_vec, pch=16, cex=2)
beta <- lin_fit$coefficients
curve(beta[1] + beta[2]*x + beta[3]*x^2 + beta[4]*x^3, add=TRUE)

#Estimate time (days) to run fit_two_layer on a N=40,000 problem
N <- 40000
(beta[1] + beta[2]*N + beta[3]*N^2 + beta[4]*N^3)/60/60/24
