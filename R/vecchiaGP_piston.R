library(tictoc)
library(quack)
library(lhs)
library(GPvecchia)
source("R/piston.R")
rmspe <- function(y, y2) sqrt(mean((y-y2)^2))

# Generate data
set.seed(111)
N <- 4000
X <- quack::smartLHS(N, 7)
y <- apply(X, 1, piston, scale01=TRUE)

# Get data for prediction
N2 <- 10000
X2 <- quack::smartLHS(N2, 7)
y2 <- apply(X2, 1, piston, scale01=TRUE)

# Fit Vecchia GP (Katzfuss et al., 2020)
tic()
fit_vecc <- vecchia_estimate(y, X, m=40) # NOTE: default m=20 does not converge.
toc_vecc <- toc()

# If predictions could be done non-sequentially, vecchia would be faster.
# But it still takes a long time and crashes R script if you try to stop it early.
#pred_batch <- vecchia_pred(fit_vecc, X2)

# Online-predictions
yht_vecc <- cov_vecc <- rep(NA, N2)
tic()
for(i in 1:1000){
  print(i)
  pred <- vecchia_pred(fit_vecc, matrix(X2[i,], nrow=1))
  pred$mean <- pred$mean.pred
  pred$var <- pred$var.pred
  yht_vecc[i] <- pred$mean[1]
  cov_vecc[i] <- (pred$mean[1] - 1.96*sqrt(pred$var[1]) < y2[i]) & (pred$mean[1] + 1.96*sqrt(pred$var[1]) > y2[i])
}
toc_vecc2 <- toc()

tmp <- rmspe(y2[1:(i-1)], yht_vecc[1:(i-1)])*60 #multiply by 60 for units of minutes/cycle

N3 <- 1000
cat("Scaled Vecchia GP on piston function with N = ", N,
    "\n   Train (s): ", toc_vecc$toc - toc_vecc$tic,
    "\n   Pred (s): ", (toc_vecc2$toc - toc_vecc2$tic)*N2/N3,
    "\n   RMSPE: ", tmp,
    "\n   Coverage: ", sum(cov_vecc, na.rm=TRUE)/N3)

save(yht_vecc, toc_vecc, toc_vecc2, cov_vecc,
     file=paste0("Data/VeccGP_Piston_results", N, "_1000.rda"))


