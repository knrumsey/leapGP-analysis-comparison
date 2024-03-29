library(lhs)              # Latin hypercube sampling
library(tictoc)           # Timing
library(quack)            # leapGP: devtools::install_github("knrumsey/quack")
library(laGP)             # laGP
source("R/twin_galaxies.R") # twin_galaxies function

set.seed(11111)
# Simulations using Twin Galaxies function
n <- 1000
X <- smartLHS(n, 2)
y <- apply(X, 1, twin_galaxies, scale01=TRUE)

# Make prediction
X2 <- smartLHS(1000, 2)
y2 <- apply(X2, 1, twin_galaxies, scale01=TRUE)

# 1. laGP
tic()
yhat_la <- rep(NA, 1000)
for(i in 1:1000){
  yhat_la[i] <- laGP(matrix(X2[i,], nrow=1), start=6, end=30, X=X, Z=y)$mean
}
toc_la <- toc()

# 2. leapGP
tic()
leap <- leapGP_build(X, y, H = 20, start=6, n=30)
toc_leap1 <- toc()

tic()
yhat_leap <- rep(NA, 1000)
for(i in 1:1000){
  yhat_leap[i] <- leapGP_predict(leap, X2[i,])
}
toc_leap2 <- toc()


# 3. slapGP(rho=0.99) - with training
rho <- 0.9
slap <- leapGP_synch(leap, rho=rho)
hubs <- slap$hubs
tic()
yhat_slap <- rep(NA, 1000)
for(i in 1:1000){
  emulator <- slapGP(X2[i,], X, y, hubs=hubs,
                     rho=rho, start=6, n=30)
  yhat_slap[i] <- emulator$pred
  hubs <- emulator$hubs
}
toc_slap <- toc()
slap_hubs <- hubs


# 4. slapGP without training
rho_vec <- c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999)
hub_vec <- rep(NA, length(rho_vec))
hubs_list <- list()
yhat_slap2 <- matrix(NA, nrow=length(y2), ncol=length(rho_vec))
toc_slap2 <- matrix(NA, nrow=length(y2)/100, ncol=length(rho_vec))
for(ii in seq_along(rho_vec)){
  rho <- rho_vec[ii]
  hubs <- list()
  tic()
  for(jj in 1:1000){
    emulator <- slapGP(X2[jj,], X, y, hubs=hubs,
                       rho=rho, start=6, n=30)
    yhat_slap2[jj, ii] <- emulator$pred
    hubs <- emulator$hubs
    if((jj %% 100) == 0){
      ttt <- toc()
      toc_slap2[jj/100, ii] <- ttt$toc - ttt$tic
      tic()
    }
  }
  hub_vec[ii] <- length(hubs)
  hubs_list[[ii]] <- hubs
  print(paste0("Finished rho = ", rho))
}

