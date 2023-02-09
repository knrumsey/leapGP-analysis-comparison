# CAUTION: These simulations can take a really long time for large N

library(stargazer)     # Latex Tables
library(tictoc)        # Timing package
library(RColorBrewer)  # R Colors
library(lhs)           # Design
library(quack)         # leapGP: devtools::install_github("knrumsey/quack")
library(laGP)          # laGP
library(BASS)          # Bayesian MARS
library(BART)          # bart
library(BayesPPR)      # bppr: devtools::install_github("gqcollins/BayesPPR")
library(gplite)        # FITC (inducing point GP)
library(tgp)           # bcart
library(GPvecchia)
source("R/piston.R")     # piston function


# SIMULATION SETTINGS
# Get data for training
set.seed(111)
N <- 4000
X <- quack::smartLHS(N, 7)
y <- apply(X, 1, piston, scale01=TRUE)

# Get data for prediction
N2 <- 10000
X2 <- quack::smartLHS(N2, 7)
y2 <- apply(X2, 1, piston, scale01=TRUE)

# laGP and leapGP parameters
if(N = 4000){
  n <- M0 <- 60
}
if(N=40000){
  n <- 120
  M <- 200
}

# 1. BASS Model
# =============================================
tic()
mod_bass <- bass(X, y)
toc_bass <- toc()

tic()
yht_bass <- cov_bass <- rep(NA, N2)
for(i in 1:N2){
  pred <- predict(mod_bass, X2[i,])
  pred2 <- pred + rnorm(1000, 0, sqrt(mod_bass$s2))
  yht_bass[i] <- mean(pred)
  cov_bass[i] <- (quantile(pred2, 0.025) < y2[i]) & (quantile(pred2, 0.975) > y2[i])
}
toc_bass2 <- toc()

save(toc_bass, toc_bass2, toc_bart, toc_bart2, toc_sgp, toc_sgp2, toc_la, toc_leap, toc_leap2, toc_slap, toc_tgp, toc_tgp2, toc_ppr, toc_ppr2,
     yht_bass, yht_bart, yht_sgp, yht_la, yht_leap, yht_slap, yht_tgp, yht_ppr,
     cov_bass, cov_bart, cov_sgp, cov_la, cov_leap, cov_slap, cov_tgp, cov_ppr,
     file=paste0("Data/piston_sim_output_N", N, ".rda"))

# 2. BART Model
# =============================================
tic()
mod_bart <- wbart(X, y)
toc_bart <- toc()

# Note, because BART is so slow at serial predictions,
# we only run this for the first 1000 predictions
tic()
yht_bart <- cov_bart <- rep(NA, 1000)
for(i in 1:1000){
  pred <- predict(mod_bart, matrix(X2[i,], nrow=1))
  pred2 <- pred + rnorm(1000, 0, mod_bart$sigma)
  yht_bart[i] <- mean(pred)
  cov_bart[i] <- (quantile(pred2, 0.025) < y2[i]) & (quantile(pred2, 0.975) > y2[i])
}
toc_bart2 <- toc()


save(toc_bass, toc_bass2, toc_bart, toc_bart2, toc_sgp, toc_sgp2, toc_la, toc_leap, toc_leap2, toc_slap, toc_tgp, toc_tgp2, toc_ppr, toc_ppr2,
     yht_bass, yht_bart, yht_sgp, yht_la, yht_leap, yht_slap, yht_tgp, yht_ppr,
     cov_bass, cov_bart, cov_sgp, cov_la, cov_leap, cov_slap, cov_tgp, cov_ppr,
     file=paste0("Data/piston_sim_output_N", N, ".rda"))

# 3. Inducing point GP
# =============================================
tic()
gp <- gp_init(cf_sexp(), method = method_fitc(num_inducing = sqrt(N)))
gp <- gp_optim(gp, X, y)
toc_sgp <- toc()

tic()
yht_sgp <- cov_sgp <- rep(NA, N2)
for(i in 1:1000){
  pred <- gp_pred(gp, matrix(X2[i,], nrow=1), var=TRUE)
  yht_sgp[i] <- pred$mean
  cov_sgp[i] <- (pred$mean - 1.96*sqrt(pred$var) < y2[i]) & (pred$mean + 1.96*sqrt(pred$var) > y2[i])
}
toc_sgp2 <- toc()

save(toc_bass, toc_bass2, toc_bart, toc_bart2, toc_sgp, toc_sgp2, toc_la, toc_leap, toc_leap2, toc_slap, toc_tgp, toc_tgp2, toc_ppr, toc_ppr2,
     yht_bass, yht_bart, yht_sgp, yht_la, yht_leap, yht_slap, yht_tgp, yht_ppr,
     cov_bass, cov_bart, cov_sgp, cov_la, cov_leap, cov_slap, cov_tgp, cov_ppr,
     file=paste0("Data/piston_sim_output_N", N, ".rda"))

# 4. laGP Model
# =============================================
tic()
yht_la <- cov_la <- rep(NA, N2)
for(i in 1:N2){
  mod <- laGP(matrix(X2[i,], nrow=1), start=30, end=n, X, y)
  yht_la[i] <- mod$mean
  cov_la[i] <- (mod$mean - 1.984*sqrt(mod$s2) < y2[i]) & (mod$mean + 1.984*sqrt(mod$s2) > y2[i])
}
toc_la <- toc()

save(toc_bass, toc_bass2, toc_bart, toc_bart2, toc_sgp, toc_sgp2, toc_la, toc_leap, toc_leap2, toc_slap, toc_tgp, toc_tgp2, toc_ppr, toc_ppr2,
     yht_bass, yht_bart, yht_sgp, yht_la, yht_leap, yht_slap, yht_tgp, yht_ppr,
     cov_bass, cov_bart, cov_sgp, cov_la, cov_leap, cov_slap, cov_tgp, cov_ppr,
     file=paste0("Data/piston_sim_output_N", N, ".rda"))

# 5. leapGP Model
# =============================================
tic()
leap <- leapGP_build(X, y, H=M0, start=30, n=n, scale=TRUE)
toc_leap <- toc()

tic()
yht_leap <- cov_leap <- rep(NA, N2)
for(i in 1:N2){
  pred <- leapGP_predict(leap, X2[i,], scale=TRUE)
  yht_leap[i] <- pred[1]
  cov_leap[i] <- (pred[1] - 1.96*sqrt(pred[2]) < y2[i]) & (pred[1] + 1.96*sqrt(pred[2]) > y2[i])
}
toc_leap2 <- toc()


save(toc_bass, toc_bass2, toc_bart, toc_bart2, toc_sgp, toc_sgp2, toc_la, toc_leap, toc_leap2, toc_slap, toc_tgp, toc_tgp2, toc_ppr, toc_ppr2,
     yht_bass, yht_bart, yht_sgp, yht_la, yht_leap, yht_slap, yht_tgp, yht_ppr,
     cov_bass, cov_bart, cov_sgp, cov_la, cov_leap, cov_slap, cov_tgp, cov_ppr,
     file=paste0("Data/piston_sim_output_N", N, ".rda"))

# 6. slapGP Model(s)
# =============================================
rho_vec <- c(0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99)
yht_slap <- cov_slap  <- matrix(NA, nrow=N2, ncol=length(rho_vec))
toc_slap <- rep(NA, length(rho_vec))
for(jj in 1:length(rho_vec)){
  tic()
  rho <- rho_vec[jj]
  slap <- leapGP_synch(leap, rho)
  hubs <- slap$hubs
  for(i in 1:N2){
    emulator <- slapGP(X2[i,], X, y, rho=rho, scale=TRUE, start=30, n=n, hubs=hubs)
    hubs <- emulator$hubs
    m <- emulator$pred
    s <- sqrt(emulator$scale)
    yht_slap[i,jj] <- m
    cov_slap[i,jj] <- (m - 1.96*s < y2[i]) & (m + 1.96*s > y2[i])
  }
  ttt <- toc()
  toc_slap[jj] <- ttt$toc - ttt$tic
}

save(toc_bass, toc_bass2, toc_bart, toc_bart2, toc_sgp, toc_sgp2, toc_la, toc_leap, toc_leap2, toc_slap, toc_tgp, toc_tgp2, toc_ppr, toc_ppr2,
     yht_bass, yht_bart, yht_sgp, yht_la, yht_leap, yht_slap, yht_tgp, yht_ppr,
     cov_bass, cov_bart, cov_sgp, cov_la, cov_leap, cov_slap, cov_tgp, cov_ppr,
     file=paste0("piston_sim_output_N", N, ".rda"))

# 6. Treed CART Model
# =============================================
tic()
mod_tgp <- bcart(X, y, pred.n=FALSE)
toc_tgp <- toc()

tic()
yht_tgp <- cov_tgp <- rep(NA, N2)
for(i in 1:N2){
  pred <- predict(mod_tgp, matrix(X2[i,], nrow=1))
  yht_tgp[i] <- pred$ZZ.mean
  cov_tgp[i] <- (pred$ZZ.mean - 1.96*sqrt(pred$ZZ.ks2) < y2[i]) & (pred$ZZ.mean + 1.96*sqrt(pred$ZZ.ks2) > y2[i])
}
toc_tgp2 <- toc()

save(toc_bass, toc_bass2, toc_bart, toc_bart2, toc_sgp, toc_sgp2, toc_la, toc_leap, toc_leap2, toc_slap, toc_tgp, toc_tgp2, toc_ppr, toc_ppr2,
     yht_bass, yht_bart, yht_sgp, yht_la, yht_leap, yht_slap, yht_tgp, yht_ppr,
     cov_bass, cov_bart, cov_sgp, cov_la, cov_leap, cov_slap, cov_tgp, cov_ppr,
     file=paste0("Data/piston_sim_output_N", N, ".rda"))


# 7. Bayesian PPR (Only use first 1000 for prediction, too slow.)
# =============================================
tic()
mod_ppr <- bppr(X, y)
toc_ppr <- toc()

tic()
yht_ppr <- cov_ppr <- rep(NA, N2)
for(i in 1:1000){
  pred <- predict(mod_ppr, matrix(X2[i,], nrow=1))
  pred2 <- pred + rnorm(1000, 0, mod_ppr$sd_resid)
  yht_ppr[i] <- mean(pred)
  cov_ppr[i] <- (quantile(pred2, 0.025) < y2[i]) & (quantile(pred2, 0.975) > y2[i])
}
toc_ppr2 <- toc()

save(toc_bass, toc_bass2, toc_bart, toc_bart2, toc_sgp, toc_sgp2, toc_la, toc_leap, toc_leap2, toc_slap, toc_tgp, toc_tgp2, toc_ppr, toc_ppr2,
     yht_bass, yht_bart, yht_sgp, yht_la, yht_leap, yht_slap, yht_tgp, yht_ppr,
     cov_bass, cov_bart, cov_sgp, cov_la, cov_leap, cov_slap, cov_tgp, cov_ppr,
     file=paste0("Data/piston_sim_output_N", N, ".rda"))







# MAKE TABLE
# =============================================
rmspe <- function(y, y2) sqrt(mean((y-y2)^2))
my_tab <- matrix(NA, nrow=14, ncol=4)
rownames(my_tab) <- c("BASS", "BART", "FITC", "BCART", "BPPR", "laGP",
                      "leapGP(M0, 0)",
                      "leapGP(M0, 0.80)",
                      "leapGP(M0, 0.90)",
                      "leapGP(M0, 0.95)",
                      "leapGP(M0, 0.96)",
                      "leapGP(M0, 0.97)",
                      "leapGP(M0, 0.98)",
                      "leapGP(M0, 0.99)")
colnames(my_tab) <- c("Train (s)", "Pred (s)", "RMSPE", "Coverage")



my_tab[1,] <- c(toc_bass$toc - toc_bass$tic,
                toc_bass2$toc - toc_bass$tic,
                rmspe(yht_bass, y2),
                sum(cov_bass, na.rm=TRUE)/length(cov_bass))

my_tab[2,] <- c(toc_bart$toc - toc_bart$tic,
                (toc_bart2$toc - toc_bart2$tic)*N2/1000,
                rmspe(yht_bart[1:1000], y2[1:1000]),
                sum(cov_bart, na.rm=TRUE)/length(cov_bart[1:1000]))

my_tab[3,] <- c(toc_sgp$toc - toc_sgp$tic,
                (toc_sgp2$toc - toc_sgp2$tic)*N2/1000,
                rmspe(yht_sgp[1:1000], y2[1:1000]),
                sum(cov_sgp[1:1000], na.rm=TRUE)/length(cov_sgp[1:1000]))

my_tab[4,] <- c(toc_tgp$toc - toc_tgp$tic,
                toc_tgp2$toc - toc_tgp2$tic,
                rmspe(yht_tgp, y2),
                sum(cov_tgp, na.rm=TRUE)/length(cov_sgp))

my_tab[5,] <- c(toc_ppr$toc - toc_ppr$tic,
                (toc_ppr2$toc - toc_ppr2$tic)*N2/1000,
                rmspe(yht_ppr[1:1000], y2[1:1000]),
                sum(cov_ppr, na.rm=TRUE)/length(cov_ppr[1:1000]))


my_tab[6,] <- c(0,
                toc_la$toc - toc_la$tic,
                rmspe(yht_la, y2),
                sum(cov_la, na.rm=TRUE)/length(cov_la))

my_tab[7,] <- c(toc_leap$toc - toc_leap$tic,
                toc_leap2$toc - toc_leap2$tic,
                rmspe(yht_leap, y2),
                sum(cov_leap, na.rm=TRUE)/length(cov_leap))


for(i in 1:7){
  my_tab[7+i,] <- c(toc_leap$toc - toc_leap$tic,
                  toc_slap[i],
                  rmspe(yht_slap[,i], y2),
                  sum(cov_slap[,i], na.rm=TRUE)/nrow(cov_slap))
}


my_tab[,3] <- my_tab[,3]*60 #RMSPE NOW HAS AN INTERPRETATION FOR UNITS OF MINUTES INSTEAD OF SECONDS
save(my_tab, file=paste0("Data/my_tab_", N, ".rda"))


# SAVE LATEX OUTPUT
# =============================================
sink(file=paste0("Data/table_", N, ".txt"))
tmp = stargazer(my_tab, digits=2)
sink()







