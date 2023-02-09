load("Data/my_tab_40000.rda")
load("Data/piston_sim_output_N40000.rda")
load("Data/VeccGP_Piston_results40000_1000.rda")


# Generate data
set.seed(111)
N <- 40000
X <- quack::smartLHS(N, 7)
y <- apply(X, 1, piston, scale01=TRUE)

# Get data for prediction
N2 <- 10000
X2 <- quack::smartLHS(N2, 7)
y2 <- apply(X2, 1, piston, scale01=TRUE)

my_tab <- rbind(my_tab,
                c(toc_vecc$toc - toc_vecc$tic,
                  (toc_vecc2$toc - toc_vecc2$tic)*N2/1000,
                  rmspe(yht_vecc[1:1000], y2[1:1000])*60,
                  mean(cov_vecc[1:1000]))
                )
rownames(my_tab) <- c(rownames(my_tab)[-nrow(my_tab)], "GPvecchia")
my_tab <- my_tab[c(1,2,5,3,15,4,6,7:14),]
rn <- rownames(my_tab)
rv <- c(0, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99)
for(i in 8:15){
  if(N == 4000){
    rn[i] <- paste0("leapGP(60, ", rv[i-7], ")")
  }else{
    if(N == 40000){
      rn[i] <- paste0("leapGP(200, ", rv[i-7], ")")
    }else{
    rn[i] <- paste0("leapGP(M0, ", rv[i-7], ")")
    }
  }
}
rownames(my_tab) <- rn
save(my_tab, file=paste0("Data/my_tab_", N, "v2.Rda"))


library(RColorBrewer)

#MAKE FIGURES
# =============================================
# RMSPE
png(paste0("Figures/piston_rmspe_", N, "v2.png"),
    width=5, height=5.5, res=300, units="in")
par(mar=c(9,4,2.5,2) + 0.1)
tmp <- my_tab[,3]
ord <- order(tmp)
plot(1:length(tmp), tmp[ord], type='h', xaxt='n', ylab='RMSPE', xlab="", col='white',
     xlim=c(.5, length(tmp)+0.5),
     ylim=c(0, min(max(tmp), 5*median(tmp))))
axis(1, 1:length(tmp), rownames(my_tab)[ord], las=2)
eps <- 0.4
bob <- rep("black", 15)
bob[c(8:15, 7)] <- brewer.pal(9, "Blues")
for(i in 1:15){
  zz <- tmp[ord[i]]
  rect(i-eps, 0, i+eps, zz, border="black", col=adjustcolor(bob[ord[i]], 0.6), lwd=2)
}
dev.off()

# Timing
png(paste0("Figures/piston_timing_", N, "v2.png"),
    width=5, height=5.5, res=300, units="in")
par(mar=c(9,4,2.5,2) + 0.1)
tmp <- my_tab[,2] + my_tab[,1]
ord <- order(my_tab[,2])
plot(1:length(tmp), log10(my_tab[ord,2]), type='p', xaxt='n', ylab='Time (seconds)', xlab="", col='white', xlim=c(.5, length(tmp)+0.5), yaxt='n')
axis(1, 1:length(tmp), rownames(my_tab)[ord], las=2)
eps <- 0.4
bob <- rep("black", 15)
bob[c(8:15, 7)] <- brewer.pal(9, "Blues")
for(i in 1:15){
  zz <- log10(my_tab[ord[i], 2])
  rect(i-eps, log10(min(my_tab[,2])/1.1), i+eps, zz, border="black", col=adjustcolor(bob[ord[i]], 0.6), lwd=2)
  #rect(i-eps, zz, i+eps, log(tmp[ord[i]]), border="black", lwd=2)
}
axis(2, pretty(range(log10(tmp)), n=4), 10^pretty(range(log10(tmp)), n=4))
dev.off()
