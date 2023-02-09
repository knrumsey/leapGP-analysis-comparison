library(RColorBrewer)
load("Data/my_tab_4000.rda")


#MAKE FIGURES
# =============================================
# RMSPE
png(paste0("Figures/tmp_piston_rmspe_", N, ".png"),
    width=5, height=5.5, res=300, units="in")
par(mar=c(9,4,2.5,2) + 0.1)
tmp <- my_tab[,3]
ord <- order(tmp)
plot(1:length(tmp), tmp[ord], type='h', xaxt='n', ylab='RMSPE', xlab="", col='white',
     xlim=c(.5, length(tmp)+0.5))
axis(1, 1:length(tmp), rownames(my_tab)[ord], las=2)
eps <- 0.4
bob <- rep("black", 15)
bob[c(7:14, 6)] <- brewer.pal(9, "Blues")
for(i in 1:15){
  zz <- tmp[ord[i]]
  rect(i-eps, 0, i+eps, zz, border="black", col=adjustcolor(bob[ord[i]], 0.6), lwd=2)
}
dev.off()

# Timing
png(paste0("Figures/tmp_piston_timing_", N, ".png"),
    width=5, height=5.5, res=300, units="in")
par(mar=c(9,4,2.5,2) + 0.1)
tmp <- my_tab[,2] + my_tab[,1]
ord <- order(my_tab[,2])
plot(1:length(tmp), log10(my_tab[ord,2]), type='p', xaxt='n', ylab='Time (seconds)', xlab="", col='white', xlim=c(.5, length(tmp)+0.5), yaxt='n')
axis(1, 1:length(tmp), rownames(my_tab)[ord], las=2)
eps <- 0.4
bob <- rep("black", 15)
bob[c(7:14, 6)] <- brewer.pal(9, "Blues")
for(i in 1:15){
  zz <- log10(my_tab[ord[i], 2])
  rect(i-eps, log10(min(my_tab[,2])/1.1), i+eps, zz, border="black", col=adjustcolor(bob[ord[i]], 0.6), lwd=2)
  #rect(i-eps, zz, i+eps, log(tmp[ord[i]]), border="black", lwd=2)
}
axis(2, pretty(range(log10(tmp)), n=4), 10^pretty(range(log10(tmp)), n=4))
dev.off()

