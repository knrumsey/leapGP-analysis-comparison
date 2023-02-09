source("R/twin_galaxies_analysis.R")    # Run twin_galaxies_analysis.R first
source("R/image_scale.R")                   # Color bar for image
library(RColorBrewer)                     # R Colors
library(plotrix)                          # For plotting circles

# Arrange data for plotting with image
yy <- matrix(0, nrow=201, ncol=201)
i <- 1
for(x1 in seq(0, 1, length.out=201)){
  j <- 1
  for(x2 in seq(0, 1, length.out=201)){
    yy[i,j] <- twin_galaxies(c(x1, x2), scale01=TRUE)
    j <- j+1
  }
  i <- i + 1
}

# Make Figure 1a
png("Figures/twin_galaxies.png", height=4.5, width=5, units="in", res=300)
layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(6, 1), heights=c(1,1))
layout.show(4)
par(mar=c(4,4,1,0) + 0.25)
image(yy, col=colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
      xlab=expression("x"[1]), ylab=expression("x"[2]), font=4)
par(mar=c(4,0,1,2) + 0.25)
image.scale(yy, c(0, 6), xlab="",
            col=colorRampPalette(rev(brewer.pal(11, "RdBu")))(100), horiz=FALSE)
dev.off()


# Make Figure 2a (note: the version in the paper is run with (M0=20), for illustrative purposes)
set.seed(111)
rho0 = 0.92
leap0 <- leapGP_build(X, y, H = 20, start=6, n=30)
slap0 <- leapGP_synch(leap0, rho=rho0)
hubs0 <- slap0$hubs
for(i in 1:1000){
  emulator0 <- slapGP(X2[i,], X, y, hubs=hubs0,
                     rho=rho0, start=6, n=30)
  hubs0 <- emulator0$hubs
}

png("Figures/twin_galaxies_hub_placement.png", height=5, width=5, units="in", res=300)
par(mfrow=c(1,1), mar=c(5,4,4,2)+0.1)
image(yy, col=adjustcolor(colorRampPalette(rev(brewer.pal(11, "RdBu")))(100), alpha.f=0.6),
      xlab=expression("x"[1]), ylab=expression("x"[2]), font=4)
rho = 0.995
bob <- rep(brewer.pal(8, "Dark2"), 10)
points(X, pch=16, col='darkgrey', cex=0.5)
for(i in 1:length(leap0$hubs)){
  hh <- leap0$hubs[[i]]
  xx <- X[hh$coord_id,]
  rr <- sqrt(-1/hh$kappa*log(rho))
  points(xx[1], xx[2], col=bob[1], pch=4, cex=2, lwd=3)
  draw.circle(xx[1], xx[2], radius = rr, border=bob[1], lwd=3, col=adjustcolor(bob[1], alpha.f=0.2))
}
my_slap_hubs <- hubs_list[[5]]
for(i in 21:length(hubs0)){
  hh <- hubs0[[i]]
  xx <- hh$coord
  rr <- sqrt(-1/hh$kappa*log(rho))
  points(xx[1], xx[2], col=bob[2], pch=3, cex=2, lwd=3)
  draw.circle(xx[1], xx[2], radius = rr, border=bob[2], lwd=3, col=adjustcolor(bob[2], alpha.f=0.2))
}
dev.off()


# Figure 2a
png("Figures/twin_galaxies_timing.png", height=5, width=5, units="in", res=300)
toc_slap2 <- toc_slap2[, c(1, 4,5,6,7,8)] # Reduce rho values plotted to remove clutter and redundancy
bob <- brewer.pal(8, "Set1")
bob[6] <- bob[8]
plot(NULL, xlim=c(0, 1000), ylim=c(0, 70), xlab='Number of Predictions', ylab='Time (seconds)')
abline(0, (toc_la$toc - toc_la$tic)/1000, lwd=2, lty=3, col=bob[1])
for(i in 1:ncol(toc_slap2)){
  lines(seq(0, 1000, by=100), c(0,cumsum(toc_slap2[,i])), col=bob[1+i])
  points(seq(0, 1000, by=100), c(0,cumsum(toc_slap2[,i])), col=bob[1+i],
         pch=i-1, cex=1.2)
}
legend("topleft", c("laGP",
                 expression(paste("leapGP(", rho, "=0.50)")),
                 expression(paste("leapGP(", rho, "=0.80)")),
                 expression(paste("leapGP(", rho, "=0.90)")),
                 expression(paste("leapGP(", rho, "=0.95)")),
                 expression(paste("leapGP(", rho, "=0.99)")),
                 expression(paste("leapGP(", rho, "=0.999)"))),
       col=bob, lwd=2, lty=c(3, rep(1, 6)), pch=-1:5, cex=0.9)
dev.off()


# Figure 2b
png("Figures/twin_galaxies_results.png", height=5, width=5, units="in", res=300)
rmspe_vals <- log(apply(yhat_slap2, 2, function(zz) sqrt(sum((y2-zz)^2))))
plot(NULL, xlim=c(0.5, 1.0), ylim=range(rmspe_vals),
     xlab=expression(rho), ylab='RMSPE', yaxt="n")
abline(h=log(sqrt(sum((y2-yhat_la)^2))), lwd=2, lty=3, col=bob[1])
lines(c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999), rmspe_vals,
      lwd=2)
points(c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999), rmspe_vals,
      pch=16)
legend("topright", c("laGP", expression(paste("leapGP(",rho,")"))),
       col=c("red", "black"), lwd=2, lty=c(3,1), pch=c(-1, 16),
       cex=1.2)
axis(2, log(c(0.05, 0.1,  0.2, 0.5, 1.0)), c(0.05, 0.1,  0.2, 0.5, 1.0))
dev.off()


