piston <- function(x, scale01=FALSE){
  if(scale01){
    RR <- cbind(c(30, 0.005, 0.002, 1000, 90000, 290, 340),
                c(60, 0.020, 0.010, 5000, 110000, 296, 360))
    x[1:7] <- x[1:7]*(RR[,2] - RR[,1]) + RR[,1]
  }

  M  <- x[1]
  S  <- x[2]
  V0 <- x[3]
  k  <- x[4]
  P0 <- x[5]
  Ta <- x[6]
  T0 <- x[7]

  Aterm1 <- P0 * S
  Aterm2 <- 19.62 * M
  Aterm3 <- -k*V0 / S
  A <- Aterm1 + Aterm2 + Aterm3

  Vfact1 <- S / (2*k)
  Vfact2 <- sqrt(A^2 + 4*k*(P0*V0/T0)*Ta)
  V <- Vfact1 * (Vfact2 - A)

  fact1 <- M
  fact2 <- k + (S^2)*(P0*V0/T0)*(Ta/(V^2))

  C <- 2 * pi * sqrt(fact1/fact2)
  return(C)
}
