lim_polynomial <- function(x, scale01=FALSE){
  x1 <- x[1]
  x2 <- x[2]

  term1 <- (5/2)*x1 - (35/2)*x2
  term2 <- (5/2)*x1*x2 + 19*x2^2
  term3 <- -(15/2)*x1^3 - (5/2)*x1*x2^2
  term4 <- -(11/2)*x2^4 + (x1^3)*(x2^2)

  y <- 9 + term1 + term2 + term3 + term4
  return(y)
}

grlee2 <- function(x, scale01=FALSE){
  if(scale01){
    x[1:2] <- 8*x[1:2] - 2
  }
  y = x[1]*exp(-x[1]^2-x[2]^2)
  return(y)
}

twin_galaxies <- function(x, scale01=FALSE){
  if(!scale01){
    x[1] <- x[1]/360
    x[2] <- (x[2] + 90)/180
  }

  y <- 22/40*lim_polynomial(x, scale01=TRUE)
  y <- y + 5*grlee2(x, scale01=TRUE)
  return(y)
}
