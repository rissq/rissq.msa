d2 <- function (n = NA){

  if (is.na(n) | n < 2 | abs(n-round(n)) != 0){
    stop("Invalid sample size (", n, ")")
  }

  f <- function(x){
    1 - ptukey(x, n, Inf)
  }

  d2 <- integrate(f, 0, Inf)

  if (d2$abs.error > 0.001)
    warning("Absolute error after numerical integration greater than 0.001")

  d2 <- d2$value

  return(d2)
}

d3 <- function (n = NA){

  if (is.na(n) | n < 2 | abs(n-round(n))!=0){
    stop("Invalid sample size")
  }

  f <- function (x){
    x * (1 - ptukey(x, n, Inf))
  }

  d3 <- integrate(f, 0, Inf)

  if (d3$abs.error > 0.001)
    warning("Absolute error after numerical integration greater than 0.001")

  d3 <- 2 * d3$value
  d2 <- d2(n)
  d3 <- sqrt(d3 - d2^2)

  return(d3)
}
