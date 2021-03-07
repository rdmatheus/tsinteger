## Distributions: ----

### Binomial -------------------------------------------------------------------
binomial_t <- function(){
  out <- list()

  ## Probability mass function
  out$d <- function(x, alpha){
    stats::dbinom(x, size = 1, prob = alpha)
  }

  ## Random generator
  out$r <- function(n, alpha){
    stats::rbinom(n, size = 1, prob = alpha)
  }

  ## Name
  out$name <- "Binomial"

  ## Variance function
  out$var <- function(alpha) alpha * (1 - alpha)

  out
}

### Negative-binomial ----------------------------------------------------------
nbinomial_t <- function(){
  out <- list()

  ## Probability mass function
  out$d <- function(x, alpha){
    stats::dgeom(x, prob = 1/(alpha + 1))
  }

  ## Random generator
  out$r <- function(n, alpha){
    stats::rgeom(n, prob =  1/(alpha + 1))
  }

  ## Name
  out$name <- "Negative binomial"

  ## Variance function
  out$var <- function(alpha) alpha * (1 + alpha)

  out
}

### Poisson --------------------------------------------------------------------
poisson_t <- function(){
  out <- list()

  ## Probability mass function
  out$d <- function(x, alpha){
    stats::dpois(x, lambda = alpha)
  }

  ## Random generator
  out$r <- function(n, alpha){
    stats::rpois(n, lambda = alpha)
  }

  ## Name
  out$name <- "Poisson"

  ## Variance function
  out$var <- function(alpha) alpha

  out
}
