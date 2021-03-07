#' @name innovation
#'
#' @title Family of Innovations Process for Fitting an INAR(p) Model
#'
#' @description Provide the current available distributions that can be
#'     used as an innovation process in a fit of the INAR(p) model.
#'
#' @param dist Character specification of the innovation process, see
#'     details.
#' @param x An \code{"innovation"} object.
#' @param ... Further arguments for other specific methods.
#'
#' @details There are some discrete distributions available for the
#'     innovation process specification. The following table display their
#'     names and their abbreviations to be passed to \code{innovation()}.
#'
#'   \tabular{lll}{
#'  \bold{Distribution}  \tab \bold{Abbreviation} \tab \bold{Number of parameters}\cr
#'  Bernoulli \tab \code{"BE"} \tab 1 \cr
#'  BerPoi    \tab \code{"BP"} \tab 2 \cr
#'  BerG    \tab \code{"BG"} \tab 2 \cr
#'  Mean-Parameterized COM-Poisson \tab \code{"CP"} \tab 2 \cr
#'  Geometric \tab \code{"GE"} \tab 1 \cr
#'  Poisson  \tab \code{"PO"}      \tab  1  \cr
#'  }
#'
#'
#' @return The function \code{innovation()} returns a list whose components
#'     are set of functions referring to the innovation process. More
#'     specifically, returns an \code{"innovation"} object with the following
#'     elements:
#'  \itemize{
#'    \item{d:}{ Probability mass function \code{function(x, mu, phi)}.}
#'    \item{r:}{ Random generator function \code{function(n, mu, phi)}.}
#'    \item{par:}{ Describes the parametric space of the distribution.}
#'    \item{npar:}{ Number of parameters of the innovation process specified.}
#'    \item{constraint:}{ Character; describes the restriction between
#'        parameters if any.}
#'    \item{disp:}{ Describes the types of dispersion that the distribution
#'        can model.}
#'    \item{name:}{ Name of the distribution.}
#'  }
#'
#'  The arguments of the returned functions are:
#'  \describe{
#'    \item{\code{x}}{ Vector of discrete non-negative quantiles.}
#'    \item{\code{n}}{ Number of observations to return.}
#'    \item{\code{mu}}{ Mean parameter of the innovation process.}
#'    \item{\code{phi}}{ Dispersion parameter of the distribution, if it
#'     does not exist, it receives \code{NULL}.}
#'  }
#'
#'
#' @author Rodrigo M. R. Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @examples
#' \dontrun{
#'  ### Specification of the Poisson innovation to 'inv' object
#'  inv <- innovation("PO")
#'
#'  ### Methods
#'  inv
#'  plot(inv)
#'
#'  ### Generating observations
#'  x <- inv$r(500, 5)
#'
#'  ### Barplot and probability mass function
#'  xaxis <- barplot(prop.table(table(x)), main = inv$name,
#'                   xlab = "x", ylab = "Proportion")
#'  points(xaxis, inv$d(sort(unique(x)), 5),
#'         type = "b", pch = 16, col = 2)
#' }
#'
#' @export
innovation <- function(dist){
  dist <- get(dist)()
  class(dist) <- "innovation"
  dist
}

#' @rdname innovation
#' @export
print.innovation <- function(x, ...){

  innovations <- c("BE", "BP", "BG", "GE", "PO")
  INNOVATIONS <- c("Bernoulli", "BerPoi", "BerG", "Geometric", "Poisson")

  cat("----------------------------------------",
      "\nInnovation process:", x$name,
      "\n----------------------------------------",
      "\nAbreviation:", innovations[INNOVATIONS == x$name],
      "\nParameters:", x$par)

      if (x$npar > 1){
       cat("\nConstraints:", x$constraint)
      }
      cat("\nDispersions:", x$disp)
      cat("\n---")
}

### Distributions: ----

### Bernoulli ------------------------------------------------------------------
BE <- function(){
  out <- list()

  ## Probability mass function
  out$d <- function(x, par){
    stats::dbinom(x, size = 1, prob = par)
  }

  ## Random generator
  out$r <- function(n, par){
    stats::rbinom(n, size = 1, prob = par)
  }

  ## Parameters
  out$parameters <- "0 < theta < 1"

  ## Number of parameters
  out$npar <- 1

  ## Name
  out$name <- "Bernoulli"

  ## Constraint
  out$constraint <- NULL

  ## Dispersions
  out$disp <- "Under-dispersion"

  out
}

### BerPoi ---------------------------------------------------------------------
BP <- function(){
  out <- list()

  ## Probability mass function
  out$d <- function(x, par){
    mu <- par[1]
    phi <- par[2]

    lambda <- mu - sqrt(mu * (1 - phi))
    alpha <- sqrt(mu * (1 - phi))
    prob <- (1 - alpha) * stats::dpois(x, lambda) +
      alpha * stats::dpois(x - 1, lambda)
  }

  ## Random generator
  out$r <- function(n, par){
    mu <- par[1]
    phi <- par[2]

    lambda <- mu - sqrt(mu * (1 - phi))
    alpha <- sqrt(mu * (1 - phi))

    stats::rpois(n, lambda) + stats::rbinom(n, size = 1, prob = alpha)
  }

  ## Parameters
  out$par <- "theta > 0; 0 < phi < 1"

  ## Number of parameters
  out$npar <- 2

  ## Name
  out$name <- "BerPoi"

  ## Constraint
  out$constraint <- "phi > 1 - min{theta, 1/theta}"

  ## Dispersions
  out$disp <- "Under-dispersion"

  out
}

### BerG -----------------------------------------------------------------------
BG <- function(){
  out <- list()

  ## Probability mass function
  out$d <- function(x, par){
    mu <- par[1]
    phi <- par[2]

    p0 <- (1 - mu + phi) / (1 + mu + phi)
    p  <- 4 * mu * ((mu + phi - 1)^(x[x > 0] - 1)) / ((mu + phi + 1)^(x[x > 0] + 1))

    prob <- c(rep(p0, sum(x == 0)),p)
    index <- c(which(x == 0), which(x > 0))

    prob[sort(index, index.return = T)$ix]
  }

  ## Random generator
  out$r <- function(n, par){
    mu <- par[1]
    phi <- par[2]

    p <- stats::runif(n)

    p0 <- (1 - mu + phi)/(1 + mu + phi)

    ifelse(length(p) > 1, p.star <- p[p > p0], p.star <- p)
    ifelse(length(mu) > 1, mu.star <- mu[p > p0], mu.star <- mu)
    ifelse(length(phi) > 1, phi.star <- phi[p > p0], phi.star <- phi)

    q <- ceiling(
      round(log((1 - p.star) * (1 + mu.star + phi.star) / (2 * mu.star)) /
              log((mu.star + phi.star - 1) / (mu.star + phi.star + 1)), 2)
    )

    quanti <- c(rep(0, sum(p <= p0)), q)
    index <- c(which(p <= p0), which(p > p0))

    quanti[sort(index, index.return = TRUE)$ix]
  }

  ## Parameters
  out$parameters <- "theta, phi > 0"

  ## Number of parameters
  out$npar <- 2

  ## Name
  out$name <- "BerG"

  ## Constraint
  out$constraint <- "phi > |theta - 1|"

  ## Dispersions
  out$disp <- "Under-, equi-, and over-dispersion"

  out
}

### Mean-Parameterized COM-Poisson ---------------------------------------------
CP <- function(){
  out <- list()

  ## Probability mass function
  out$d <- function(x, par){
    mpcmp::dcomp(x, mu = par[1], nu = 1/par[2])
  }

  ## Random generator
  out$r <- function(n, par){
    mpcmp::rcomp(n, mu = par[1], nu = 1/par[2])
  }

  ## Parameters
  out$parameters <- "theta, phi > 0"

  ## Number of parameters
  out$npar <- 2

  ## Name
  out$name <- "Mean-Parameterized COM-Poisson"

  ## Constraint
  out$constraint <- NULL

  ## Dispersions
  out$disp <- "Under-, equi-, and over-dispersion"

  out

}

### Geometric ------------------------------------------------------------------
GE <- function(){
  out <- list()

  ## Probability mass function
  out$d <- function(x, par){
    stats::dgeom(x, prob = 1/(par + 1))
  }

  ## Random generator
  out$r <- function(n, par){
    stats::rgeom(n, prob =  1/(par + 1))
  }

  ## Parameters
  out$parameters <- "theta > 0"

  ## Number of parameters
  out$npar <- 1

  ## Name
  out$name <- "Geometric"

  ## Constraint
  out$constraint <- NULL

  ## Dispersions
  out$disp <- "Over-dispersion"

  out
}

### Poisson --------------------------------------------------------------------
PO <- function(){
  out <- list()

  ## Probability mass function
  out$d <- function(x, par){
    stats::dpois(x, lambda = par)
  }

  ## Random generator
  out$r <- function(n, par){
    stats::rpois(n, lambda = par)
  }

  ## Parameters
  out$parameters <- "theta > 0"

  ## Number of parameters
  out$npar <- 1

  ## Name
  out$name <- "Poisson"

  ## Constraint
  out$constraint <- NULL

  ## Dispersions
  out$disp <- "Equi-dispersion"

  out
}

