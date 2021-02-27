#' @name innovation
#'
#' @title Family of Innovations Process for Fitting an INAR(p) Model
#'
#' @description Provide the current available distributions that can be
#'     used as an innovation process in a fit of the INAR(p) model.
#'
#' @param dist Character specification of the innovation process, see
#'     details.
#' @param x A \code{"innovation"} object.
#' @param ... Further arguments for other specific methods.
#'
#' @details There are some discrete distributions available for the
#'     innovation process specification. The following table display their
#'     names and their abbreviations to be passed to \code{innovation()}.
#'
#'   \tabular{lll}{
#'  \bold{Distribution}  \tab \bold{Abbreviation} \tab \bold{Number of parameters}\cr
#'  BerPoi    \tab \code{"BP"} \tab 2 \cr
#'  BerG    \tab \code{"BG"} \tab 2 \cr
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
  dist <- match.fun(dist)
  dist <- eval(dist())
  class(dist) <- "innovation"
  dist
}

#' @rdname innovation
#' @export
print.innovation <- function(x, ...){

  innovations <- c("BP", "BG", "GE", "PO")
  INNOVATIONS <- c("BerPoi", "BerG", "Geometric", "Poisson")

  cat("----------------------------------------",
      "\nInnovation process:", x$name,
      "\n----------------------------------------",
      "\nAbreviation:", innovations[INNOVATIONS == x$name],
      "\nNumber of parameters:", x$npar)

      if (x$npar > 1){
       cat("\nConstraints:", x$constraint)
      }
      cat("\nDispersions:", x$disp)
      cat("\n---")
}

#' @rdname innovation
#' @export
plot.innovation <- function(x, ...){

  ### BerPoi distribution ------------------------------------------------------
  if (x$name == "BerPoi"){
    op <- graphics::par(mfrow = c(3, 1), mar = c(2, 10, 2, 10))
    xvals <- 0:15

    ## Varying mu
    graphics::plot(xvals, x$d(xvals, 0.5, 0.9), type = "n",
                   ylab = "Density", xlab = "x")
    mu <- seq(0.5, 5, length.out = 30)
    col <- RColorBrewer::brewer.pal(n = 5, name = "Spectral")[5:1]
    col <- grDevices::colorRampPalette(colors = col)(length(mu))
    for (a in seq_along(mu)) {
      graphics::lines(xvals, x$d(xvals, mu[a], 0.9), col = col[a], type = "b", pch = 16)
    }

    ## Varying phi
    graphics::plot(xvals, x$d(xvals, 1, 0.1), type = "n",
                   ylab = "Density", xlab = "x")
    phi <- seq(0.1, 0.9, length.out = 30)
    col <- RColorBrewer::brewer.pal(n = 5, name = "Spectral")[5:1]
    col <- grDevices::colorRampPalette(colors = col)(length(phi))
    for (a in seq_along(phi)) {
      graphics::lines(xvals, x$d(xvals, 1, phi[a]), col = col[a], type = "b", pch = 16)
    }

    ## Constraints
    f <- function(mu){
      p1 <- mu
      p2 <- 1/mu

      out <- vector()
      out[p1 <= p2] <- 1 - p1[p1 <= p2]
      out[p2 <= p1] <- 1 - p2[p2 <= p1]
      out
    }

    graphics::curve(f(x), xlab = "Mean", ylab = "Dispersion index",
          xlim = c(0, 15), col = "white")
    graphics::polygon(x = c(seq(0, 20, 0.001), seq(20, 0, -0.001)),
            y = c(f(seq(0, 20, 0.001)), rep(1, length(seq(0, 20, 0.001)))),
            col = "darkgray", border = F)
    graphics::abline(h = 1, lty = 2)
    graphics::legend("bottomright", legend = "Poisson",
                     lty = 2, bty = "n")
    graphics::box()
  }

  ### BerG distribution --------------------------------------------------------
  if (x$name == "BerG"){
    op <- graphics::par(mfrow = c(3, 1), mar = c(2, 10, 2, 10))
    xvals <- 0:15

    ## Varying mu
    graphics::plot(xvals, x$d(xvals, 0.5, 6), type = "n",
                   ylab = "Density", xlab = "x")
    mu <- seq(0.5, 4.9, length.out = 30)
    col <- RColorBrewer::brewer.pal(n = 5, name = "Spectral")[5:1]
    col <- grDevices::colorRampPalette(colors = col)(length(mu))
    for (a in seq_along(mu)) {
      graphics::lines(xvals, x$d(xvals, mu[a], 4), col = col[a], type = "b", pch = 16)
    }

    ## Varying phi
    graphics::plot(xvals, x$d(xvals, 1, 0.1), type = "n",
                   ylab = "Density", xlab = "x")
    phi <- seq(0.1, 5, length.out = 30)
    col <- RColorBrewer::brewer.pal(n = 5, name = "Spectral")[5:1]
    col <- grDevices::colorRampPalette(colors = col)(length(phi))
    for (a in seq_along(phi)) {
      graphics::lines(xvals, x$d(xvals, 1, phi[a]), col = col[a], type = "b", pch = 16)
    }

    ## Constraints
    f1 <- function(mu) return(mu-1)
    f2 <- function(mu) return(1 - mu)

    graphics::curve(f1(x), xlab = "Mean", ylab = "Dispersion index",
                    xlim = c(0, 15), ylim = c(0, 15),col="white")
    graphics::polygon(x = c(seq(0, 1, 0.001), seq(1, 20, 0.001), seq(20, 0, -0.001)),
            y = c(f2(seq(0, 1, 0.001)), f1(seq(1, 20, 0.001)),
                  rep(20, length(seq(20, 0, -0.001)))),
            col="gray80",border = F)
    graphics::abline(1, 0, lty = 2)
    graphics::legend("right", legend = "Poisson",
                     lty = 2, bty = "n")

    graphics::box()
  }

  ### Geometric distribution ---------------------------------------------------
  if (x$name == "Geometric"){
    op <- graphics::par(mfrow = c(1, 2), mar = c(5.5, 4, 4.5, 1))
    xvals <- 0:20
    graphics::plot(xvals, x$d(xvals, 0.5), type = "n",
                   ylab = "Density", xlab = "x")
    mu <- seq(0.5, 10, length.out = 30)
    col <- RColorBrewer::brewer.pal(n = 5, name = "Spectral")[5:1]
    col <- grDevices::colorRampPalette(colors = col)(length(mu))
    for (a in seq_along(mu)) {
      graphics::lines(xvals, x$d(xvals, mu[a], phi), col = col[a],
                      type = "b", pch = 16)
    }

    graphics::curve(1 + x, xlab = "Mean", xlim = c(0, 10),
                    ylab = "Dispersion index", lwd = 1.5)
    graphics::abline(1, 0, lty = 2)
    graphics::legend("topleft", legend = c("Geometric", "Poisson"),
                     lty = c(1, 2), bty = "n")
  }

  ### Poisson distribution -----------------------------------------------------
  if (x$name == "Poisson"){
    op <- graphics::par(mfrow = c(1, 2), mar = c(5.5, 4, 4.5, 1))
    xvals <- 0:20
    graphics::plot(xvals, x$d(xvals, 0.5), type = "n",
                   ylab = "Density", xlab = "x")
    mu <- seq(0.5, 10, length.out = 30)
    col <- RColorBrewer::brewer.pal(n = 5, name = "Spectral")[5:1]
    col <- grDevices::colorRampPalette(colors = col)(length(mu))
    for (a in seq_along(mu)) {
      graphics::lines(xvals, x$d(xvals, mu[a], phi), col = col[a],
                      type = "b", pch = 16)
    }

    graphics::plot(0:10, 0:10, type = "n", xlab = "Mean",
                   ylab = "Dispersion index", lwd = 1.5)
    graphics::abline(1, 0)
  }

  graphics::par(op)
}


### Distributions: ----

### BerPoi ---------------------------------------------------------------------
BP <- function(){
  out <- list()

  ## Generating function
  out$d <- function(x, mu, phi){
    lambda <- mu - sqrt(mu * (1 - phi))
    alpha <- sqrt(mu * (1 - phi))
    prob <- (1 - alpha) * stats::dpois(x, lambda) +
      alpha * stats::dpois(x - 1, lambda)
  }

  ## Random generator
  out$r <- function(n, mu, phi){
    lambda <- mu - sqrt(mu * (1 - phi))
    alpha <- sqrt(mu * (1 - phi))

    stats::rpois(n, lambda) + stats::rbinom(n, size = 1, prob = alpha)
  }

  ## Number of parameters
  out$npar <- 2

  ## Name
  out$name <- "BerPoi"

  ## Constraint
  out$constraint <- "phi > 1 - min{mu, 1/mu}"

  ## Dispersions
  out$disp <- "Under-dispersion"

  out
}

### BerG -----------------------------------------------------------------------
BG <- function(){
  out <- list()

  ## Generating function
  out$d <- function(x, mu, phi){
    p0 <- (1 - mu + phi) / (1 + mu + phi)
    p  <- 4 * mu * ((mu + phi - 1)^(x[x > 0] - 1)) / ((mu + phi + 1)^(x[x > 0] + 1))

    prob <- c(rep(p0, sum(x == 0)),p)
    index <- c(which(x == 0), which(x > 0))

    prob[sort(index, index.return = T)$ix]
  }

  ## Random generator
  out$r <- function(n, mu, phi){

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

  ## Number of parameters
  out$npar <- 2

  ## Name
  out$name <- "BerG"

  ## Constraint
  out$constraint <- "phi > |mu - 1|"

  ## Dispersions
  out$disp <- "Under-, equi-, and over-dispersion"

  out
}


### Geometric ------------------------------------------------------------------
GE <- function(){
  out <- list()

  ## Generating function
  out$d <- function(x, mu, phi = NULL){
    stats::dgeom(x, prob = 1/(mu + 1))
  }

  ## Random generator
  out$r <- function(n, mu, phi = NULL){
    stats::rgeom(n, prob =  1/(mu + 1))
  }

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

  ## Generating function
  out$d <- function(x, mu, phi = NULL){
    stats::dpois(x, lambda = mu)
  }

  ## Random generator
  out$r <- function(n, mu, phi = NULL){
    stats::rpois(n, lambda = mu)
  }

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

