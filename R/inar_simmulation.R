#' @name simulation
#'
#' @title Simulate From an INAR(p) Model
#'
#' @description Simulate from an INAR(p) model.
#'
#' @param n A strictly positive integer given the length of the output series.
#' @param alpha A vector of INAR coefficients.
#' @param par A vector with the innovation process parameters.
#' @param inv Character specification of the innovation process, see
#'     details.
#' @param n.start The length of 'burn-in' period. If \code{NA},
#'      the default, a reasonable value (500) is computed.
#'
#' @return A time-series object of class \code{"ts"}.
#'
#' @details There are some discrete distributions available for the
#'     innovation process specification. The following table display their
#'     names and their abbreviations to be passed to \code{\link[tsinteger]{innovation}()}.
#'
#'   \tabular{llll}{
#'  \bold{Distribution}  \tab \bold{Abbreviation} \tab \tab \bold{Parameters}\cr
#'  Bernoulli \tab \code{"BE"} \tab \tab \code{0 < theta < 1} \cr
#'  BerPoi    \tab \code{"BP"} \tab \tab \code{theta > 0; 0 < phi < 1} \cr
#'  BerG    \tab \code{"BG"} \tab \tab \code{theta, phi > 0} \cr
#'  Geometric \tab \code{"GE"} \tab \tab \code{theta > 0} \cr
#'  Mean-Parameterized COM-Poisson \tab  \code{"CP"} \tab \tab \code{theta, phi > 0} \cr
#'  Negative Binomial \tab \code{"NB"} \tab \tab \code{theta, phi > 0} \cr
#'  Poisson  \tab \code{"PO"}      \tab \tab  \code{theta > 0}  \cr
#'  }
#'
#' @references
#'   Du, J.G. and Li,Y. (1991).
#'   The integer-valued autorregressive (INAR(p)) model.
#'   \emph{Journal of time series analysis}. \bold{12}, 129--142.
#'
#' @examples
#' # A Poisson INAR(3) simulation
#' y <- inar.sim(n = 1000, alpha = c(0.2, 0.3, 0.3), par = 5)
#' layout(matrix(c(1,2,1,2,1,3,1,3), ncol = 4))
#' plot(y)
#' acf(y, main = "")
#' pacf(y, main = "")
#' layout(1)
#'
#' @export
inar.sim <- function(n, alpha, par, inv = "PO", n.start = NA){

  ## Autoregressive order
  p <- length(alpha)

  ## Innovation process
  inv <- innovation(inv)

  ## Length of burn-in period
  if (is.na(n.start)) n.start <- 500

  ## Main body
  nn <- n.start + n

  e <- inv$r(nn, par)
  y <- vector()
  y[1:p] <- e[1:p]

  ## Thinning binomial operation
  thinning <- function(alpha, y){
    stats::rbinom(1, y, alpha)
  }

  # Random generation
  for (t in (p + 1):nn){
    y[t] <- sum(mapply(thinning, alpha, y[t - p:1])) + e[t]
  }

  stats::ts(y[(n.start + 1):nn])
}
