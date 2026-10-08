#' nmeweibull: The New Modified Exponential-Weibull Distribution
#'
#' Tools for modelling positive survival-time data with the New Modified
#' Exponential-Weibull (NME-Weibull) distribution of Alshanbari et al.
#' (2023). The CDF is
#' \deqn{F(x) = 1 - \frac{\alpha^2 e^{-\lambda x^\beta}}{(\alpha + 1 -
#'   e^{-\lambda x^\beta})^2}, \quad x \ge 0,}
#' with \eqn{\alpha, \beta, \lambda > 0}.
#'
#' The package provides
#' \itemize{
#'   \item distribution functions: [dnmeweibull()], [pnmeweibull()],
#'     [snmeweibull()], [hnmeweibull()], [qnmeweibull()], and
#'     [rnmeweibull()];
#'   \item maximum likelihood estimation for complete data:
#'     [fitnmeweibull()];
#'   \item quantile-based risk measures: [varnmeweibull()],
#'     [tvarnmeweibull()], and [rvarnmeweibull()].
#' }
#'
#' @references Alshanbari, H. M., Odhah, O. H., Ahmad, Z., Khan, F., and
#'   El-Bagoury, A. A.-A. H. (2023). A new probability distribution: Model,
#'   theory and analyzing the recovery time data. *Axioms*, 12(5), 477.
#'   \doi{10.3390/axioms12050477}
#'
#' @keywords internal
"_PACKAGE"
