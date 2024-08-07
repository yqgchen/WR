#' Boundary projection for the prediction of log mapped distributions
#' @param Logfit A vector holding a log map to be projected.
#' @param sup A vector holding the support grid of \code{Logfit}.
#' @return The projection parameter.
#' @references
#' \cite{Chen, Y., Lin, Z., & Müller, H.-G. (2021). "Wasserstein regression." Journal of the American Statistical Association, in press.}
#' @export
map2boundary <- function(Logfit, sup) {
  if (!is.unsorted(Logfit + sup))
    return (1)
  eta0 <- 1
  eta1 <- 0.5
  while (abs(eta0 - eta1) > 1e-6 | is.unsorted(eta1 * Logfit + sup)) {
    if (is.unsorted(eta1 * Logfit + sup)) {
      tmp <- eta1
      eta1 <- eta1 - (eta0 - eta1) / 2
      eta0 <- tmp
    } else {
      eta1 <- (eta0 + eta1) / 2
    }
    #print(c(eta0,eta1))
  }
  return (eta1)
}