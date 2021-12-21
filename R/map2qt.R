#' Projecting functions to the space of non-decreasing functions.
#' @param y A function to be projected.
#' @param lower,upper Lower and upper bounds of the image space of the function. 
#' Defaults are both \code{NULL}, corresponding to an unbounded image space, i.e., the entire real line.
#' @return The projected function.
#' @references
#' \cite{Petersen, A., & Müller, H.-G. (2019). "Fréchet regression for random objects with Euclidean predictors." The Annals of Statistics, 47(2), 691--719.}
#' @export
map2qt <- function ( y, lower = NULL, upper = NULL ) {
  m <- length(y)
  A <- cbind(diag(m), rep(0,m)) + cbind(rep(0,m), -diag(m))
  if (!is.null(upper) & !is.null(lower)) {
    b0 <- c(lower, rep(0,m-1), -upper)
  } else if(!is.null(upper)) {
    A <- A[,-1]
    b0 <- c(rep(0,m-1), -upper)
  } else if(!is.null(lower)) {
    A <- A[,-ncol(A)]
    b0 <- c(lower,rep(0,m-1))
  } else {
    A <- A[,-c(1,ncol(A))]
    b0 <- rep(0,m-1)
  }
  Pmat <- as(diag(m), "sparseMatrix")
  Amat <- as(t(A), "sparseMatrix")
  res <- do.call(
    osqp::solve_osqp,
    list(P=Pmat, q= -y, A=Amat, l=b0, pars = osqp::osqpSettings(verbose = FALSE))
  )
  sort(res$x)
}
