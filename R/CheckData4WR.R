#' Check data format for Wasserstein regression
#' @description Check if there are problems with the form and structure of input data for \code{\link{WR}}. 
#' @param X,Y The input data of predictors and responses for \code{\link{WR}}, if both are specified;  
#' Otherwise, if \code{X} is not specified and hence is \code{NULL}, \code{Y} is the input \code{Y} of \code{\link{WAR}}.
#' @param qSup Support grid of quantile functions for distributional variables; see \code{\link{WR}} for details. 
#' @return A list of the follow elements if \code{X} is specified:
#' \item{scalarRps}{Logical; whether the responses in \code{Y} are scalars (\code{TRUE}) or distributions (\code{FALSE}).}
#' \item{anyScalarPdt}{Logical; \code{TRUE} if there exists scalar predictors and \code{FALSE} otherwise.}
#' \item{anyDistnlPdt}{Logical; \code{TRUE} if there exists distributional predictors and \code{FALSE} otherwise.}
#' \item{indScalarPdt}{A vector holding the indices of elements in \code{X} corresponding to scalar predictors.}
#' \item{indDistnlPdt}{A vector holding the indices of elements in \code{X} corresponding to distributional predictors.}
#' \item{n}{Number of realizations.}
#' \item{q}{Number of grid points in \code{qSup}.}
#' @export
#' 
CheckData4WR <- function ( X = NULL, Y, qSup ) {
  # check the format of X and Y ----
  if ( !is.vector(Y) & !is.matrix(Y) ) stop ("Y is not a vector or a matrix.")
  scalarRps <- is.vector(Y)
  if ( !is.null(X) ) {
    if ( !is.list(X) ) stop ("X is not a list.")
    scalarPdt <- sapply( X, is.vector )
    distnlPdt <- sapply( X, is.matrix )
    if ( any( !( scalarPdt | distnlPdt ) ) ) stop ("Some elements of X are neither a vector nor a matrix.")
    if ( all(scalarPdt) & scalarRps ) stop ("All the predictors and responses are scalar. Please use lm() instead.")
    anyScalarPdt <- any(scalarPdt)
    anyDistnlPdt <- any(distnlPdt)
    indScalarPdt <- which(scalarPdt)
    indDistnlPdt <- which(distnlPdt)
  }
  
  # check the number of realizations ----
  n <- ifelse( scalarRps, length(Y), ncol(Y) )
  if ( !is.null(X) ) {
    nX <- numeric( length(X) )
    if ( anyScalarPdt ) {
      nX[ indScalarPdt ] <- sapply( X[ indScalarPdt ], length )
    }
    if ( anyDistnlPdt ) {
      nX[ indDistnlPdt ] <- sapply( X[ indDistnlPdt ], ncol )
    }
    if ( sum( abs( diff( nX ) ) ) > 0 ) stop ("Numbers of realizations of predictors in X are not the same.")
    if ( sum( abs( nX[1] - n ) ) > 0 ) stop ("Numbers of realizations of predictors in X and responses in Y are not the same.")
  }
  
  # check if quantile functions are non-decreasing ----
  if ( !is.null(X) ) {
    if ( anyDistnlPdt ) {
      for ( j in indDistnlPdt ) {
        if ( any( apply( X[[j]], 2, diff ) < 0 ) ) {
          stop ("There exists a column of X[[",j,"]] that is not non-descreasing.")
        }
      }
    }
  }
  if ( !scalarRps ) {
    if ( any( apply( Y, 2, diff ) < 0 ) ) {
      stop ("There exists a column of Y that is not non-descreasing.")
    }
  }
  
  # check the support grid of quantile functions ----
  if ( !is.vector(qSup) ) stop ("qSup is not a vector.")
  if ( any( diff(qSup) <= 0 ) ) stop ("qSup is not in a strictly increasing order.")
  if ( min(qSup) < 0 ) stop ("qSup has minimum less than 0.")
  if ( max(qSup) > 1 ) stop ("qSup has maximum greater than 1.")
  q <- length( qSup )
  if ( !scalarRps ) q <- c( q, nrow(Y) )
  if ( !is.null(X) ) {
    if ( anyDistnlPdt ) q <- c( q, sapply( X[ indDistnlPdt ], nrow ) )
  }
  if ( any( abs( q[-1] - q[1] ) > 0 ) ) stop ("Numbers of support grid points of the quantile functions are not the same as the length of qSup for all distributional variables.")
  q <- q[1]
  
  if ( !is.null(X) ) {
    return ( list(
      scalarRps = scalarRps,
      anyScalarPdt = anyScalarPdt, anyDistnlPdt = anyDistnlPdt,
      indScalarPdt = indScalarPdt, indDistnlPdt = indDistnlPdt,
      n = n, q = q
    ) )
  }
}
