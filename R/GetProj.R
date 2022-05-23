#' Project to the log image space
#' @description Obtain projections of fitted/predicted log mapped distributions to the log image space if they are not in the space.
#' @param LogYhat A \eqn{p}-by-\eqn{n} matrix holding \eqn{n} fitted/predicted log mapped distributions 
#' evaluated on a common grid given by \code{LogSup} of length \eqn{p}.
#' @param LogSup A vector of length \eqn{p} holding the common support grid.
#' @param optns A list of control parameters specified by \code{list(name=value)}. See \code{\link{WR}} for details.
#' @return A list of the following elements:
#' \item{Yhat}{A \eqn{p}-by\eqn{n} matrix holding the quantile functions of fitted/predicted distributions.}
#' \item{outOfLogSpace}{A logical vector holding whether fitted log mapped distributions 
#' lie out of the log image space (\code{TRUE}) or not (\code{FALSE}).}
#' @export
#' 
GetProj <- function ( LogYhat, LogSup, optns ) {
  Yhat <- LogYhat + LogSup
  outOfLogSpace <- apply( Yhat, 2, is.unsorted )
  if ( any( outOfLogSpace ) ) {
    if ( optns$methodProj == 'qt' ) {
      Yhat[,which(outOfLogSpace)] <- apply(
        Yhat[,which(outOfLogSpace),drop=FALSE], 2, map2qt,
        lower = optns$lower, upper = optns$upper
      )
    } else if ( optns$methodProj == 'log' ) {
      LogYhat.bo <- LogYhat %*% diag( apply( LogYhat, 2, map2boundary, sup = LogSup ) )
      Yhat <- LogYhat.bo + LogSup
    }
  }
  return(list(
    Yhat = Yhat,
    outOfLogSpace = outOfLogSpace
  ))
}
