#' Get the weights used for numerical integration given a grid
#' @param x A vector holding the grid points where numerical integration is to be taken, the elements of which are in an increasing order.
#' @return A vector holding the weights.
#' @export
#' 
GetIntgWts <- function ( x ) {
  x <- unique(x)
  if ( is.unsorted (x) ) x <- sort(x)
  len <- length(x)
  mids <- ( x[-1] + x[-len] ) / 2
  
  c( mids, x[len] ) - c( x[1], mids )
}
