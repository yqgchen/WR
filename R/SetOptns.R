#' Set option list for WR and WAR
#' @param optns A list of optional control parameters specified by \code{list(name=value)} for \code{\link{WR}} or \code{\link{WAR}}.
#' @param Ymin,Ymax The minimum and maximum of the input \code{Y} for \code{\link{WAR}} or \code{\link{WR}}, 
#' the latter if the responses are scalar. 
#' @return A list of optional control parameters.
#' @export
#' 
SetOptns <- function ( optns = list(), Ymin = NULL, Ymax = NULL ) {
  if ( is.null(optns$methodProj) ) {
    optns$methodProj <- 'log'
  } else {
    if ( !optns$methodProj %in% c('qt','log') ) {
      stop('Unrecognized value for optns$methodProj.')
    }
  }
  if ( !is.null(Ymin) ) {
    if ( !is.null(optns$lower) ) {
      if ( optns$lower > Ymin ) stop ("Some entries of Y is less than optns$lower.")
    }
  }
  if ( !is.null(Ymax) ) {
    if ( !is.null(optns$upper) ) {
      if ( optns$upper < Ymax ) stop ("Some entries of Y is greater than optns$upper.")
    }
  }
  return (optns)
}
