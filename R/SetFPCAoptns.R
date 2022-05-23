#' Set FPCA option list
#' @param FPCAoptns A list of optional control parameters specified by \code{list(name=value)} for \code{\link[fdapace]{FPCA}}.
#' @return A list of optional control parameters.
#' @export
#' 
SetFPCAoptns <- function ( FPCAoptns = list() ) {
  if ( is.null(FPCAoptns$methodSelectK) ) {
    FPCAoptns$methodSelectK <- 'FVE'
  }
  if ( FPCAoptns$methodSelectK == 'FVE' ) {
    if ( is.null(FPCAoptns$FVEthreshold) ) {
      FPCAoptns$FVEthreshold <- 0.95
    } else {
      if ( FPCAoptns$FVEthreshold > 1 | FPCAoptns$FVEthreshold <= 0 ) {
        stop ('Wrong values for FVEthreshold in FPCA options.')
      }
    }
  }
  if ( is.null(FPCAoptns$useBinnedData) ) {
    FPCAoptns$useBinnedData <- 'OFF'
  }
  # set userMu to be constantly zero.
  FPCAoptns$userMu <- list(
    t = seq( 0, 1, length.out = 51 ),
    mu = rep( 0, 51 )
  )
  
  return (FPCAoptns)
}
