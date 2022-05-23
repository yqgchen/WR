#' Check prediction set format for Wasserstein regression
#' @description Check if there are problems with the form and structure of predictors in the prediction set. 
#' @param Xpred A list holding the values of distributional and/or scalar predictors for a held-out set awaiting prediction. 
#' @param resCheckData A list holding the output of \code{\link{CheckData4WR}}.
#' @return The number of predictions to be made.
#' @export
#' 
CheckXpred4WR <- function ( Xpred, resCheckData ) {
  
  anyScalarPdt <- resCheckData$anyScalarPdt
  anyDistnlPdt <- resCheckData$anyDistnlPdt
  indScalarPdt <- resCheckData$indScalarPdt
  indDistnlPdt <- resCheckData$indDistnlPdt
  n <- resCheckData$n
  q <- resCheckData$q
  
  if ( !is.list(Xpred) ) stop ("Xpred is not a list.")
  scalarPdtPred <- sapply( Xpred, is.vector )
  distnlPdtPred <- sapply( Xpred, is.matrix )
  
  anyScalarPdtPred <- any(scalarPdtPred)
  anyDistnlPdtPred <- any(distnlPdtPred)
  indScalarPdtPred <- which(scalarPdtPred)
  indDistnlPdtPred <- which(distnlPdtPred)
  
  # check the formats of each predictor in Xpred ----
  if ( sum( abs( indScalarPdt - indScalarPdtPred ) ) + sum( abs( indDistnlPdt - indDistnlPdtPred ) ) > 0 ) {
    stop ("Formats of elements in Xpred do not match those in X.")
  }
  
  # check the number of realizations ----
  nPred <- numeric( length(Xpred) )
  if ( anyScalarPdt ) {
    nPred[ indScalarPdtPred ] <- sapply( Xpred[ indScalarPdtPred ], length )
  }
  if ( anyDistnlPdt ) {
    nPred[ indDistnlPdtPred ] <- sapply( Xpred[ indDistnlPdtPred ], ncol )
  }
  if ( sum( abs( diff( nPred ) ) ) > 0 ) stop ("Numbers of realizations of predictors in Xpred are not the same.")
  nPred <- nPred[1]
  
  # check if quantile functions are non-decreasing ----
  if ( anyDistnlPdt ) {
    for ( j in indDistnlPdtPred ) {
      if ( any( apply( Xpred[[j]], 2, diff ) < 0 ) ) {
        stop ("There exists a column of X[[",j,"]] that is not non-descreasing.")
      }
    }
  }
  
  # check the support grid of quantile functions ----
  if ( anyDistnlPdtPred ) {
    qPred <- sapply( Xpred[ indDistnlPdtPred ], nrow )
    if ( any( abs( qPred - q ) > 0 ) ) stop ("Numbers of support grid points of the quantile functions are not the same as the length of qSup for all distributional variables in Xpred.")
  }
  return ( nPred )
}
