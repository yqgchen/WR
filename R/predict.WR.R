#' Predict responses for a new sample given a WR object
#' @param object A WR object, i.e., an output of \code{\link{WR}}.
#' @param Xpred A list holding the values of distributional and/or scalar predictors 
#' for a held-out set awaiting prediction. The form and structure are similar to 
#' the input \code{X} of \code{\link{WR}}, except that the number of realizations 
#' awaiting prediction can be different from the number of realizations \eqn{n} used to fit the model.
#' @param ... Not used.
#' @return A list of the following
#' \item{Ypred}{Predicted responses; either a vector of length \eqn{n} for scalar responses 
#' or a \eqn{q}-by-\eqn{n} matrix holding the quantile functions for distributional responses, 
#' evaluated on \code{qSup}, where \eqn{q} is the length of \code{qSup}.}
#' \item{qSup}{The input \code{qSup}, only returned when the responses are distributional.}
#' \item{outOfLogSpace}{A logical vector holding whether initial predicted log mapped distributions 
#' lie out of the log image space (\code{TRUE}) or not (\code{FALSE}), only returned when the responses are distributional.}
#' @import stats
#' @method predict WR
#' @export
#' 
predict.WR <- function ( object, Xpred, ... ) {
  
  beta <- object$beta
  workGridY <- object$workGridY
  fpcaLogX <- object$fpcaLogX
  fpcaLogY <- object$fpcaLogY
  X <- object$X
  Y <- object$Y
  qSup <- object$qSup
  optns <- object$optns
  
  checkXY <- CheckData4WR( X = X, Y = Y, qSup = qSup )
  
  scalarRps <- checkXY$scalarRps
  anyScalarPdt <- checkXY$anyScalarPdt
  anyDistnlPdt <- checkXY$anyDistnlPdt
  indScalarPdt <- checkXY$indScalarPdt
  indDistnlPdt <- checkXY$indDistnlPdt
  n <- checkXY$n
  q <- checkXY$q
  
  nPred <- CheckXpred4WR( Xpred, resCheckData = checkXY )
  
  # obtain input data for prediction: scalar Xpred and Log(dist'nal Xpred) ----
  inputXpred <- list()
  if ( anyScalarPdt ) {
    inputXpred[ indScalarPdt ] <- lapply( indScalarPdt, function (j) {
      Xpred[[j]] - mean( X[[j]] )
    })
  }
  if ( anyDistnlPdt ) {
    QXmean <- lapply( X[ indDistnlPdt ], rowMeans )
    LogXpred <- lapply( indDistnlPdt, function (j) { Xpred[[j]] - QXmean[[j]] } )
    
    inputXpred[ indDistnlPdt ] <- 
      lapply( seq_along( indDistnlPdt ), function (j) {
        if ( all.equal( fpcaLogX[[j]]$workGrid, qSup ) ) {
          t( LogXpred[[j]] )
        } else {
          t( apply( LogXpred[[j]], 2, function (y) {
            approx( x = qSup, y = y, xout = fpcaLogX[[j]]$workGrid )$y
          }) )
          # predict(
          #   object = fpcaLogX[[j]],
          #   newLy = lapply( seq_len(nPred), function (i) LogXpred[[j]][,i] ),
          #   newLt = rep( list(qSup), nPred ),
          # )$predCurves # evaluated on fpcaLogX[[j]]$workGrid
        }
      })
  }
  
  # obtain predicted responses ----
  Ypred <- list()
  # Ypred: initially a list holding the fit/prediction relevant to each predictor. 
  if ( scalarRps ) {
    ## for scalar responses ----
    if ( anyScalarPdt ) {
      Ypred[ indScalarPdt ] <- lapply( indScalarPdt, function (j) {
        inputXpred[[j]] * beta[[j]]
      })
    }
    if ( anyDistnlPdt ) {
      Ypred[ indDistnlPdt ] <- lapply( indDistnlPdt, function (j) {
        as.vector( inputXpred[[j]] %*% ( beta[[j]] * GetIntgWts( fpcaLogX[[j]]$workGrid ) ) )
      })
    }
  } else {
    ## for distributional responses ----
    if ( anyScalarPdt ) {
      Ypred[ indScalarPdt ] <- lapply( indScalarPdt, function (j) {
        matrix( beta[[j]], ncol = 1 ) %*% matrix( inputXpred[[j]], nrow = 1 )
      })
    }
    if ( anyDistnlPdt ) {
      Ypred[ indDistnlPdt ] <- lapply( indDistnlPdt, function (j) {
        t( inputXpred[[j]] %*% ( beta[[j]] * GetIntgWts( fpcaLogX[[j]]$workGrid ) ) )
      })
    }
  }
  # combine elements of Ypred to obtain predicted scalar responses or log distributional responses
  Ypred <- Reduce( '+', Ypred )
  
  if ( scalarRps ) {
    ## add back mean of Y for scalar responses ----
    Ypred <- Ypred + mean(Y)
  } else {
    ## map fitted log distributional responses back to distributions (quantile functions) ----
    LogYSup <- workGridY
    LogYpred <- Ypred
    Ypred <- GetProj( LogYhat = LogYpred, LogSup = LogYSup )
    outOfLogSpace <- Ypred$outOfLogSpace
    Ypred <- Ypred$Yhat # length(fpcaLogY$workGrid) x n
    
    ## obtain fitted quantile functions evaluated on qSup for distributional responses ----
    if ( !isTRUE( all.equal(fpcaLogY$workGrid, qSup) ) ) {
      Ypred <- apply( Ypred, 2, function (y) {
        approx( x = fpcaLogY$workGrid, y = y, xout = qSup )$y
      })
    }
  }
  
  # summary output ----
  res <- list(
    Ypred = Ypred
  )
  if ( !scalarRps ) {
    res$qSup <- qSup
    res$outOfLogSpace <- outOfLogSpace
  }
  return (res)
}
