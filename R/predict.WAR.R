#' Predict distributions given a WAR object
#' @param object A WAR object, i.e., an output of \code{\link{WAR}}.
#' @param n.ahead Number of steps ahead at which to predict. Default: 1.
#' @param Ynew A \eqn{q}-by-\eqn{m} matrix holding the quantile functions of 
#' \eqn{m} distributions evaluated on \code{object$qSup} to which to apply the prediction, 
#' where \eqn{m} must be greater than or equal to \code{object$order}.
#' Default: \code{object$Y}.
#' @param ... Not used.
#' @return A list of the following elements:
#' \item{Ypred}{A \eqn{q}-by-\eqn{l} matrix holding the predicted quantile functions of 
#' \eqn{l} distributions ahead of \code{Ynew}, where \eqn{l} is given by \code{n.ahead}.}
#' \item{qSup}{The input \code{qSup}}
#' \item{outOfLogSpace}{A logical vector holding whether initial predicted log mapped distributions 
#' lie out of the log image space (\code{TRUE}) or not (\code{FALSE}).}
#' @import stats
#' @method predict WAR
#' @export
#' 
predict.WAR <- function ( object, n.ahead = 1, Ynew = NULL, ... ) {
  arFPCs <- object$arFPCs
  fpcaLogY <- object$fpcaLogY
  qSup <- object$qSup
  Y <- object$Y
  
  nFPCs <- length( arFPCs$x.mean )
  QYmean <- rowMeans(Y)
  
  # obtain FPCs time series to which prediction is to be applied ----
  if ( is.null(Ynew) ) {
    fpcsYnew <- fpcaLogY$xiEst[,1:nFPCs,drop = FALSE]
  } else {
    CheckData4WR( X = NULL, Y = Ynew, qSup = qSup )
    nPred <- ncol(Ynew)
    LogYnew <- Ynew - QYmean
    fpcsYnew <- predict(
      fpcaLogY,
      newLy = lapply( seq_len(nPred), function(i) LogYnew[,i] ),
      newLt = rep( list(qSup), nPred ),
      K = nFPCs
    )$scores
  }
  fpcsYnew <- ts( data = fpcsYnew, frequency = 1, start = 1 )
  
  # obtain predicted FPCs ----
  fpcsYpred <- predict( arFPCs, newdata = fpcsYnew, n.ahead = n.ahead, se.fit = FALSE )
  # obtain predicted log maps and their support grid ----
  LogYpred <- fpcaLogY$phi[,1:nFPCs, drop = FALSE] %*% t(fpcsYpred) + fpcaLogY$mu
  if ( all.equal( qSup, fpcaLogY$workGrid ) ) {
    LogSup <- QYmean
  } else {
    LogSup <- approx( x = qSup, y = QYmean, xout = fpcaLogY$workGrid )$y
  }
  # map predicted log maps back to distributions ----
  Ypred <- GetProj( LogYhat = LogYpred, LogSup = LogSup )
  outOfLogSpace <- Ypred$outOfLogSpace
  Ypred <- Ypred$Yhat
  # obtain predicted quantile functions on qSup ----
  if ( !all.equal(fpcaLogY$workGrid, qSup) ) {
    Ypred <- apply( Ypred, 2, function (y) {
      approx( x = fpcaLogY$workGrid, y = y, xout = qSup )$y
    })
  }
  
  # summary output ----
  return ( list(
    Ypred = Ypred,
    qSup = qSup,
    outOfLogSpace = outOfLogSpace
  ))
}
