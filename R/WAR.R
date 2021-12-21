#' Wasserstein autoregressive models.
#' @param Q A \eqn{q}-by-\eqn{n} matrix holding the quantile functions of 
#' \eqn{n} distributions evaluated on a common grid on \eqn{[0,1]} consisting of \eqn{q} points.
#' @param qSup A vector of length \eqn{q} holding the support grid of quantile functions.
#' @param optns A list of options control parameters specified by list(name=value). See 'Details'.
#' @param n.ahead Number of steps ahead at which to predict. Default: 1.
#' @param aic,order.max,method,... arguments for \code{\link[stats]{ar}}. 
#' Default: \code{aic = FALSE, order.max = 1, method = "yule-walker"}. 
#' Note that \code{demean} should not be input and is set to be \code{FALSE}.
#' @details Available control options are
#' \describe{
#' \item{methodNFPC}{The method for choosing the number of FPCs; \code{'FVE'} (default), \code{'CV'}.}
#' \item{methodProj}{The projection method if the predictions of log mapped distributional responses lie out of the log image space: 
#' \code{'log'} (the method as per Chen, Lin, and Müller, 2021); 
#' \code{'qt'} (default, the method as per Petersen & Müller, 2019).}
#' \item{FVEthreshold}{Fraction-of-Variance-Explained threshold used for choosing the number of FPCs, 
#' relevant only if \code{methodNFPC} is \code{'FVE'}; numeric \eqn{(0,1]}, default: 0.95.}
#' \item{CVfold}{The number of folds for cross validation, relevant only if \code{methodNFPC} is \code{'CV'}; 
#' Any positive integer up to \eqn{n}; Default: \eqn{n} (leave-one-out) if \eqn{n\le 30} and \eqn{10} (10-fold) if \eqn{n > 30}.}
#' \item{lower}{Lower bound of the support of distributions; default: \code{NULL}, i.e., no finite lower bound.}
#' \item{upper}{Upper bound of the support of distributions; default: \code{NULL}, i.e., no finite upper bound.}
#' }
#' @return A list of class \code{'WAR'} with the following elements:
#' \item{QPred}{A \eqn{q}-by-\code{n.ahead} matrix holding the quantile functions of the \code{n.ahead} predicted distributions.}
#' \item{qSup}{The support of quantile functions, same as the input.}
#' \item{order}{The order of the fitted Wasserstein autoregressive model.}
#' \item{optns}{The \code{optns} used.}
#' @references
#' \cite{Chen, Y., Lin, Z., & Müller, H.-G. (2021). "Wasserstein regression." Journal of the American Statistical Association, in press.}
#' \cite{Petersen, A., & Müller, H.-G. (2019). "Fréchet regression for random objects with Euclidean predictors." The Annals of Statistics, 47(2), 691--719.}
#' @export
#' 
WAR <- function(Q, qSup, optns = list(), n.ahead = 1, aic = FALSE, order.max = 1, method = "yule-walker", ...){
  n <- ncol(Q)
  q <- nrow(Q)
  if ( length(qSup) != q ) {
    stop('Length of qSup does not match the number of rows in Q.')
  }
  if ( is.null(optns$methodNFPC) ) {
    optns$methodNFPC <- 'FVE'
  } else {
    if ( !optns$methodNFPC %in% c('FVE','CV') ) {
      stop ('Unrecognized value for optns$methodNFPC.')
    }
  }
  if ( optns$methodNFPC == 'FVE' ) {
    if ( is.null(optns$FVEthreshold) ) {
      optns$FVEthreshold <- 0.95
    } else {
      if ( optns$FVEthreshold > 1 | optns$FVEthreshold <= 0 ) {
        stop ('Wrong value for optns$FVEthreshold.')
      }
    }
  } else if ( optns$methodNFPC == 'CV' ) {
    if ( is.null(optns$CVfold) ) {
      optns$CVfold <- ifelse( n > 30, 10, n )
    } else {
      if ( optns$CVfold > n | optns$CVfold < 0 ) {
        stop ('Misspecified optns$CVfold.')
      }
    }
  }
  
  if ( is.null(optns$methodProj) ) {
    optns$methodProj <- 'log'
  } else {
    if ( !optns$methodProj %in% c('qt',"log") ) {
      stop('Unrecognized value for optns$methodProj.')
    }
  }
  
  qf <- Q
  qfMean <- rowMeans(qf)
  Log <- qf - qfMean
  colnames(Log) <- 1:ncol(Log)
  
  fpcaMan <- FPCA(
    Ly = as.list(data.frame(Log)),
    Lt = lapply(1:ncol(Log), function(i) qSup),
    optns = list(useBinnedData='OFF', FVEthreshold = 0.9999)
  )
  if ( optns$methodNFPC == 'FVE' ) {
    nFPCs_man <- min( which( fpcaMan$cumFVE >= optns$FVEthreshold ) )
  } else if ( optns$methodNFPC == 'CV' ) {
    cv_partition <- sample(seq_len(n-1), n-1)
    cv_partition <- findInterval( x = seq_len(n), seq(1,n,length.out = nFPCs_man + 1)[1:nFPCs_man] )
    
    MSPE_cv_man <- foreach(nFPCs = 1:length(fpcaMan$cumFVE), .combine = c) %do% {
      require(MASS)
      X <- fpcaMan$xiEst[,1:nFPCs,drop = FALSE]
      Cov1 <- lapply(2:nrow(X), function(t) {
        X[t,] %*% t(X[t-1,])
      })
      Cov0 <- lapply(1:(nrow(X)-1), function(t) {
        X[t,] %*% t(X[t,])
      })
      Log.pred <- lapply(1:length(cv_partition), function(k) {
        phi1 <- Reduce("+", Cov1[-cv_partition[[k]]]) %*% ginv(Reduce("+", Cov0[-cv_partition[[k]]]))
        X.pred <- t(phi1 %*% t(X[cv_partition[[k]],drop = FALSE]))
        return(fpcaMan$phi[,1:nFPCs,drop = FALSE] %*% t(X.pred) + fpcaMan$mu)
      })
      Log.pred <- Reduce(cbind,Log.pred)
      Log.pred <- Log.pred[, order(unlist(cv_partition))]
      Log.pred <- sapply(1:ncol(Log.pred), function(i) {
        approx(x = fpcaMan$workGrid, y = Log.pred[,i], xout = qSup, rule = 2)$y
      })
      return(mean(apply((Log.pred - Log[,-1])^2, 2, pracma::trapz, x = qSup)))
    }
    nFPCs_man <- which.min(MSPE_cv_man)
  }
  
  ## treated as multivariate time series
  fpcManTS <- ts(data = fpcaMan$xiEst[,1:nFPCs_man, drop = FALSE], frequency = 1, start = 1)
  reMan <- ar(fpcManTS, aic = aic, order.max = order.max, method = method, demean = FALSE, ...)
  predMan <- predict(reMan, newdata = fpcManTS, n.ahead = n.ahead, se.fit = FALSE)
  
  # predicted Log functions on fpcaMan$workGrid
  LogPred <- fpcaMan$phi[,1:nFPCs_man, drop = FALSE] %*% t(predMan) + fpcaMan$mu
  LogSup <- approx( x = qSup, y = qfMean, xout = fpcaMan$workGrid )$y
  
  QPred <- LogPred + LogSup
  outOfLogSpace <- apply( QPred, 2, is.unsorted )
  # pull back after one-on-one prediction
  if ( any( outOfLogSpace ) ) {
    if ( optns$methodProj == 'qt' ) {
      QPred[,which(outOfLogSpace)] <- apply( QPred[,which(outOfLogSpace),drop=FALSE], 2, map2qt, lower = optns$lower, upper = optns$upper )
    } else if ( optns$methodProj == 'log' ) {
      LogPred.bo <- LogPred %*% diag( apply( LogPred, 2, map2boundary, sup = LogSup ) )
      QPred <- LogPred.bo + LogSup
    }
  }
  
  QPred <- apply( QPred, 2, function(y) {
    approx( x = fpcaMan$workGrid, y = y, xout = qSup )$y
  })
  
  res <- list(
    QPred = QPred,
    qSup = qSup,
    order = reMan$order,
    optns = optns
  )
  class(res) <- "WAR"
  return (res)
}
