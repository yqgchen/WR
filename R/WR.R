#' Wasserstein regression
#' @description Wasserstein regression with 
#' (1) distributional and/or scalar predictors and distributional responses, 
#' (2) distributional (and scalar) predictors and scalar responses. 
#' @param X A list whose elements consist of distributional predictors and/or scalar predictors. 
#' Each distributional predictor has to be represented by a \eqn{q}-by-\eqn{n} matrix 
#' whose columns consist of quantile functions of \eqn{n} realizations of this distributional predictor, 
#' which are evaluated on a common grid \code{qSup} on \eqn{[0,1]} consisting of \eqn{q} points. 
#' Each scalar predictor (if any) has to be represented by a vector of length \eqn{n}. 
#' @param Y Distributional or scalar responses. 
#' A distributional response has to be represented by a \eqn{q}-by-\eqn{n} matrix 
#' whose columns consist of quantile functions of \eqn{n} realizations, 
#' which are evaluated on a common grid \code{qSup} on \eqn{[0,1]} consisting of \eqn{q} points. 
#' A scalar response has to be represented by a vector of length \eqn{n}. 
#' @param qSup A vector of length \eqn{q} holding the common support grid of quantile functions, 
#' of which the elements are in a strictly increasing order and lie between 0 and 1.
#' @param optns A list of optional control parameters for regression specified by \code{list(name=value)}. See 'Details'.
#' @param FPCAoptnsX,FPCAoptnsY Lists of optional control parameters for functional principal 
#' component analysis (FPCA) of log mapped distributional predictors and responses, respectively. 
#' See 'Details' of \code{\link[fdapace]{FPCA}}. 
#' For example, one can specify the methods of choosing numbers of functional principal components (FPCs) 
#' through \code{methodSelectK} and set up Fraction-of-Variance-Explained (FVE) thresholds 
#' through \code{FVEthreshold} when using FVE to choose numbers of FPCs. 
#' Default of \code{methodSelectK}: \code{'FVE'}; 
#' Default of \code{FVEthreshold}: 0.95 (different from \code{\link[fdapace]{FPCA}}).
#' Default of \code{useBinnedData}: \code{'OFF'} (different from \code{\link[fdapace]{FPCA}}).
#' @details Available control options are
#' \describe{
#' \item{methodProj}{
#' The projection method if the fitted values of log mapped distributional responses lie out of the log image space: 
#' \code{'log'} (the method as per Chen, Lin, and Müller, 2021); 
#' \code{'qt'} (default, the method as per Petersen & Müller, 2019).
#' }
#' \item{lower}{A scalar specifying the lower bound of the support of responses if 
#' they are distributional, only relevant if \code{methodProj} is \code{'qt'}. 
#' Default: \code{NULL}, i.e., no finite lower bound.}
#' \item{upper}{A scalar specifying the upper bound of of the support of responses if 
#' they are distributional, only relevant if \code{methodProj} is \code{'qt'}. 
#' Default: \code{NULL}, i.e., no finite upper bound.}
#' }
#' @return A list of the following:
#' \item{Yfit}{Fitted responses; either a vector of length \eqn{n} for scalar responses 
#' or a \eqn{q}-by-\eqn{n} matrix holding the quantile functions for distributional responses, 
#' evaluated on \code{qSup}.}
#' \item{qSup}{The input \code{qSup}.}
#' \item{beta}{A list of fitted coefficient functions, 
#' \eqn{j}-th element corresponding to \eqn{j}-th predictor in \code{X}. 
#' Elements take the following form: 
#' (1) a scalar, if the response and the corresponding predictor are both scalar;
#' (2) a vector holding the fitted coefficient function evaluated on \code{workGridX}, 
#' if the response is scalar and the corresponding predictor is distributional;
#' (3) a vector holding the fitted coefficient function evaluated on \code{workGridY}, 
#' if the response is distributional and the corresponding predictor is scalar;
#' (4) a matrix holding the fitted coefficient function, where (\code{j},\code{k})-th entry 
#' holds the value evaluated at (\code{workGridX[j]},\code{workGridY[k]}), 
#' if the response and the corresponding predictor are both distributional.}
#' \item{workGridX}{A list of elements, 
#' which are either vectors holding a working grid for elements of \code{beta} for a distributional predictor, 
#' or \code{NULL} for a scalar predictor.}
#' \item{workGridY}{A vector holding a working grid for elements of \code{beta} for a distributional response; 
#' \code{NULL} for a scalar response.}
#' \item{fpcaLogX}{A list of FPCA objects holding the FPCA output for each log mapped distributional predictors; 
#' Or \code{NULL} if all the predictors are scalar.}
#' \item{fpcaLogY}{An FPCA object holding the FPCA output for log mapped responses if the responses are distributional;
#' \code{NULL} if the responses are scalar.}
#' \item{X}{The input \code{X}.}
#' \item{Y}{The input \code{Y}.}
#' \item{outOfLogSpace}{A logical vector holding whether initial fitted log mapped distributional responses 
#' lie out of the log image space (\code{TRUE}) or not (\code{FALSE}); 
#' \code{NULL} if the responses are scalar.}
#' \item{optns}{A list of options actually used.}
#' \item{FPCAoptnsX}{A list of FPCA options actually used for distributional predictors; 
#' \code{NULL} if all the predictors are scalar.}
#' \item{FPCAoptnsY}{A list of FPCA options actually used for distributional responses; 
#' \code{NULL} if the responses are scalar.}
#' @references
#' \cite{Chen, Y., Lin, Z., & Müller, H.-G. (2021). "Wasserstein regression." Journal of the American Statistical Association, in press.}
#' \cite{Petersen, A., & Müller, H.-G. (2019). "Fréchet regression for random objects with Euclidean predictors." The Annals of Statistics, 47(2), 691--719.}
#' @examples 
#' \donttest{
#' # X = N( mu1, sigma1^2 ), Y = N( mu2, sigma2^2 )
#' # E( mu2 | mu1, sigma1 ) = E(mu2) + b11 ( mu1 - E(mu1) ) + b12 ( sigma1 - E(sigma1) )
#' # E( sigma2 | mu1, sigma1 ) = E(sigma2) + b21 ( mu1 - E(mu1) ) + b22 ( sigma1 - E(sigma1) )
#' # E(sigma2) + b21 ( mu1 - E(mu1) ) + b22 ( sigma1 - E(sigma1) ) >= 0 a.s.
#' # Then beta(s,t) = b11 + b12 ( s - E(mu1) ) / E(sigma1)  + b21 ( t - E(mu2) ) / E(sigma2) + 
#' # b22 ( ( s - E(mu1) ) / E(sigma1) ) ( ( t - E(mu2) ) / E(sigma2) )
#' bMat <- matrix( c( 0.5, 0.5, 0.5, -0.5 ), ncol = 2, byrow = TRUE )
#' n <- 100
#' # mu1 ~ Beta(2,2)
#' mu1 <- rbeta( n, 2, 2 ); Emu1 <- 0.5
#' # sigma1 ~ Uniform(0.5,1.5)
#' sigma1 <- runif( n, 0.5, 1.5 ); Esigma1 <- 1
#' # mu2 = E(mu2) + b11 ( mu1 - E(mu1) ) + b12 ( sigma1 - E(sigma1) ) + err_{mu2}
#' Emu2 <- 2
#' err_mu2 <- rnorm( n, mean = 0, sd = 0.1 ) # unexplained stochastic part
#' mu2True <- Emu2 + bMat[1,1] * ( mu1 - Emu1 ) + bMat[1,2] * ( sigma1 - Esigma1 )
#' mu2 <- mu2True + err_mu2
#' # sigma2 = E(sigma2) + b21 ( mu1 - E(mu1) ) + b22 ( sigma1 - E(sigma1) ) + err_{sigma2}
#' Esigma2 <- 2
#' err_sigma2 <- ( rbeta( n, 2, 2 ) - 0.5 ) / 5 # unexplained stochastic part
#' sigma2True <- Esigma2 + bMat[2,1] * ( mu1 - Emu1 ) + bMat[2,2] * ( sigma1 - Esigma1 )
#' sigma2 <- sigma2True + err_sigma2
#' 
#' nqSup <- 1000
#' qSup <- seq( 0, 1, length.out = nqSup+1 )
#' qSup <- ( qSup[-1] + qSup[-(nqSup+1)] ) / 2
#' X <- list(
#'   X1 = sapply( seq_len(n), function (i) { 
#'     qnorm( qSup, mean = mu1[i], sd = sigma1[i] ) 
#'   })
#' )
#' QXmean <- list(
#'   qnorm( qSup, mean = Emu1, sd = Esigma1 )
#' )
#' QYmean <- qnorm( qSup, mean = Emu2, sd = Esigma2 )
#' Y <- sapply( seq_len(n), function (i) { 
#'   qnorm( qSup, mean = mu2[i], sd = sigma2[i] ) 
#' })
#' res <- WR( X = X, Y = Y, qSup = qSup )
#' }
#' 
#' @importFrom fdapace FPCA
#' @export
#' 

WR <- function ( X, Y, qSup, optns = list(), FPCAoptnsX = list(), FPCAoptnsY = list() ) {
  # check X, Y, qSup----
  checkXY <- CheckData4WR( X = X, Y = Y, qSup = qSup )
  
  scalarRps <- checkXY$scalarRps
  anyScalarPdt <- checkXY$anyScalarPdt
  anyDistnlPdt <- checkXY$anyDistnlPdt
  indScalarPdt <- checkXY$indScalarPdt
  indDistnlPdt <- checkXY$indDistnlPdt
  n <- checkXY$n
  q <- checkXY$q
  
  # check and set up default values of control parameters in optns ----
  if ( scalarRps ) {
    optns <- SetOptns( optns = optns, Ymin = NULL, Ymax = NULL )
  } else {
    optns <- SetOptns( optns = optns, Ymin = min(Y), Ymax = max(Y) )
  }
  
  # set up default in FPCAoptnsX, FPCAoptnsY ----
  if ( anyDistnlPdt ) {
    FPCAoptnsX <- SetFPCAoptns( FPCAoptnsX )
  }
  if ( !scalarRps ) {
    FPCAoptnsY <- SetFPCAoptns( FPCAoptnsY )
  }
  
  # FPCA of Log(X) ----
  if ( anyDistnlPdt ) {
    QXmean <- lapply( X[ indDistnlPdt ], rowMeans )
    LogX <- lapply( indDistnlPdt, function (j) { X[[j]] - QXmean[[j]] })
    
    fpcaLogX <- lapply( LogX, function (LogXj) {
      FPCA(
        Ly = as.list(data.frame(LogXj)),
        Lt = rep( list(qSup), n ),
        optns = FPCAoptnsX
      )
    })
  }
  
  # FPCA of Log(Y) ----
  if ( !scalarRps ) {
    QYmean <- rowMeans(Y)
    LogY <- Y - QYmean
    
    fpcaLogY <- FPCA(
      Ly = as.list(data.frame(LogY)),
      Lt = rep( list(qSup), n ),
      optns = FPCAoptnsY
    )
  }
  
  # obtain input data for regression: scalar X and FPCs of LogX & scalar Y or FPCs of LogY ----
  if ( scalarRps ) {
    Ymean <- mean(Y)
    inputY <- list( Y - Ymean )
  } else {
    inputY <- as.list( as.data.frame( fpcaLogY$xiEst ) )
  }
  inputX <- list() # a list containing scalar predictors and FPCs of log mapped distributional predictors
  indPdt <- list()
  if ( anyScalarPdt ) {
    inputX[ indScalarPdt ] <- lapply( indScalarPdt, function (j) {
      X[[j]] - mean( X[[j]] )
    })
    indPdt[ indScalarPdt ] <- lapply( indScalarPdt, function(j) j )
  }
  if ( anyDistnlPdt ) {
    inputX[ indDistnlPdt ] <- lapply( fpcaLogX, function (res) res$xiEst )
    nFPCsX <- sapply( inputX[ indDistnlPdt ], ncol )
    indPdt[ indDistnlPdt ] <- lapply( seq_along( indDistnlPdt ), function (j) {
      rep( indDistnlPdt[j], nFPCsX[j] )
    })
  }
  indPdt <- Reduce( c, indPdt ) # a vector holding to which predictor each entry of coefficient estimates obtain by lm() corresponds;
  # i.e., to which predictor each row of \code{coefs} corresponds
  
  # regress scalar Y or FPCs of LogY on scalar X and FPCs of LogX ----
  regFPC <- lapply( inputY, function(y) {
    dfReg <- as.data.frame( c( list(y), inputX ) )
    names(dfReg)[1] <- 'y'
    lm( y ~ 0 + ., data = dfReg )
  })
  
  # obtain regression coefficient functions and fitted responses for X ----
  coefs <- sapply( regFPC, function (res) res$coefficients )
  beta <- Yfit <- list()
  # Yfit: initially a list, each element of which holds the part of fitted scalar responses or log distributional responses relevant to each predictor. 
  if ( scalarRps ) {
    ## for scalar responses ----
    if ( anyScalarPdt ) {
      tmp_beta <- coefs[ match( indScalarPdt, indPdt ) ]
      beta[ indScalarPdt ] <- as.list( as.data.frame( t( tmp_beta ) ) )
      for ( j in indScalarPdt ) {
        Yfit[[j]] <- inputX[[j]] * beta[[j]]
      }
    }
    if ( anyDistnlPdt ) {
      for ( j in seq_along(indDistnlPdt) ) {
        tmp_ind <- indDistnlPdt[j]
        tmp_coefs <- coefs[ indPdt == tmp_ind ]
        beta[[ tmp_ind ]] <- as.vector( fpcaLogX[[j]]$phi %*% tmp_coefs ) # support grid: fpcaLogX[[j]]$workGrid
        Yfit[[ tmp_ind ]] <- as.vector( inputX[[ tmp_ind ]] %*% tmp_coefs )
      }
    }
  } else {
    ## for distributional responses ----
    if ( anyScalarPdt ) {
      for ( j in indScalarPdt ) {
        tmp_coefs <- coefs[ which( indPdt == j ), ]
        beta[[j]] <- as.vector( fpcaLogY$phi %*% tmp_coefs ) # support grid: fpcaLogY$workGrid
        Yfit[[j]] <- matrix( beta[[j]], ncol = 1 ) %*% matrix( inputX[[j]], nrow = 1 ) # dim = length(fpcaLogY$workGrid) x n
      }
    }
    if ( anyDistnlPdt ) {
      for ( j in seq_along(indDistnlPdt) ) {
        tmp_ind <- indDistnlPdt[j]
        tmp_coefs <- coefs[ which( indPdt == tmp_ind ), , drop = FALSE ] # (l,k)-th entry: coefficient associated with l-th eigenfunction of LogXj and k-th eigenfunction of LogY
        beta[[ tmp_ind ]] <- fpcaLogX[[j]]$phi %*% tmp_coefs %*% t( fpcaLogY$phi ) # (s,t)-th entry: beta function evaluated at ( fpcaLogX[[j]]$workGrid[s], fpcaLogY$workGrid[t] )
        Yfit[[ tmp_ind ]] <- t( inputX[[ tmp_ind ]] %*% tmp_coefs %*% t( fpcaLogY$phi ) ) # (t,i)-th entry: part of fits relevant to this predictor evaluated at fpcaLogY$workGrid[t] for the i-th realization
      }
    }
  }
  # combine elements of Yfit to obtain fitted scalar responses or log distributional responses
  Yfit <- Reduce( '+', Yfit ) # (t,i)-th entry: part of fits relevant to this predictor evaluated at fpcaLogY$workGrid[t] for the i-th realization
  
  if ( scalarRps ) {
    ## add back mean of Y for scalar responses ----
    Yfit <- Yfit + Ymean
  } else {
    ## obtain support grid of log distributional responses ----
    if ( all.equal(qSup,fpcaLogY$workGrid) ) {
      LogYSup <- QYmean
    } else {
      LogYSup <- approx( x = qSup, y = QYmean, xout = fpcaLogY$workGrid )$y
    }
    ## map fitted log distributional responses back to distributions (quantile functions) ----
    LogYfit <- Yfit
    Yfit <- GetProj( LogYhat = LogYfit, LogSup = LogYSup )
    outOfLogSpace <- Yfit$outOfLogSpace
    Yfit <- Yfit$Yhat # quantile functions of fitted responses; length(fpcaLogY$workGrid) x n
    ## obtain fitted quantile functions evaluated on qSup for distributional responses ----
    if ( !all.equal(fpcaLogY$workGrid, qSup) ) {
      Yfit <- apply( Yfit, 2, function (y) {
        approx( x = fpcaLogY$workGrid, y = y, xout = qSup )$y
      })
    }
  }
  
  # working grid of beta ----
  if ( anyDistnlPdt ) {
    workGridX <- lapply( seq_along(fpcaLogX), function (j) {
      if ( all.equal( qSup, fpcaLogX[[j]]$workGrid ) ) {
        QXmean[[j]]
      } else {
        approx( x = qSup, y = QXmean[[j]], xout = fpcaLogX[[j]]$workGrid )$y
      }
    })
  } else {
    workGridX <- NULL
    fpcaLogX <- NULL
    FPCAoptnsX <- NULL
  }
  if ( !scalarRps ) {
    workGridY <- LogYSup
  } else {
    workGridY <- NULL
    fpcaLogY <- NULL
    FPCAoptnsY <- NULL
    outOfLogSpace <- NULL
  }
  
  res <- list(
    Yfit = Yfit,
    qSup = qSup,
    beta = beta,
    workGridX = workGridX, workGridY = workGridY,
    fpcaLogX = fpcaLogX, fpcaLogY = fpcaLogY,
    X = X, Y = Y,
    outOfLogSpace = outOfLogSpace,
    optns = optns, 
    FPCAoptnsX = FPCAoptnsX, FPCAoptnsY = FPCAoptnsY
  )
  class(res) <- 'WR'
  
  return(res)
}
