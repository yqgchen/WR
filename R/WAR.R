#' Fit Wasserstein autoregressive models to distributional time series
#' @param Y A \eqn{q}-by-\eqn{n} matrix holding the quantile functions of 
#' \eqn{n} distributions evaluated on a common grid on \eqn{[0,1]} consisting of \eqn{q} points.
#' @param qSup A vector of length \eqn{q} holding the support grid of quantile functions.
#' @param optns A list of control parameters specified by \code{list(name=value)}. See 'Details'.
#' @param FPCAoptns A list of optional control parameters for functional principal 
#' component analysis (FPCA) of log mapped distributional predictors and responses, respectively. 
#' See 'Details' of \code{\link[fdapace]{FPCA}}. For example, one can specify 
#' the methods of choosing numbers of functional principal components (FPCs) through \code{methodSelectK}; 
#' set up Fraction-of-Variance-Explained (FVE) thresholds through \code{FVEthreshold} 
#' when using FVE to choose numbers of FPCs. 
#' Default of \code{methodSelectK}: \code{'FVE'}; 
#' Default of \code{FVEthreshold}: 0.95 (different from \code{\link[fdapace]{FPCA}}).
#' Default of \code{useBinnedData}: \code{'OFF'} (different from \code{\link[fdapace]{FPCA}}).
#' In addition, if \code{optns$CVfold} is specified (i.e., not \code{NULL}), 
#' \code{methodSelectK = 'FVE'} and \code{FVEthreshold = 0.9999}. 
#' @param aic,order.max,method,... arguments for \code{\link[stats]{ar}}. 
#' Default: \code{aic = FALSE}, \code{order.max = 1}, \code{method = "yule-walker"}. 
#' Note that \code{demean} should not be input and is set to be \code{FALSE}.
#' @details Available control options are
#' \describe{
#' \item{methodProj}{The projection method if the predictions of log mapped distributional responses lie out of the log image space: 
#' \code{'log'} (default, the method as per Chen, Lin, and Müller, 2021); 
#' \code{'qt'} (the method as per Petersen & Müller, 2019).}
#' \item{CVfold}{Either the number of folds if using cross validation to choose the number of FPCs, 
#' which can be any positive integer up to \eqn{n}, 
#' with suggested values: \eqn{n} (leave-one-out) if \eqn{n\le 30} and \eqn{10} (10-fold) if \eqn{n > 30}; 
#' Or \code{NULL} if using other methods to choose the number of FPCs.}
#' \item{lower}{Lower bound of the support of distributions, only relevant if \code{methodProj} is \code{'qt'}. 
#' Default: \code{NULL}, i.e., no finite lower bound.}
#' \item{upper}{Upper bound of the support of distributions, only relevant if \code{methodProj} is \code{'qt'}. 
#' Default: \code{NULL}, i.e., no finite upper bound.}
#' }
#' @return A list of class \code{'WAR'} with the following elements:
#' \item{Yfit}{A \eqn{q}-by-\eqn{n} matrix holding the fitted quantile functions of the \eqn{n} distributions.}
#' \item{qSup}{The support of quantile functions, same as the input.}
#' \item{beta}{A matrix holding the fitted coefficient function, where the (\code{j},\code{k})-th entry 
#' holds the value evaluated at (\code{workGrid[j]},\code{workGrid[k]}).}
#' \item{workGrid}{A vector holding a working grid for \code{beta}.}
#' \item{arFPCs}{An ar object returned by \code{\link[stats]{ar}}, holding the fitted autoregressive model to the FPCs.}
#' \item{order}{The order of the fitted Wasserstein autoregressive model.}
#' \item{fpcaLogY}{An FPCA object holding the FPCA output for log mapped distributions.}
#' \item{Y}{The input \code{Y}.}
#' \item{outOfLogSpace}{A logical vector holding whether initial fitted log mapped distributions 
#' lie out of the log image space (\code{TRUE}) or not (\code{FALSE}).}
#' \item{optns}{A list of options actually used.}
#' \item{FPCAoptns}{A list of FPCA options actually used.}
#' @references
#' \cite{Chen, Y., Lin, Z., & Müller, H.-G. (2021). "Wasserstein regression." Journal of the American Statistical Association, in press.}
#' \cite{Petersen, A., & Müller, H.-G. (2019). "Fréchet regression for random objects with Euclidean predictors." The Annals of Statistics, 47(2), 691--719.}
#' @examples 
#' \donttest{
#' # X_{t} = N( mu_{t}, sigma_{t}^2 )
#' # mu_{t+1} = Emu + b11 * ( mu_{t} - Emu ) + b12 * ( sigma_{t} - Esigma ) + err_{mu,t+1}
#' # sigma_{t+1} = Esigma + b21 * ( mu_{t} - Emu ) + b22 * ( sigma_{t} - Esigma ) + err_{sigma,t+1}
#' bMat <- matrix( c( 1, 1, 1, -1 ) * 0.4, ncol = 2, byrow = TRUE )
#' # mu_{1} ~ Beta(2,2), sigma_{1} ~ Uniform(0.5,1.5) + Esigma - 1
#' # err ~ Uniform(-1,1) * M_err; M_err = 0.05
#' # | mu_{t+1} - Emu | < | mu_{t} - Emu |
#' # | sigma_{t+1} - Esigma | < | sigma_{t} - Esigma |
#' set.seed(1)
#' Emu <- 2
#' Esigma <- 2
#' M_err <- 0.05
#' mu_c <- rbeta( 1, 2, 2 ) - Emu
#' sigma_c <- runif( 1, 0.5, 1.5 ) - 1
#' dat_c <- matrix( c(mu_c,sigma_c), ncol = 1 )
#' n <- 100
#' for ( i in 2:n ) {
#'   dat_c <- cbind( dat_c, bMat %*% dat_c[,i-1] + runif(2,-1,1) * M_err )
#' }
#' mu_c <- dat_c[1,]
#' mu <- mu_c + Emu
#' sigma_c <- dat_c[2,]
#' sigma <- sigma_c + Esigma
#' 
#' nqSup <- 1000
#' qSup <- seq( 0, 1, length.out = nqSup+1 )
#' qSup <- ( qSup[-1] + qSup[-(nqSup+1)] ) / 2
#' 
#' Y <- sapply( seq_len(n), function (i) {
#'   qnorm( qSup, mean = mu[i], sd = sigma[i] )
#' })
#' res <- WAR( Y = Y, qSup = qSup )
#' }
#' 
#' @importFrom fdapace FPCA
#' @importFrom MASS ginv
#' @import stats
#' @export
#' 
WAR <- function (Y, qSup, optns = list(), FPCAoptns = list(), aic = FALSE, order.max = 1, method = "yule-walker", ...){
  CheckData4WR( X = NULL, Y = Y, qSup = qSup )
  n <- ncol(Y)
  q <- nrow(Y)
  
  optns <- SetOptns( optns = optns, Ymin = min(Y), Ymax = max(Y) )
  
  FPCAoptns <- SetFPCAoptns( FPCAoptns )
  if ( !is.null( optns$CVfold ) ) {
    if ( optns$CVfold > n | optns$CVfold < 0 ) {
      stop ('Mis-specified optns$CVfold.')
    }
    FPCAoptns$methodSelectK <- 'FVE'
    FPCAoptns$FVEthreshold <- 0.9999
  }
  
  QYmean <- rowMeans(Y)
  LogY <- Y - QYmean
  colnames(LogY) <- 1:ncol(LogY)
  
  fpcaLogY <- FPCA(
    Ly = as.list(data.frame(LogY)),
    Lt = lapply(1:ncol(LogY), function(i) qSup),
    optns = FPCAoptns
  )
  
  nFPCs <- length( fpcaLogY$cumFVE )
  
  # use CV to choose # of FPCs ----
  if ( !is.null(optns$CVfold) ) {
    CVfold <- optns$CVfold
    ind_pert <- sample(seq_len(n-1), n-1)
    ind_partition <- findInterval( x = seq_len(n-1), seq(1,n-1,length.out = CVfold + 1)[1:CVfold] )
    cv_partition <- lapply( seq_len( CVfold ), function (k) {
      ind_pert[ ind_partition == k ]
    })
    
    MSPE_cv <- sapply( seq_len(nFPCs), function(K) {
      # require(MASS)
      fpcsTS_cv <- fpcaLogY$xiEst[,1:K,drop = FALSE]
      Cov1 <- lapply(2:n, function(t) {
        fpcsTS_cv[t,] %*% t(fpcsTS_cv[t-1,])
      })
      Cov0 <- lapply(1:(n-1), function(t) {
        fpcsTS_cv[t,] %*% t(fpcsTS_cv[t,])
      })
      LogYpred_cv <- lapply(1:length(cv_partition), function(k) {
        coefs_cv <- Reduce("+", Cov1[-cv_partition[[k]]]) %*% ginv(Reduce("+", Cov0[-cv_partition[[k]]]))
        fpcsYpred_cv <- t(coefs_cv %*% t(fpcsTS_cv[cv_partition[[k]],,drop = FALSE]))
        return( fpcaLogY$phi[,1:K,drop = FALSE] %*% t(fpcsYpred_cv) + fpcaLogY$mu )
      })
      LogYpred_cv <- Reduce(cbind,LogYpred_cv)
      LogYpred_cv <- LogYpred_cv[, order(unlist(cv_partition))]
      if ( !all.equal( qSup, fpcaLogY$workGrid) ) {
        LogYpred_cv <- sapply(1:ncol(LogYpred_cv), function(i) {
          approx(x = fpcaLogY$workGrid, y = LogYpred_cv[,i], xout = qSup, rule = 2)$y
        })
      }
      return( mean( apply( ( LogYpred_cv - LogY[,-1] )^2, 2, pracma::trapz, x = qSup ) ) )
    })
    nFPCs <- which.min(MSPE_cv)
    
    ## update fpcaLogY according to nFPCs ----
    # if ( nFPCs < length( fpcaLogY$cumFVE ) ) {
    #   fpcaLogY$lambda <- fpcaLogY$lambda[1:nFPCs]
    #   fpcaLogY$phi <- fpcaLogY$phi[,1:nFPCs,drop = FALSE]
    #   fpcaLogY$xiEst <- fpcaLogY$xiEst[,1:nFPCs,drop = FALSE]
    #   fpcaLogY$cumFVE <- fpcaLogY$cumFVE[1:nFPCs]
    #   fpcaLogY$FVE <- fpcaLogY$cumFVE[nFPCs]
    #   fpcaLogY$criterionValue <- NULL
    # }
  }
  
  # treat FPCs as multivariate time series ----
  fpcsTS <- ts( data = fpcaLogY$xiEst[,1:nFPCs,drop = FALSE], frequency = 1, start = 1 )
  # fit VAR on fpcsTS ----
  arFPCs <- ar( fpcsTS, aic = aic, order.max = order.max, method = method, demean = FALSE, ... )
  # obtain fitted FPCs starting from the second time point ----
  fpcsYfit <- fpcsTS[-1,,drop = FALSE] - arFPCs$resid[-1,,drop = FALSE]
  # obtain fitted log maps and their support grid ----
  LogYfit <- fpcaLogY$phi[,1:nFPCs,drop = FALSE] %*% t(fpcsYfit) + fpcaLogY$mu
  if ( all.equal( qSup, fpcaLogY$workGrid ) ) {
    LogSup <- QYmean
  } else {
    LogSup <- approx( x = qSup, y = QYmean, xout = fpcaLogY$workGrid )$y
  }
  # map fitted log maps back to distributions ----
  Yfit <- GetProj( LogYhat = LogYfit, LogSup = LogSup )
  outOfLogSpace <- Yfit$outOfLogSpace
  Yfit <- Yfit$Yhat
  # obtain fitted quantile functions on qSup ----
  if ( !all.equal(fpcaLogY$workGrid, qSup) ) {
    Yfit <- apply( Yfit, 2, function (y) {
      approx( x = fpcaLogY$workGrid, y = y, xout = qSup )$y
    })
  }
  
  # obtain coefficient function ----
  coefs <- arFPCs$ar[1,,]
  beta <- fpcaLogY$phi[,1:nFPCs,drop = FALSE] %*% coefs %*% t( fpcaLogY$phi[,1:nFPCs,drop = FALSE] )
  workGrid <- LogSup
  
  res <- list(
    Yfit = Yfit,
    qSup = qSup,
    beta = beta,
    workGrid = workGrid,
    arFPCs = arFPCs,
    order = arFPCs$order,
    fpcaLogY = fpcaLogY,
    Y = Y,
    outOfLogSpace = outOfLogSpace,
    optns = optns, FPCAoptns = FPCAoptns
  )
  class(res) <- 'WAR'
  return (res)
}
