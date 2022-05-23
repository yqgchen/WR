# library(testthat)

test_that("D2D regression between Gaussian distributions w/ unexplained stochastic parts", {
  # X = N( mu1, sigma1^2 ), Y = N( mu2, sigma2^2 )
  # E( mu2 | mu1, sigma1 ) = E(mu2) + b11 ( mu1 - E(mu1) ) + b12 ( sigma1 - E(sigma1) )
  # E( sigma2 | mu1, sigma1 ) = E(sigma2) + b21 ( mu1 - E(mu1) ) + b22 ( sigma1 - E(sigma1) )
  # E(sigma2) + b21 ( mu1 - E(mu1) ) + b22 ( sigma1 - E(sigma1) ) >= 0 a.s.
  # Then beta(s,t) = b11 + b12 ( s - E(mu1) ) / E(sigma1)  + b21 ( t - E(mu2) ) / E(sigma2) + b22 [ ( s - E(mu1) ) / E(sigma1) ] [ ( t - E(mu2) ) / E(sigma2) ]
  
  bMat <- matrix( c( 0.5, 0.5, 0.5, -0.5 ), ncol = 2, byrow = TRUE )
  
  set.seed(1)
  n <- 100
  nPred <- 100
  # mu1 ~ Beta(2,2)
  Emu1 <- 0.5
  mu1 <- rbeta( n, 2, 2 )
  mu1pred <- rbeta( nPred, 2, 2 )
  # sigma1 ~ Uniform(0.5,1.5)
  Esigma1 <- 1
  sigma1 <- runif( n, 0.5, 1.5 )
  sigma1pred <- runif( nPred, 0.5, 1.5 )
  # mu2 = E(mu2) + b11 ( mu1 - E(mu1) ) + b12 ( sigma1 - E(sigma1) ) + err_{mu2}
  Emu2 <- 2
  err_mu2 <- rnorm( n, mean = 0, sd = 0.1 ) # unexplained stochastic part
  mu2True <- Emu2 + bMat[1,1] * ( mu1 - Emu1 ) + bMat[1,2] * ( sigma1 - Esigma1 )
  mu2predTrue <- Emu2 + bMat[1,1] * ( mu1pred - Emu1 ) + bMat[1,2] * ( sigma1pred - Esigma1 )
  mu2 <- mu2True + err_mu2
  # sigma2 = E(sigma2) + b21 ( mu1 - E(mu1) ) + b22 ( sigma1 - E(sigma1) ) + err_{sigma2}
  Esigma2 <- 2
  err_sigma2 <- ( rbeta( n, 2, 2 ) - 0.5 ) / 5 # unexplained stochastic part
  sigma2True <- Esigma2 + bMat[2,1] * ( mu1 - Emu1 ) + bMat[2,2] * ( sigma1 - Esigma1 )
  sigma2predTrue <- Esigma2 + bMat[2,1] * ( mu1pred - Emu1 ) + bMat[2,2] * ( sigma1pred - Esigma1 )
  sigma2 <- sigma2True + err_sigma2
  
  nqSup <- 1000
  qSup <- seq( 0, 1, length.out = nqSup+1 )
  qSup <- ( qSup[-1] + qSup[-(nqSup+1)] ) / 2
  X <- list(
    X1 = sapply( seq_len(n), function (i) { 
      qnorm( qSup, mean = mu1[i], sd = sigma1[i] ) 
    })
  )
  Y <- sapply( seq_len(n), function (i) { 
    qnorm( qSup, mean = mu2[i], sd = sigma2[i] ) 
  })
  res <- WR( X = X, Y = Y, qSup = qSup )
  
  Xpred <- list(
    X1 = sapply( seq_len(nPred), function (i) { 
      qnorm( qSup, mean = mu1pred[i], sd = sigma1pred[i] ) 
    })
  )
  resPred <- predict( res, Xpred = Xpred )
  YpredTrue <- sapply( seq_len(nPred), function (i) { 
    qnorm( resPred$qSup, mean = mu2predTrue[i], sd = sigma2predTrue[i] ) 
  })
  MISE_Ypred <- mean( apply( ( resPred$Ypred - YpredTrue )^2, 2, pracma::trapz, x = qSup ) )
  
  expect_lt( MISE_Ypred, 1e-3 )
})

test_that("D2S regression", {
  # E(Y|X,Z) = EY + < LogX, beta > + ( Z - EZ ) alpha
  
  # X = N( mu1, sigma1^2 )
  set.seed(1)
  n <- 100
  nPred <- 100
  # mu1 ~ Beta(2,2)
  Emu1 <- 0.5
  mu1 <- rbeta( n, 2, 2 )
  mu1pred <- rbeta( nPred, 2, 2 )
  # sigma1 ~ Uniform(0.5,1.5)
  Esigma1 <- 1
  sigma1 <- runif( n, 0.5, 1.5 )
  sigma1pred <- runif( nPred, 0.5, 1.5 )
  
  nqSup <- 1000
  qSup <- seq( 0, 1, length.out = nqSup+1 )
  qSup <- ( qSup[-1] + qSup[-(nqSup+1)] ) / 2
  X <- list(
    X1 = sapply( seq_len(n), function (i) { 
      qnorm( qSup, mean = mu1[i], sd = sigma1[i] ) 
    })
  )
  QXmean <- list(
    qnorm( qSup, mean = Emu1, sd = Esigma1 )
  )
  
  # beta(s) = b1 + b2 * ( s - Emu1 ) / Esigma1
  # beta( F^{-1}_{Xmean}(u) ) = b1 + b2 Phi^{-1}(u)
  bVec <- c( 2, -1 )
  beta_QXmean <- bVec[1] + bVec[2] * qnorm(qSup)
  
  betaLogX <- lapply( seq_along(X), function (j) {
    colMeans( ( X[[j]] - QXmean[[j]] ) * beta_QXmean )
  })
  betaLogX <- Reduce( `+`, betaLogX )
  
  # Z - EZ ~ N( 0, 2^2 )
  Zc <- rnorm( n, mean = 0, sd = 2 )
  Zpredc <- rnorm( nPred, mean = 0, sd = 2 )
  EZ <- 1
  X$Z <- Zc + EZ
  
  alpha <- -1.5
  EY <- 3
  sd_err <- 0.2
  err <- rnorm( n, mean = 0, sd = sd_err )
  Ytrue <- EY + Zc * alpha + betaLogX
  Y <- EY + Zc * alpha + betaLogX + err
  
  res <- WR( X = X, Y = Y, qSup = qSup )
  
  Xpred <- list(
    X1 = sapply( seq_len(nPred), function (i) { 
      qnorm( qSup, mean = mu1pred[i], sd = sigma1pred[i] ) 
    })
  )
  betaLogXpred <- lapply( seq_along(Xpred), function (j) {
    colMeans( ( Xpred[[j]] - QXmean[[j]] ) * beta_QXmean )
  })
  betaLogXpred <- Reduce( `+`, betaLogXpred )
  Xpred$Z <- Zpredc + EZ
  YpredTrue <- EY + Zpredc * alpha + betaLogXpred
  
  resPred <- predict( res, Xpred = Xpred )
  
  ISE_Ypred <- mean( ( resPred$Ypred - YpredTrue )^2 )
  expect_lt( ISE_Ypred, sd_err/100 )
  
})

