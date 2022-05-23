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
  # mu1 ~ Beta(2,2)
  Emu1 <- 0.5
  mu1 <- rbeta( n, 2, 2 )
  # sigma1 ~ Uniform(0.5,1.5)
  Esigma1 <- 1
  sigma1 <- runif( n, 0.5, 1.5 )
  # mu2 = E(mu2) + b11 ( mu1 - E(mu1) ) + b12 ( sigma1 - E(sigma1) ) + err_{mu2}
  Emu2 <- 2
  err_mu2 <- rnorm( n, mean = 0, sd = 0.1 ) # unexplained stochastic part
  mu2True <- Emu2 + bMat[1,1] * ( mu1 - Emu1 ) + bMat[1,2] * ( sigma1 - Esigma1 )
  mu2 <- mu2True + err_mu2
  # sigma2 = E(sigma2) + b21 ( mu1 - E(mu1) ) + b22 ( sigma1 - E(sigma1) ) + err_{sigma2}
  Esigma2 <- 2
  err_sigma2 <- ( rbeta( n, 2, 2 ) - 0.5 ) / 5 # unexplained stochastic part
  sigma2True <- Esigma2 + bMat[2,1] * ( mu1 - Emu1 ) + bMat[2,2] * ( sigma1 - Esigma1 )
  sigma2 <- sigma2True + err_sigma2
  
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
  QYmean <- qnorm( qSup, mean = Emu2, sd = Esigma2 )
  Y <- sapply( seq_len(n), function (i) { 
    qnorm( qSup, mean = mu2[i], sd = sigma2[i] ) 
  })
  res <- WR( X = X, Y = Y, qSup = qSup )
  
  Ytrue <- sapply( seq_len(n), function (i) { 
    qnorm( qSup, mean = mu2True[i], sd = sigma2True[i] ) 
  })
  ASWD_Yfit <- mean( apply( ( res$Yfit - Ytrue )^2, 2, pracma::trapz, x = qSup ) )
  expect_lt( ASWD_Yfit, 1e-3 )
  
  # # parallel tranported hat{beta} evaluated on QXmean(qSup) x QYmean(qSup)
  # betaEst_QXmean_QYmean <- t( apply(
  #   apply( res$beta[[1]], 2, function (y) {
  #     approx( x = res$workGridX[[1]], y = y, xout = QXmean[[1]], rule = 2 )$y
  #   }),
  #   1,
  #   function (y) {
  #     approx( x = res$workGridY, y = y, xout = QYmean, rule = 2 )$y
  #   }
  # ) )
  # betaTrue_QXmean_QYmean <- bMat[1,1] + bMat[1,2] * qnorm(qSup) +
  #   bMat[2,1] * matrix( qnorm(qSup), ncol = nqSup, nrow = nqSup, byrow = TRUE ) +
  #   bMat[2,2] * qnorm(qSup) %*% t( qnorm(qSup) )
  # ISE_beta <- pracma::trapz(
  #   x = qSup,
  #   y = apply( ( betaEst_QXmean_QYmean - betaTrue_QXmean_QYmean )^2, 2, pracma::trapz, x = qSup )
  # )
  # expect_lt( ISE_beta, 0.03 )
})

test_that("D2S regression", {
  # E(Y|X,Z) = EY + < LogX, beta > + ( Z - EZ ) alpha
  
  # X = N( mu1, sigma1^2 )
  set.seed(1)
  n <- 100
  # mu1 ~ Beta(2,2)
  Emu1 <- 0.5
  mu1 <- rbeta( n, 2, 2 )
  # sigma1 ~ Uniform(0.5,1.5)
  Esigma1 <- 1
  sigma1 <- runif( n, 0.5, 1.5 )
  
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
  EZ <- 1
  X$Z <- Zc + EZ
  
  alpha <- -1.5
  EY <- 3
  sd_err <- 0.2
  err <- rnorm( n, mean = 0, sd = sd_err )
  Ytrue <- EY + Zc * alpha + betaLogX
  Y <- EY + Zc * alpha + betaLogX + err
  
  res <- WR( X = X, Y = Y, qSup = qSup )
  
  ISE_Yfit <- mean( ( res$Yfit - Ytrue )^2 )
  expect_lt( ISE_Yfit, sd_err/100 )
  
})

