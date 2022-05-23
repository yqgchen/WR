# library(testthat)

test_that("Gaussian distributional time series - default", {
  # X_{t} = N( mu_{t}, sigma_{t}^2 )
  # mu_{t+1} = Emu + b11 * ( mu_{t} - Emu ) + b12 * ( sigma_{t} - Esigma ) + err_{mu,t+1}
  # sigma_{t+1} = Esigma + b21 * ( mu_{t} - Emu ) + b22 * ( sigma_{t} - Esigma ) + err_{sigma,t+1}
  bMat <- matrix( c( 1, 1, 1, -1 ) * 0.4, ncol = 2, byrow = TRUE )
  # mu_{1} ~ Beta(2,2), sigma_{1} ~ Uniform(0.5,1.5) + Esigma - 1
  # err ~ Uniform(-1,1) * M_err; M_err = 0.05
  # | mu_{t+1} - Emu | < | mu_{t} - Emu |
  # | sigma_{t+1} - Esigma | < | sigma_{t} - Esigma |
  set.seed(1)
  Emu <- 2
  Esigma <- 2
  M_err <- 0.05
  mu_c <- rbeta( 1, 2, 2 ) - Emu
  sigma_c <- runif( 1, 0.5, 1.5 ) - 1
  dat_c <- matrix( c(mu_c,sigma_c), ncol = 1 )
  n <- 100
  for ( i in 2:n ) {
    dat_c <- cbind( dat_c, bMat %*% dat_c[,i-1] + runif(2,-1,1) * M_err )
  }
  mu_c <- dat_c[1,]
  mu <- mu_c + Emu
  sigma_c <- dat_c[2,]
  sigma <- sigma_c + Esigma
  
  nqSup <- 1000
  qSup <- seq( 0, 1, length.out = nqSup+1 )
  qSup <- ( qSup[-1] + qSup[-(nqSup+1)] ) / 2
  
  Y <- sapply( seq_len(n), function (i) {
    qnorm( qSup, mean = mu[i], sd = sigma[i] )
  })
  res <- WAR( Y = Y, qSup = qSup )
  
  nPred <- 100
  datPred_c <- bMat %*% dat_c[,n]
  for ( i in 2:nPred ) {
    datPred_c <- cbind( datPred_c, bMat %*% datPred_c[,i-1] + runif(2,-1,1) * M_err )
  }
  muPred_c <- datPred_c[1, ]
  muPred <- muPred_c + Emu
  sigmaPred_c <- datPred_c[2,]
  sigmaPred <- sigmaPred_c + Esigma
  
  YpredTrue <- sapply( seq_len(nPred), function (i) {
    qnorm( qSup, mean = muPred[i], sd = sigmaPred[i] )
  })
  resPred <- predict( res, n.ahead = nPred )
  ASWD_Ypred <- mean( apply( ( resPred$Ypred - YpredTrue ) ^2, 2, pracma::trapz, x = qSup ) )
  expect_lt( ASWD_Ypred, 0.005 )
})

