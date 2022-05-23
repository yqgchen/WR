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
  
  ASWD_Yfit <- mean( apply( ( res$Yfit - Y[,-1] ) ^2, 2, pracma::trapz, x = qSup ) )
  expect_lt( ASWD_Yfit, 0.002 )
  
})

test_that("Gaussian distributional time series - using CV to choose numbers of FPCs", {
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
  res <- WAR( Y = Y, qSup = qSup, optns = list( CVfold = 10 ) )
  
  ASWD_Yfit <- mean( apply( ( res$Yfit - Y[,-1] ) ^2, 2, pracma::trapz, x = qSup ) )
  expect_lt( ASWD_Yfit, 0.002 )
})
