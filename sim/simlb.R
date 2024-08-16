source("code/lrem.R")
source("code/lcm.R")
source("code/kerFctn.R")
source("code/bwCV.R")

# Parallel computing
library(doSNOW)
cl <- makeCluster(50) 
registerDoSNOW(cl)
progress <- function(q) {
  if(q %% 10 == 0){
    cat(sprintf("%d runs are complete\n", q))
  }
}

Q <- 100
nVec <- c(50, 100, 200, 500, 1000)
nOut <- 200
eta0 <- 0
sigma0 <- 3
alpha <- 3
beta <- 0.5
tau <- 0.5
kappa <- 1

bw0 <- nVec^(-0.2)
bw <- nVec^(-0.2)
bwp <- nVec^(-0.2)
ise <- foreach(n = nVec, .combine = 'cbind') %:%
  foreach(icount(Q), .combine = 'c', .options.snow = list(progress = progress)) %dopar% {
    id <- which(nVec == n)
    lambdan <- 0.3 * n
    N <- NULL
    while (length(N) < n) {
      N <- c(N, rpois(1, lambda = lambdan))
      if (N[length(N)] < 2) {# == 0 enough for grem. Petersen and Müller (2019) won't work if N_i=1
        N <- N[-length(N)]
      }
    }
    M <- min(plcm(N), n * max(N), 5000)
    x <- runif(n, min = -1, max = 1)
    y0 <- list()
    y <- list()
    y0T <- list()
    yT <- list()
    yMean <- matrix(nrow = n, ncol = M - 1)
    for(i in 1:n) {
      mu <- rnorm(1, mean = eta0 + alpha * sin(pi * x[i]), sd = tau)
      sigma <- rgamma(1, shape = (sigma0 + beta * sin(pi * x[i]))^2 / kappa,
                      scale = kappa / (sigma0 + beta * sin(pi * x[i])))
      y0[[i]] <- mu + sigma * qbinom(1:(M - 1) / M, 5, 0.5)
      y[[i]] <- mu + sigma * rbinom(N[i], 5, 0.5)
      yMean[i, ] <- eta0 + alpha * sin(pi * x[i]) + (sigma0 + beta * sin(pi * x[i])) * qbinom(1:(M - 1) / M, 5, 0.5)
      
      k <- sample(c(-2, -1, 1, 2), 1)
      y0T[[i]] <- y0[[i]] - sin(k * y0[[i]]) / abs(k)
      yT[[i]] <- y[[i]] - sin(k * y[[i]]) / abs(k)
    }
    res0 <- lrem(y0, x, optns = list(bw = bw0[id]))
    res <- lrem(y, x, optns = list(bw = bw[id]))
    resp <- frechet::LocDenReg(xin = x, yin = y, optns = list(qSup = c(0, res$qfSupp), bwReg = bwp[id]))
    res0T <- lrem(y0T, x, optns = list(bw = bw0[id]))
    resT <- lrem(yT, x, optns = list(bw = bw[id]))
    respT <- frechet::LocDenReg(xin = x, yin = yT, optns = list(qSup = c(0, resT$qfSupp), bwReg = bwp[id]))
    2 * c(sum((res0$qf[, -M] - yMean)^2) / ((M - 1) * n), 
          sum((res$qf[, -M] - yMean)^2) / ((M - 1) * n), 
          sum((resp$qout[, -c(1, M)] - yMean)^2) / ((M - 1) * n), 
          sum((res0T$qf[, -M] - yMean)^2) / ((M - 1) * n), 
          sum((resT$qf[, -M] - yMean)^2) / ((M - 1) * n), 
          sum((respT$qout[, -c(1, M)] - yMean)^2) / ((M - 1) * n))# no 2 * if x ~ unif(0, 1)
  }
isel0 <- ise[seq(1, 6*Q, by = 6), ]
isel <- ise[seq(2, 6*Q, by = 6), ]
iselp <- ise[seq(3, 6*Q, by = 6), ]
isel0T <- ise[seq(4, 6*Q, by = 6), ]
iselT <- ise[seq(5, 6*Q, by = 6), ]
iselpT <- ise[seq(6, 6*Q, by = 6), ]
save(isel0, isel, iselp, isel0T, iselT, iselpT, file = 'data/iselb.RData')
stopCluster(cl)
