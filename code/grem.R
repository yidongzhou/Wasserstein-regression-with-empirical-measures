#' @title Global Regression with Empirical Measures (REM)
#' @description Global regression for empirical measures with Euclidean predictors.
#' @param y A list of \eqn{n} empirical measures, represented as vectors of
#' observed values. The \eqn{i}th vector contains \eqn{N_i} observations for
#' the \eqn{i}th individual.
#' @param x An n by p matrix or data frame of predictors. It can be a vector of 
#' length n if p = 1.
#' @param xOut An nOut by p matrix or data frame of output predictor levels.
#' It can be a vector of length p if nOut = 1 or a vector of length nOut if p = 1
#' or a scalar if both p and nOut are equal to 1. Default is \code{NULL}.
#' @param optns A list of options control parameters specified by \code{list(name = value)}. 
#' See `Details'.
#' @details Available control options are
#' \describe{
#' \item{lower}{The lower bound of the support of the measure. Default is \code{NULL}.}
#' \item{upper}{The upper bound of the support of the measure. Default is \code{NULL}.}
#' }
#' @return A \code{rem} object, which is a list containing the following components:
#' \item{qf}{A matrix holding the quantile functions corresponding to \code{x}. 
#' Each row gives a quantile and the domain grid is given in \code{qfSupp}.}
#' \item{qfSupp}{A vector giving the domain grid of \code{qf}.}
#' \item{qp}{A matrix holding the quantile functions corresponding to 
#' \code{xOut}. Each row gives a quantile and the domain grid is given in 
#' \code{qpSupp}. Included if \code{xOut} is not \code{NULL}.}
#' \item{qpSupp}{A vector giving the domain grid of \code{qp}. 
#' Included if \code{xOut} is not \code{NULL}.}
#' \item{RSquare}{Fréchet coefficient of determination.}
#' \item{adjRSquare}{Adjusted Fréchet coefficient of determination.}
#' \item{residuals}{Wasserstein distances between the empirical and fitted measures.}
#' \item{y}{The empirical measures used.}
#' \item{x}{The predictors used.}
#' \item{xOut}{The output predictor levels used.}
#' \item{optns}{The control options used.}
#' @references
#' \itemize{
#' \item \cite{Zhou, Y. and Müller, H.G., 2023. Wasserstein Regression with Empirical Measures and Density Estimation for Sparse Data. arXiv preprint arXiv:2308.12540.}
#' \item \cite{Petersen, A. and Müller, H.-G. (2019). Fréchet regression for random objects with Euclidean predictors. The Annals of Statistics, 47(2), 691--719.}
#' }
#' @export

grem <- function(y = NULL,
                 x = NULL,
                 xOut = NULL,
                 optns = list()) {
  if (is.null(y) | is.null(x)) {
    stop("requires the input of both y and x")
  }
  if (!is.matrix(x)) {
    if (is.data.frame(x) | is.vector(x)) {
      x <- as.matrix(x)
    } else {
      stop("x must be a matrix or a data frame or a vector")
    }
  }
  n <- nrow(x) # number of observations
  p <- ncol(x) # number of covariates
  if (!is.list(y)) {
    stop("y must be a list")
  }
  if (length(y) != n) {
    stop("the number of rows in x must be the same as the number of empirical measures in y")
  }
  if (!is.null(xOut)) {
    if (!is.matrix(xOut)) {
      if (is.data.frame(xOut)) {
        xOut <- as.matrix(xOut)
      } else if (is.vector(xOut)) {
        if (p == 1) {
          xOut <- as.matrix(xOut)
        } else {
          xOut <- t(xOut)
        }
      } else {
        stop("xOut must be a matrix or a data frame or a vector")
      }
    }
    if (ncol(xOut) != p) {
      stop("x and xOut must have the same number of columns")
    }
    nOut <- nrow(xOut) # number of predictions
  } else {
    nOut <- 0
  }

  N <- sapply(y, length)
  y <- lapply(1:n, function(i) {
    sort(y[[i]])
  }) # sort observed values
  
  M <- min(plcm(N), n * max(N), 5000) # least common multiple of N_i
  yM <- t(sapply(1:n, function(i) {
    residual <- M %% N[i]
    if(residual) {
      sort(c(rep(y[[i]], each = M %/% N[i]), sample(y[[i]], residual)))
    } else {
      rep(y[[i]], each = M %/% N[i])
    }
  })) # n by M

  # initialization of OSQP solver
  A <- cbind(diag(M), rep(0, M)) + cbind(rep(0, M), -diag(M))
  if (!is.null(optns$upper) &
    !is.null(optns$lower)) {
    # if lower & upper are neither NULL
    l <- c(optns$lower, rep(0, M - 1), -optns$upper)
  } else if (!is.null(optns$upper)) {
    # if lower is NULL
    A <- A[, -1]
    l <- c(rep(0, M - 1), -optns$upper)
  } else if (!is.null(optns$lower)) {
    # if upper is NULL
    A <- A[, -ncol(A)]
    l <- c(optns$lower, rep(0, M - 1))
  } else {
    # if both lower and upper are NULL
    A <- A[, -c(1, ncol(A))]
    l <- rep(0, M - 1)
  }
  # P <- as(diag(M), "sparseMatrix")
  # A <- as(t(A), "sparseMatrix")
  P <- diag(M)
  A <- t(A)
  q <- rep(0, M)
  u <- rep(Inf, length(l))
  model <-
    osqp::osqp(
      P = P,
      q = q,
      A = A,
      l = l,
      u = u,
      osqp::osqpSettings(max_iter = 1e05, eps_abs = 1e-05, eps_rel = 1e-05, verbose = FALSE)
    )

  xMean <- colMeans(x)
  invVa <- solve(var(x) * (n - 1) / n)
  wc <-
    t(apply(x, 1, function(xi) {
      t(xi - xMean) %*% invVa
    })) # n by p
  if (nrow(wc) != n) {
    wc <- t(wc)
  } # for p=1
  
  qf <- matrix(nrow = n, ncol = M)
  residuals <- rep.int(0, n)
  totVa <- sum((scale(yM, scale = FALSE))^2) / M
  for (i in 1:n) {
    w <- apply(wc, 1, function(wci) {
      1 + t(wci) %*% (x[i, ] - xMean)
    })
    qNew <- apply(yM, 2, weighted.mean, w) # M
    if (any(w < 0)) {
      # if negative weights exist, project
      model$Update(q = -qNew)
      qNew <- sort(model$Solve()$x)
    }
    if (!is.null(optns$upper)) {
      qNew <- pmin(qNew, optns$upper)
    }
    if (!is.null(optns$lower)) {
      qNew <- pmax(qNew, optns$lower)
    }
    qf[i, ] <- qNew
    residuals[i] <- sqrt(sum((yM[i, ] - qf[i, ])^2) / M)
  }
  qfSupp <- 1:M / M
  resVa <- sum(residuals^2)
  RSquare <- 1 - resVa / totVa
  adjRSquare <- RSquare - (1 - RSquare) * p / (n - p - 1)

  if (nOut > 0) {
    qp <- matrix(nrow = nOut, ncol = M)
    for (i in 1:nOut) {
      w <- apply(wc, 1, function(wci) {
        1 + t(wci) %*% (xOut[i, ] - xMean)
      })
      qNew <- apply(yM, 2, weighted.mean, w) # M
      if (any(w < 0)) {
        # if negative weights exist
        model$Update(q = -qNew)
        qNew <- sort(model$Solve()$x)
      }
      if (!is.null(optns$upper)) {
        qNew <- pmin(qNew, optns$upper)
      }
      if (!is.null(optns$lower)) {
        qNew <- pmax(qNew, optns$lower)
      }
      qp[i, ] <- qNew
    }
    qpSupp <- 1:M / M
    
    res <-
      list(
        qf = qf,
        qfSupp = qfSupp,
        qp = qp,
        qpSupp = qpSupp,
        RSquare = RSquare,
        adjRSquare = adjRSquare,
        residuals = residuals,
        y = y,
        x = x,
        xOut = xOut,
        optns = optns
      )
  } else {
    res <- list(
      qf = qf,
      qfSupp = qfSupp,
      RSquare = RSquare,
      adjRSquare = adjRSquare,
      residuals = residuals,
      y = y,
      x = x,
      optns = optns
    )
  }

  class(res) <- "rem"
  res
}
