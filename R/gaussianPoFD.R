#' Gaussian Partially Observed Functional Data
#'
#' Generates samples of partially observed gaussian functions following different censoring regimes.
#'
#' @param n total number of functional observations
#' @param p total number of points observed for each function
#' @param type type of partially observed data. Options are "sparse", "interval" and "common". See Elías et al (2020).
#' @param observability mean observed proportion of the domain where each function is observed.
#' @param ninterval if type = "interval", n_interval is an integer with the number of observed intervals 1, 2, 3...
#' Large values of this parameter requires a large parameter p to guarantee the observability level.
#'
#' @return a list containing two elements 1) a functional sample and 2) the same sample of functions but
#' partially observed following one of the schemes described in the argument type.
#'
#' @importFrom MASS mvrnorm
#' @importFrom FastGP rcpp_rmvnorm
#'
#' @references Elías, Antonio, Jiménez, Raúl, Paganoni, Anna M. and Sangalli, Laura M. (2020). Integrated Depths for Partially Observed Functional Data.
#' @examples
#'
#' gaussian_pofd <- gaussian_PoFD(n=100, p=200, type="sparse", observability=0.5)
#'
#' @export
gaussian_PoFD <- function(n, p, type, observability, ninterval){
  #parameters
  time_grid <- seq(0, 1, length.out = p)
  sigmaPeriodic <- Cov_Periodic(time_grid, time_grid, sigma = 3, p = 1 , l = 0.5)
  sigmaExpo <- Cov_exponential(time_grid, time_grid, alpha = 0.5, beta = 5)

  #provides the random mean
  centerline <- MASS::mvrnorm(1, rep(0, p), sigmaPeriodic)
  #provides the random sample
  dataY <- FastGP::rcpp_rmvnorm(n, sigmaExpo, centerline)
  data <- t(dataY)
  colnames(data) <- as.character(c(1:n))
  rownames(data) <- round(time_grid, digits=5)

  podata <- switch(type, "sparse" = randomPointsObserved(data, observability) , "common" = extremesCensoring(data, observability), "interval" = nObservedAtRandom(data, observability, ninterval))

  return(list(fd = data , pofd = podata))
}

Cov_exponential <- function(X1, X2, alpha = 1, beta = 1){
  x.aux <- expand.grid(i = X1, j = X2)

  cov <- alpha*exp(-beta*abs(x.aux$i - x.aux$j))

  Sigma <- matrix(cov, nrow = length(X1))
  return(Sigma)
}

Cov_Periodic <- function(X1, X2, sigma = 1, p = 1, l = 1) {
  #p = period, l = wiggles, sigma = noise
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      Sigma[i,j] <- sigma*exp(-(2*(sin(pi*abs(X1[i]-X2[j])/(p)))^2)/(l^2))
    }
  }
  return(Sigma)
}

randomPointsObserved <- function(data, observedProportion){
  P <- dim(data)[1]
  N <- dim(data)[2]
  pObserved <- ceiling(P*observedProportion)

  sparseGrid <- matrix(NA, ncol = N, nrow = P)
  for(i in c(1:N)){
    auxLocation <- sample(c(1:P), size = pObserved)
    sparseGrid[auxLocation, i]<- data[auxLocation,i]
  }

  colnames(sparseGrid) <- c(1:N)
  rownames(sparseGrid) <- rownames(data)

  return(sparseGrid)
}

extremesCensoring <- function(data, observedProportion){
  P <- dim(data)[1]
  N <- dim(data)[2]

  sparseGrid <- matrix(NA, ncol = N, nrow = P)

  A_aux <- c(1:round(P/2-P*observedProportion/4, digits = 0))
  B_aux <- c(round(P/2+P*observedProportion/4, digits = 0):P)

  A <- sample(A_aux, N, replace = TRUE);
  B <- sample(B_aux, N, replace = TRUE)

  for(i in c(1:N)){
    sparseGrid[A[i]:B[i],i] <- data[A[i]:B[i],i]
  }

  colnames(sparseGrid) <- c(1:N)
  rownames(sparseGrid) <- rownames(data)

  return(sparseGrid)
}

nObservedAtRandom <- function(data, observedProportion, nInterval){
  P <- dim(data)[1]
  N <- dim(data)[2]

  sparseGrid <- matrix(NA, ncol = N, nrow = P)
  for(i in c(1:N)){

    auxRightExtreme <- sample(c(ceiling(P*observedProportion), floor(P*observedProportion)), size = 1) #if even

    observedIntervals <- c(1, sort(sample(2:(auxRightExtreme-1), size = nInterval-1, replace = FALSE)), auxRightExtreme)

    while(sum(diff(observedIntervals) !=1 ) != nInterval){
      observedIntervals <- c(1, sort(sample(2:(auxRightExtreme-1), size = nInterval-1, replace = FALSE)), auxRightExtreme)
    }

    missingIntervals <- c(auxRightExtreme+1, sort(sample(c(auxRightExtreme+1):P, size = nInterval, replace = FALSE)), P)

    locationObservedIntervals <- seq(2, 2*nInterval, by = 2)
    locationMissingIntervals <- seq(1, 2*(nInterval+1), by = 2)


    missing <-  c(1:c(1+diff(missingIntervals)[1]))
    observed <- c((max(missing)+1):(diff(observedIntervals)[1] + max(missing)+1))

    for(j in c(2:max(length(locationObservedIntervals), length(locationMissingIntervals)))){
      if(length(locationMissingIntervals) >= j){
        missing <- c(missing, max(observed+1):(diff(missingIntervals)[j] + max(observed)))
      }
      if(length(locationObservedIntervals) >= j){
        observed <- c(observed, (max(missing)+1):(diff(observedIntervals)[j] + max(missing)))
      }
    }

    sparseGrid[observed,i] <- data[observed,i]
  }

  colnames(sparseGrid) <- colnames(data)
  rownames(sparseGrid) <- rownames(data)
  return(sparseGrid)
}
