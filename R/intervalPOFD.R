#' Random Interval Observability
#'
#' Generates samples of functions observed in different intervals. See Elías et al (2020).
#'
#' @param data functional data completely observed. pxn matrix being n the number of curves and p the number og evaluation points.
#' @param observability mean observed proportion of the domain where each function is observed.
#' @param ninterval if type = "interval", n_interval is an integer with the number of observed intervals 1, 2, 3...
#' Large values of this parameter requires a large parameter p to guarantee the observability level.
#' @param pIncomplete number between 0 and 1 related to the proportion of curves that suffers partially observability.
#' The default is 1 meaning that all the sample curves are partially observed.
#' @return a list containing two elements 1) a functional sample and 2) the same sample of functions but
#' partially observed following one of the schemes described in the argument type.
#'
#' @references Elías, Antonio, Jiménez, Raúl, Paganoni, Anna M. and Sangalli, Laura M. (2020). Integrated Depths for Partially Observed Functional Data.
#' @examples
#'
#' data <- sapply(1:100, function(x) runif(1)*sin(seq(0, 2*pi, length.out = 200)) +
#' runif(1)*cos(seq(0, 2*pi, length.out = 200)))
#'
#' data_pofd <- intervalPOFD(data, observability = 0.5, ninterval = 2, pIncomplete = 1)
#'
#' @export
intervalPOFD <- function(data, observability = NULL, ninterval = NULL, pIncomplete = NULL){

  P <- dim(data)[1]
  N <- dim(data)[2]

  time_grid <- seq(0, 1, length.out = P)

  if(pIncomplete == 1){
    podata <- nObservedAtRandom(data, observability, ninterval)

    rownames(podata) <- round(time_grid, digits = 5)

  }else{
    whichToCensor <- colnames(data)[sample(1:N, size = round(N*pIncomplete), replace = FALSE)]

    podata <- nObservedAtRandom(data[,whichToCensor], observability, ninterval)

    colnames(podata) <- whichToCensor
    dataCensoredMatrix <- cbind(data[,colnames(data)[!colnames(data) %in% whichToCensor]], podata)

    podata <- dataCensoredMatrix[,colnames(data)]
    rownames(podata) <- round(time_grid, digits = 5)
  }


  return(list(fd = data , pofd = podata))
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
