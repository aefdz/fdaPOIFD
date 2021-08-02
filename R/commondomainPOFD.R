#' Common Domain Observability
#'
#' Generates samples of functions observed in a common domain in the center part of the domain. See Elías et al (2020).
#'
#' @param data functional data completely observed. pxn matrix being n the number of curves and p the number og evaluation points.
#' @param observability mean observed proportion of the domain where each function is observed.
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
#' data_pofd <- commondomainPOFD(data, observability = 0.5, pIncomplete = 1)
#'
#' @export
commondomainPOFD <- function(data, observability = NULL, pIncomplete = NULL){

  P <- dim(data)[1]
  N <- dim(data)[2]

  time_grid <- seq(0, 1, length.out = P)

  if(pIncomplete == 1){
    podata <- extremesCensoring(data, observability)

    rownames(podata) <- round(time_grid, digits = 5)

  }else{
    whichToCensor <- colnames(data)[sample(1:N, size = round(N*pIncomplete), replace = FALSE)]

    podata <- extremesCensoring(data[,whichToCensor], observability)

    colnames(podata) <- whichToCensor
    dataCensoredMatrix <- cbind(data[,colnames(data)[!colnames(data) %in% whichToCensor]], podata)

    podata <- dataCensoredMatrix[,colnames(data)]
    rownames(podata) <- round(time_grid, digits = 5)
  }


  return(list(fd = data , pofd = podata))
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

