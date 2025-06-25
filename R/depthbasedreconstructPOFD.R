#' Depth-based reconstruction of partially observed functional data
#'
#' This function implements the reconstruction procedure [1]
#' which is based on the depth measure [2] for partially observed functional
#' data. Missing trajectories are imputed by the mean of the k nearest
#' neighbors within the envelope. The parameter k is tuned minimizing the Mean Squared Error of the
#' reconstruction in the observed part of the curve.
#'
#' [1] Elías, A., Jiménez, R., & Shang, H. L. (2023). Depth-based reconstruction
#' method for incomplete functional data. Computational Statistics, 38(3),
#' 1507-1535.
#'
#' [2] Elías, A., Jiménez, R., Paganoni, A. M., & Sangalli, L. M. (2023).
#' Integrated depths for partially observed functional data. Journal of
#' Computational and Graphical Statistics, 32(2), 341-352.
#'
#' @param data Data matrix `p` by `n`, being `n` the number of functions and `p` the number of grid points.
#' The row names of the matrix should be the common evaluation grid and the column names the identifiers of each functional data.
#' @param id_recons Vector indicating functions to be reconstructed. By default, all functions are reconstructed.
#'
#' @return The reconstructed data matrix 'recons_data'.
#'
#' @examples
#' data <- exampleData$PoFDintervals
#' recons_data <- depthbasedreconstructionPOFD(data, id_recons = 1:2)
#'
#' @export

depthbasedreconstructionPOFD <- function(data, id_recons = 1:dim(data)[2]) {
  p <- dim(data)[1]
  n <- dim(data)[2]

  X.hat <- matrix(NA, nrow = p, ncol = length(id_recons))

  for(t in id_recons) {
    J <- envelope(t(data), t)

    # Determine k and reconstruct the curve using the mean of the k nearest
    # neighbors within the envelope.
    X.hat[,t] <- rowMeans(data[,J[1:learningK(t(data), t, J)], drop = FALSE], na.rm=TRUE)

    X.hat[!is.na(data[,t]),t] <- data[!is.na(data[,t]),t]
  }

  return(X.hat)
}

#' Envelope algorithm
#'
#' Code for obtaining the envelope Ji of a curve specified by the index i. The
#' implementation is based on Algorithm 1 in [1]. References:
#'
#' [1] Elías, A., Jiménez, R., & Shang, H. L. (2023). Depth-based reconstruction
#' method for incomplete functional data. Computational Statistics, 38(3),
#' 1507-1535.
#'
#' @param data    Data matrix.
#' @param i        Index of curve.
#' @param max_iter Maximum number of nearest curves considered in for loop. By
#'                 default, max_iter = 5.
#'
#' @return Envelope (set of indices).
#' @export

envelope <- function(data, i, max_iter = 10) {
  T <- dim(data)[1]
  N <- dim(data)[2]

  # LINE 1: Initialization
  f <- data[i,]
  obs <- !is.na(f)
  J <- c()

  indices <- 1:T
  rownames(data) <- indices
  indices <- indices[-i]

  depth <- 0

  while(length(indices) > 1) {
    # LINE 2: Find nearest curve y' to f
    dist <- numeric(length(indices))

    # Compute distances (P3)
    for(k in 1:length(indices)) {
      obsk <- which(!is.na(data[indices[k],]) & !is.na(data[i,]))
      if(length(obsk) == 0) {
        dist[k] <- Inf
      } else {
        dist[k] <- sqrt(sum((f[obsk] - data[indices[k],obsk])^2)/N)/(length(obsk)/N)
      }
    }

    nearest <- order(dist)
    Nset <- c(nearest[1])

    # LINE 3: Iterate from nearest curve to farthest
    for(k in 2:min(length(indices), max_iter)) {
      j <- nearest[k]
      y <- data[indices[j],]
      Jstar <- c(J, indices[Nset])

      # LINE 6: If f is more enveloped by union(N, y) than by N (P2)
      lNy <- which(apply(data[indices[c(Nset, j)],obs], 2, min) <= f[obs] & apply(data[indices[c(Nset, j)],obs], 2, max) >= f[obs])
      lN  <- which(apply(matrix(data[indices[c(Nset)],obs], nrow = length(Nset)), 2, min) <= f[obs]
                   & apply(matrix(data[indices[c(Nset)],obs], nrow = length(Nset)), 2, max) >= f[obs])
      more_enveloped <- length(lNy) > length(lN)

      # LINE 6: If y enlarges the cover
      obsj <- !is.na(data[indices[j],])
      obsunion <- colSums(!is.na(matrix(data[Jstar,], nrow = length(Jstar), ncol = N))) > 0
      more_covered <- length(which((obsj) - sum(obsj[obsunion]) == 1)) > 0

      if(more_enveloped | more_covered) {
        Nset <- union(Nset, j)
      }
    }

    # LINE 10: Check depth (P1)
    Yaux <- t(matrix(data[c(indices[Nset], J, i),], ncol = N))
    colnames(Yaux) <- c(indices[Nset], J, i)
    poifd <- POIFD(Yaux, type = c("MBD"))[as.character(i)]
    if(poifd > depth) {
      depth <- poifd

      # LINE 11: Update J
      J <- union(J, indices[Nset])
    }

    indices <- indices[-Nset]
  }

  # Sort J according to depth
  Yaux <- t(matrix(data[J,], ncol = N))
  colnames(Yaux) <- J
  J <- as.integer(names(POIFD(Yaux, type = c("MBD"))))

  return(J)
}

#' Learning k
#'
#' This function chooses an optimal parameter k (denoting the number of curves
#' within the envelope that are considered for reconstruction).
#'
#' @param data Data matrix.
#' @param i     Index of curve.
#' @param J     Envelope.
#'
#' @return Optimal parameter k.
#' @export

learningK <- function(data, i, J){
  N <- dim(data)[2]

  observed <- rep(FALSE, N)
  sum_observed <- c()

  error <- c()

  for(k in 1:length(J)){
    observed[!is.na(data[J[k],])] <- TRUE
    sum_observed[k] <- sum(observed)
    reconstructedK <- colMeans(data[J[1:k], , drop = FALSE], na.rm=TRUE)
    error[k] <- sqrt(sum((data[i,!is.na(data[i,])] - reconstructedK[!is.na(data[i,])])^2,
                         na.rm = TRUE))
  }

  error[which(sum_observed < max(sum_observed))] <- Inf
  kopt <- which.min(error)

  return(kopt)
}
