#' Integrated Depth for Partially Observed Functional Data
#'
#' Compute the depth measure of a partially observed functional data set proposed in [1].
#' If the functions are not observed in a common partially observed domain, the code first
#' estimates the observation domain using the proposal in [2].
#'
#' [1] Elías, A., Jiménez, R., Paganoni, A. M., & Sangalli, L. M. (2023).
#' Integrated depths for partially observed functional data. Journal of
#' Computational and Graphical Statistics, 32(2), 341-352.
#'
#' [2] Elías, A., Nagy, S. (2025).
#' Statistical Properties of Partially Observed Integrated Funcional Depths.
#' TEST, 34, 125-150.
#'
#' @param data If functions are observed in a partially observed common grid `data` is a matrix `p` by `n`, being `n` the number of functions and `p` the number of grid points.
#' The row names of the matrix should be the common evaluation grid and the column names the identifiers of each functional data.
#' If functions do not have a common grid, `data`must be a list of length `n` where each element contains the values and evaluation points of each function.
#' I.e. each list element must contain two vectors `x`, including the evaluation points, and `y`, the evaluated function values.
#' For functions without a common grid, the function \code{"POIFD"}) first apply the procedure to estimate the observation domain proposed in Elías, A., Nagy, S. (2025), TEST.
#'
#' @param type chosen depth measure. Halfspace depth (\code{"HS"}), Fraiman and Muniz depth (\code{"FMD"}), Modified band depth (\code{"MBD"})
#' or Modified Half Region Depth and Modified Epigraph/Hipograph Index \code{"MHRD"})
#'
#' @param phi phi function of weights for the POIFD. The default value is as in [1]: the proportion of observed functions at each time point.
#'
#' @param t If functions do not have a common grid, `t` represents the final common grid of evaluation points to apply the
#' procedure to estimate the observation domain proposed in [2].
#'
#' @return Ordered vector of depths. The names are the functions names (if provided) or the column position.
#'
#' @examples
#' data <- exampleData$PoFDintervals
#' poifd <- POIFD(data,  type = c("FMD"))
#'
#' data <- exampleData$PoFDintervals_list
#' poifd <- POIFD(data,  type = c("FMD"), t = seq(0, 1, length.out = 100))
#'
#' @importFrom igraph graph_from_adjacency_matrix components groups
#' @importFrom stats ecdf dist na.omit approx
#'
#' @export
#'

POIFD <- function(data, type = c("HD", "FMD", "MBD", "MHRD"), phi, t = NULL){

  if(is.matrix(data)){
    n <- dim(data)[2]

    if(is.null(colnames(data))){colnames(data) <- c(1:n)}

    data <- data[rowMeans(is.na(data)) != 1,] #No full NAs

  }

  # Estimation of O
  if(is.list(data)){

    n <- length(data)

    if(is.null(names(data))){names(data) <- c(1:n)}

    data <- estimate_O_interpolate(data, t, rho = function(x) 1/sqrt(length(x)))

    data <- data[,colMeans(is.na(data)) != 1]
  }


  depth <- switch(type,
                  HD = POIFD_HD(t(data), phi),
                  FMD = POIFD_FMD(t(data), phi),
                  MBD = POIFD_MBD(t(data), phi),
                  MHRD = POIFD_MHRD(t(data), phi))

  return(depth)
}

POIFD_FMD <- function(data, phi){

  O_ind <- !is.na(data)
  n <- nrow(data)
  d <- ncol(data)

  pFMD <- apply(data, 2, function(x) (1 - abs(2*ecdf(x)(x) - 1)))

  if(missing(phi)){
    # Elias et al, 2022, w = Q(t)/int(Q(s))

    q = apply(O_ind,2,mean) # estimated Q(t) = P(t in O)
    q_mat = t(matrix(q, nrow=d, ncol=n)) # matrix version of q
    q_mat[O_ind==FALSE] = NA

    wPO_pDepth <- q_mat*pFMD
    POIFD <- apply(wPO_pDepth,1,function(x) sum(x,na.rm=TRUE)/d)
    denoms <- apply(q_mat,1,function(x) sum(x,na.rm=TRUE)/d)
    POIFD <- POIFD/denoms
  }else{
    phi_mat <- matrix(rep(phi, n), ncol = d)
    wPO_pDepth <- phi_mat*pFMD
    nums <- apply(wPO_pDepth,1,function(x) sum(x,na.rm=TRUE)/d)
    denoms <- apply(phi_mat,1,function(x) sum(x,na.rm=TRUE)/d)
    POIFD <- nums/denoms
  }

  if(!is.null(rownames(data))){
    names(POIFD) <- rownames(data)
  }else{
    names(POIFD) <- 1:n
  }

  return(POIFD)
}

POIFD_HD <- function(data, phi){

  O_ind <- !is.na(data)
  n <- nrow(data)
  d <- ncol(data)

  eps <- min(diff(sort(data)))/10 # Exact computation
  pHD <- apply(data, 2, function(x) pmin(ecdf(x)(x),1-ecdf(x)(x-eps)))

  if(missing(phi)){
    # Elias et al, 2022, w = Q(t)/int(Q(s))

    q = apply(O_ind,2,mean) # estimated Q(t) = P(t in O)
    q_mat = t(matrix(q, nrow=d, ncol=n)) # matrix version of q
    q_mat[O_ind==FALSE] = NA

    wPO_pDepth <- q_mat*pHD
    POIFD <- apply(wPO_pDepth,1,function(x) sum(x,na.rm=TRUE)/d)
    denoms <- apply(q_mat,1,function(x) sum(x,na.rm=TRUE)/d)
    POIFD <- POIFD/denoms
  }else{
    phi_mat <- matrix(rep(phi, n), ncol = d)
    wPO_pDepth <- phi_mat*pHD
    nums <- apply(wPO_pDepth,1,function(x) sum(x,na.rm=TRUE)/d)
    denoms <- apply(phi_mat,1,function(x) sum(x,na.rm=TRUE)/d)
    POIFD <- nums/denoms
  }

  if(!is.null(rownames(data))){
    names(POIFD) <- rownames(data)
  }else{
    names(POIFD) <- 1:n
  }

  return(POIFD)
}

POIFD_MBD <- function(data, phi){

  O_ind <- !is.na(data)
  n <- nrow(data)
  d <- ncol(data)

  pMBD <- apply(data, 2, function(x) 2*ecdf(x)(x)*(1-ecdf(x)(x)))

  if(missing(phi)){
    # Elias et al, 2022, w = Q(t)/int(Q(s))

    q = apply(O_ind,2,mean) # estimated Q(t) = P(t in O)
    q_mat = t(matrix(q, nrow=d, ncol=n)) # matrix version of q
    q_mat[O_ind==FALSE] = NA

    wPO_pDepth <- q_mat*pMBD
    POIFD <- apply(wPO_pDepth,1,function(x) sum(x,na.rm=TRUE)/d)
    denoms <- apply(q_mat,1,function(x) sum(x,na.rm=TRUE)/d)
    POIFD <- POIFD/denoms
  }else{
    phi_mat <- matrix(rep(phi, n), ncol = d)
    wPO_pDepth <- phi_mat*pMBD
    nums <- apply(wPO_pDepth,1,function(x) sum(x,na.rm=TRUE)/d)
    denoms <- apply(phi_mat,1,function(x) sum(x,na.rm=TRUE)/d)
    POIFD <- nums/denoms
  }

  if(!is.null(rownames(data))){
    names(POIFD) <- rownames(data)
  }else{
    names(POIFD) <- 1:n
  }

  return(POIFD)
}

POIFD_MHRD <- function(data, phi){

  O_ind <- !is.na(data)
  n <- nrow(data)
  d <- ncol(data)

  PO_pRANK <- apply(data, 2,function(x)  rank(x, na.last = "keep" , ties.method = "max"))

  PO_pEPI <- apply(PO_pRANK, 2, function(x) (max(x, na.rm = TRUE) - x)/max(x, na.rm = TRUE)) #curves above at p

  PO_pHIPO <- apply(PO_pRANK, 2, function(x) (x - 1)/max(x, na.rm = TRUE)) #curves below at p

  if(missing(phi)){
    # Elias et al, 2022, w = Q(t)/int(Q(s))

    q = apply(O_ind,2,mean) # estimated Q(t) = P(t in O)
    q_mat = t(matrix(q, nrow=d, ncol=n)) # matrix version of q
    q_mat[O_ind==FALSE] = NA

    wPO_pEPI <- q_mat*PO_pEPI
    MEPI <- apply(wPO_pEPI, 1, function(x) sum(x,na.rm=TRUE)/d)
    denoms <- apply(q_mat, 1, function(x) sum(x,na.rm=TRUE)/d)
    MEPI <- (MEPI/denoms)

    wPO_pHIPO <- q_mat*PO_pHIPO
    MHIPO <- apply(wPO_pHIPO, 1, function(x) sum(x,na.rm=TRUE)/d)
    denoms <- apply(q_mat, 1, function(x) sum(x,na.rm=TRUE)/d)
    MHIPO <- (MHIPO/denoms)

  }else{

    phi_mat <- matrix(rep(phi, n), ncol = d)

    wPO_pEPI <- phi_mat*PO_pEPI
    MEPI <- apply(wPO_pEPI, 1, function(x) sum(x,na.rm=TRUE)/d)
    denoms <- apply(phi_mat,1, function(x) sum(x,na.rm=TRUE)/d)
    MEPI <- (MEPI/denoms)

    wPO_pHIPO <- phi_mat*PO_pHIPO
    MHIPO <- apply(wPO_pHIPO, 1, function(x) sum(x,na.rm=TRUE)/d)
    denoms <- apply(phi_mat,1, function(x) sum(x,na.rm=TRUE)/d)
    MHIPO <- (MHIPO/denoms)
  }

  MHRD <- apply(cbind(MEPI, MHIPO), 1, min)

  if(!is.null(rownames(data))){
    names(MHIPO) <- names(MEPI) <- names(MHRD) <- rownames(data)
  }else{
    names(MHIPO) <- names(MEPI) <- names(MHRD) <- 1:n
  }

  return(list(mepi = MEPI, mhipo = MHIPO, hrd = MHRD))
}

estimate_O_interpolate_one <- function(fo_list, x, rho){
  x_o <- fo_list$x
  f_o <- fo_list$y

  # get spheres and connection
  centers <- cbind(x_o)
  dist_centers <- as.matrix(dist(centers, method = "euclidean"))

  #adjm <- (dist_centers <= rho)*1
  adjm <- (dist_centers <= rho(na.omit(x_o)))*1

  diag(adjm) <- 0

  # get circles that intersect
  g1 <- graph_from_adjacency_matrix(adjm, mode = "undirected")
  clu <- components(g1)
  O_ind <- groups(clu)

  # get the observed intervals
  O_i <- list()
  for(i in 1:length(O_ind)){
    O_i[[i]] <- cbind(x_o[as.numeric(O_ind[[i]])], f_o[as.numeric(O_ind[[i]])])
  }

  # interpolation step:
  O_i_inter <- list()
  for(i in 1:length(O_ind)){
    if(nrow(O_i[[i]]) > 1){
      O_i_inter[[i]] <- approx(O_i[[i]][,1], O_i[[i]][,2], x)
    }else{ #asign the point to the closes x value
      aux_y <- rep(NA, length(x))
      aux_y[which.min(abs(x-O_i[[i]][,1]))] <- O_i[[i]][,2]

      O_i_inter[[i]] <- list(x = x, y = aux_y)
    }
  }

  # return a vector
  f_o_inter <- rep(NA, length = length(x))
  for(i in 1:length(O_i_inter)){
    f_o_inter[which(!is.na(O_i_inter[[i]]$y))] <- na.omit(O_i_inter[[i]]$y)
  }

  return(f_o_inter)
}

estimate_O_interpolate <- function(fo_list, x, rho){

  f_o_inter <- sapply(fo_list, function(y) estimate_O_interpolate_one(y, x, rho))

  return(t(f_o_inter))
}
