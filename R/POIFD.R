#' Integrated Depth for Partially Observed Functional Data
#'
#' Compute the depth measures of a partially observed functional data set evaluated in a common grid.
#'
#' @param data matrix p by n, being n the number of functions and p the number of grid points.
#' Rownames are the dense grid x and colnames the identificator of each functional data.
#' @param type choosen depth measure. Fraiman and Muniz depth (\code{"FMD"}), Modified band depth (\code{"MBD"})
#' or Modified Half Region Depth and Modified Epigraph/Hipograph Index \code{"MHRD"})
#' @param phi phi function of weights for the POIFD. The default value is as in the paper, i.e. the proportion of observed functions
#' at each time point.
#'
#' @return Ordered vector of depths from the deepest to outward. The names are the functions names (if provided)
#' or the column position.
#'
#' @examples
#' data <- gaussian_PoFD(n=100, p=200, type="sparse", observability=0.5)$pofd
#' poifd <- POIFD(data,  type = c("MBD"))
#'
#' @export
POIFD <- function(data, type = c("MBD", "FMD", "MHRD"), phi){

  # if(!missing(xgrid)){
  #   data <- getBinnedData_PACE
  # }

  n <- dim(data)[2]
  if(is.null(colnames(data))){colnames(data) <- c(1:n)}

  depth <- switch(type, MBD = mbd(data, phi), FMD = fmd(data, phi), MHRD = mhrd(data, phi))

  return(depth)
}

mbd <- function(data, phi){
  n <- dim(data)[2]
  p <- dim(data)[1]

  if(missing("phi") ){
    phi_num <- (n - rowSums(is.na(data)))
  }else{
    phi_num <- phi*n
  }

  MBD_it <- t(mapply(aux.mbd, as.data.frame(t(data)), phi_num, USE.NAMES = FALSE))

  phi_den <- apply(data, 2, function(x) sum(phi_num[!is.na(x)]) )

  depth <- colSums(MBD_it, na.rm = TRUE)/phi_den
  names(depth) <- colnames(data)

  MBD <- sort(depth, decreasing=TRUE, index.return=TRUE); names(MBD$x) <- colnames(data)[MBD$ix]
  return(MBD$x)
}

aux.mbd <- function(data.t, phi_num){
  if(phi_num != 0){
    empiricalDistribution <- ecdf(data.t)(data.t)

    2*empiricalDistribution*(1 - empiricalDistribution)*phi_num
    }else{rep(NA, length(data.t))}
}

fmd <- function(data, phi){
  n <- dim(data)[2]
  p <- dim(data)[1]

  if( missing("phi") ){
    phi_num <- (n - rowSums(is.na(data)))
  }else{
    phi_num <- phi*n
  }

  FM_it <- t(mapply(aux.fmd, as.data.frame(t(data)), phi_num, USE.NAMES = FALSE))

  phi_den <- apply(data, 2, function(x) sum(phi_num[!is.na(x)]) )

  depth <- colSums(FM_it, na.rm = TRUE)/phi_den
  names(depth) <- colnames(data)

  FMd <- sort(depth, decreasing=TRUE, index.return=TRUE); names(FMd$x) <- colnames(data)[FMd$ix]
  return(FMd$x)
}

aux.fmd <- function(data.t, phi_num){
  if(phi_num != 0){(1 - abs(1/2 - ecdf(data.t)(data.t)))*phi_num}else{rep(NA, length(data.t))}
}

mhrd <- function(data, phi){
  n <- dim(data)[2]
  p <- dim(data)[1]

  if( missing("phi") ){
    phi_num <- (n - rowSums(is.na(data))) #observed functions at each t
  }else{
    phi_num <- phi*n
  }

  HRD_it <- mapply(aux.hrd, as.data.frame(t(data)), phi_num, USE.NAMES = FALSE)

  phi_den <- apply(data, 2, function(x) sum(phi_num[!is.na(x)]) )

  mhipo <- colSums(matrix(unlist(HRD_it["hipo",]), nrow = p, ncol = n, byrow = TRUE) , na.rm = TRUE)/phi_den
  mepi <- colSums(matrix(unlist(HRD_it["epi",]), nrow = p, ncol = n, byrow = TRUE) , na.rm = TRUE)/phi_den
  hrd <- apply(cbind(mepi, mhipo), 1, min)

  mepi_pofd <- sort(mepi, decreasing=TRUE, index.return=TRUE); names(mepi_pofd$x) <- colnames(data)[mepi_pofd$ix]
  mhipo_pofd <- sort(mhipo, decreasing=TRUE, index.return=TRUE); names(mhipo_pofd$x) <- colnames(data)[mhipo_pofd$ix]
  hrd_pofd <- sort(hrd, decreasing=TRUE, index.return=TRUE); names(hrd_pofd$x) <- colnames(data)[hrd_pofd$ix]

  return(list(mepi = mepi_pofd$x, mhipo = mhipo_pofd$x, hrd = hrd_pofd$x))
}

aux.hrd <- function(data.t, phi_num.t){
  if(phi_num.t != 0){
    rmat <- rank(data.t, na.last = "keep" , ties.method = "max")  #in ties, both curves get higher rank

    up <- phi_num.t - rmat
    down <- rmat - 1
    mepi <- up/phi_num.t
    mhipo <- down/phi_num.t

    output <- list(epi = mepi*phi_num.t, hipo = mhipo*phi_num.t)
  }else{
    output <- list(epi = rep(NA, length(data.t)), hipo = rep(NA, length(data.t)))
  }
  return(output)
}
