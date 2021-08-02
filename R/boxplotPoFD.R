#' Functional Boxplot for Partially Observed Functional Data
#'
#' Plots the Functional Boxplot for PoFD and returns the magnitude and domain outliers.
#' Magnitude outliers in blue, a dotted red indicates that the outlier situation occurs
#' in a region with less than \code{fdom} proportion of the central region.
#'
#' @param data matrix p by n, being n the number of functions and p the number of grid points.
#' @param centralRegion number between 0 and 1 determining the proportion of the deepest functions that builds the central region.
#' @param fmag factor to enhance the functional central region and
#'  determine the functional whiskers. Default is equal to 1.5. The whiskers provide the rule to unmask magnitude outliers.
#' @param fdom factor that provides the maximum proportion of observed functions in the central region to consider a magnitude outlier
#' as a domain outlier also. A value equals to 0 means that domain outliers are those functions that are observed
#' on the domain where any of the functions building the central region are observed.
#' A value equals to 1 determine as domain outlier any magnitude outlier out of the region where the central region
#' is completely observed.
#' @return a list with the functional boxplot for PoDF the magnitude outliers and the domain outliers.
#'
#' @examples
#' data(exampleData)
#' boxplotPOFD(exampleData$PoFDextremes_outliers, centralRegion = 0.5, fmag = 1.5, fdom = 0)
#'
#' @references Sun, Y. and Genton, M. G. (2011). Functional boxplots. Journal of Computational & Graphical Statistics, 20(2):316–334.
#' @import ggplot2
#' @importFrom magrittr "%>%"
#' @importFrom reshape2 melt
#' @importFrom stats ecdf fivenum median quantile
#'
#' @export
boxplotPOFD <- function(data, centralRegion = 0.5, fmag = 1.5, fdom = 0)
{
  N <- dim(data)[2]
  P <- dim(data)[1]
  w <- (N-rowSums(is.na(data)))/(N)

  if(is.null(rownames(data))){rownames(data) <- x <- c(1:P)}else{x <- as.numeric(rownames(data))}
  if(is.null(colnames(data))){colnames(data) <- ids <- c(1:N)}else{ids <- colnames(data)}

  mbd.ordered <- POIFD(data, type = "MBD")
  mbd.ids <- mbd.ordered[ids]
  median <- names(mbd.ordered)[1]

  # Magnitude outliers using PoFD MBD
  idCR <- names(mbd.ordered[1:ceiling(N*centralRegion)]) # bagId
  center <- data[,idCR]
  propCR <- rowMeans(!is.na(center))

  infCR <- apply(center, 1,function(x) suppressWarnings(min(x, na.rm=TRUE))); infCR[infCR == Inf | infCR == -Inf] = NA
  supCR <- apply(center, 1,function(x) suppressWarnings(max(x, na.rm=TRUE))); supCR[supCR == Inf | supCR == -Inf] = NA
  upperWhisker <- supCR + fmag*(supCR-infCR)
  lowerWhisker <- infCR - fmag*(supCR-infCR)

  #outliers rule if fdom = 1 dom.out = mag.out if fdom = 0 dom.out only where there is not CR
  #mag.out includes all the outliers dom.out only the ones determine by a region poorly observed

  if(fdom < 1){
  dom.out <- which(colSums(data <= lowerWhisker & propCR <= fdom |
                           data >= upperWhisker & propCR <= fdom, na.rm = TRUE) != 0)
  }else{
  dom.out <- which(colSums(data <= lowerWhisker & propCR < fdom |
                               data >= upperWhisker & propCR < fdom, na.rm = TRUE) != 0)
  }

  mag.out <- which(colSums( data <= lowerWhisker | data >= upperWhisker , na.rm = TRUE) != 0)

  #plotting
  fdata <- data.frame(id = rep(ids, each = P),
                  x = rep(x, N),
                  y = c(data)
                  )

  fdataMag <- data.frame(id = rep(mag.out, each = P),
                           x = rep(x, length(mag.out)),
                           y = c(data[,mag.out])
  )

  fdataDom <- data.frame(id = rep(dom.out, each = P),
                     x = rep(x, length(dom.out)),
                     y = c(data[,dom.out])
  )

  dataWinkler <- data.frame(id = c(rep(1, P), rep(2, P)),
                        x = c(x, rev(x)),
                        y = c(upperWhisker, rev(lowerWhisker)),
                        wcol = c(w, rev(w)))

  dataWinklerBars <- data.frame(
    colBar = c(w[ceiling(P/2)], w[ceiling(P/2)]),
    x1 = c(x[ceiling(P/2)], x[ceiling(P/2)]),
    x2 = c(x[ceiling(P/2)], x[ceiling(P/2)]),
    y1 = c(supCR[ceiling(P/2)], infCR[ceiling(P/2)]),
    y2 = c(upperWhisker[ceiling(P/2)], lowerWhisker[ceiling(P/2)])
  )

  auxPolygonX <- c(matrix(rep(c(1,2,2,1), P-1), nrow = 4) + matrix(rep(seq(0:c(P-2))-1, each = 4), nrow = 4))
  auxPolygonLowHigh <- c(matrix(rep(c(1,2,P+2,P+1), P-1), nrow = 4) + matrix(rep(seq(0:c(P-2))-1, each = 4), nrow = 4))

  dataBand <- data.frame(
    id = rep(factor(1:c(P-1)), each = 4),
    x = x[auxPolygonX],
    y = c(infCR, supCR)[auxPolygonLowHigh],
    wcol = rep(w[2:P], each = 4)
  )

  ggplot2::ggplot() +
    ggplot2::geom_polygon(data = dataBand, aes(x = .data$x, y = .data$y, group = .data$id, fill = .data$wcol, color = .data$wcol)) +
    ggplot2::geom_line(data = fdata[fdata$id == median,], aes(x = .data$x, y = .data$y), color = 'yellow', cex = 0.75) +
    ggplot2::geom_line(data = dataWinkler, aes(x = .data$x, y = .data$y, group = .data$id, color = .data$wcol), cex = 1) +
    ggplot2::geom_segment(data = dataWinklerBars, aes(x = .data$x1, y = .data$y1, xend = .data$x2, yend = .data$y2, color =.data$colBar), size = 0.75) +
    ggplot2::geom_line(data = fdataMag, aes(x = .data$x, y = .data$y, group = .data$id), color='blue', cex = 0.75)+
    ggplot2::geom_line(data = fdataDom, aes(x = .data$x, y = .data$y, group = .data$id), color='red', linetype = "dashed", cex = 0.75)+
    ggplot2::scale_fill_gradient(low = "white", high = "black", limits=c(0, 1) , aesthetics = c('fill', "color")) +
    ggplot2::labs(x = "", y ="")+
    ggplot2::theme(legend.position = "none",
          legend.key.size = unit(0.75, "cm"),
          legend.title = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.line = element_line(colour = "black",
                                   size = rel(1))
    ) -> boxplotPartially

  return(list(fboxplot = boxplotPartially, all.out = unique(c(mag.out, dom.out)), magnitude = mag.out, domain = dom.out))
}
