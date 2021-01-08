#' Outliergram for Partially Observed Functional Data
#'
#' Plots the Outliergram for PoFD and returns the shape outliers.
#'
#' @param data matrix p by n, being n the number of functions and p the number of grid points.
#' @param fshape inflation of the outliergram that determine the shape outlier rule.
#' @param p1 parameter of the outliergram for resampling method. Default = 1.
#' @param p2 parameter of the outliergram for resampling method. Default = 0.
#' @return a list with the functional outliergram for PoDF and the shape outliers.
#'
#' @import ggplot2
#' @importFrom magrittr "%>%"
#' @importFrom reshape2 melt
#'
#' @examples
#' data(exampleData)
#' outliergram_PoFD(exampleData$PoFDextremes_outliers, fshape = 1.5, p1 = 1, p2 = 0)
#'
#' @references Arribas-Gil, A. and Romo, J. (2014).  Shape outlier detection and visualization for functional data:  the outliergram.Biostatistics, 15(4):603–619.
#'
#' @export
outliergram_PoFD <- function(data, fshape = 1.5, p1 = 1, p2 = 0)
{
  N <- dim(data)[2]
  P <- dim(data)[1]

  if(is.null(rownames(data))){rownames(data) <- x <- c(1:P)}else{x <- as.numeric(rownames(data))}
  if(is.null(colnames(data))){colnames(data) <- ids <- c(1:N)}else{ids <- colnames(data)}

  mbd.ordered <- POIFD(data, type = "MBD")
  mbd.ids <- mbd.ordered[ids]
  epi.ids <- POIFD(data, type = "MHRD")$mepi[ids]

  #
  # Shape outlier detection Arribas-Gil and Romo J.
  Parabola <- 2*epi.ids*(1-epi.ids)
  Ps <- 2*seq(0,1,length.out = 100)*(1-seq(0,1,length.out = 100))
  d <- Parabola - mbd.ids
  limit <- fivenum(d)[4]+fshape*(fivenum(d)[4]-fivenum(d)[2])
  shape.out <- which(as.vector(d>=limit))

  #Level variation for further detection of shape outliers
  shape.out.further <- c()
  mbd.final <- mbd.ids
  epi.final <- epi.ids
  d.final <- d

  for (i in setdiff(1:N, shape.out)){
    data.aux <- data
    if (epi.ids[i] < 0.5){ #prop. mean curves above

      s <- sort(data[,i] - apply(data[, -i, drop = F], 1, function(x) quantile(x, na.rm = TRUE, probs = p1)))  #p1 percentile of each row compared to values in x[,i]

      if(max(s) > 0){
        data.aux[,i] <- data[,i] - max(s)

        mbd.ids.aux <- POIFD(data.aux, type = 'MBD')[ids]
        epi.ids.aux <- POIFD(data.aux, type = 'MHRD')$mepi[ids]

        if(mbd.ids.aux[i] < 2*epi.ids.aux[i]*(1-epi.ids.aux[i])-limit){ #if now it is an outlier
          shape.out.further <- append(shape.out.further, i)

          mbd.final[i] <- mbd.ids.aux[i]
          epi.final[i] <- epi.ids.aux[i]
          d.final[i] <- mbd.ids.aux[i]-2*epi.ids.aux[i]*(1-epi.ids.aux[i])
        }
      }
    }

    if (epi.ids[i] >= 0.5){

      s <- sort(data[,i] - apply(data[, -i, drop = F], 1, function(x) quantile(x, na.rm =TRUE, probs = p2)))

      if ( min(s) < 0){

        data.aux[,i] <- data[,i] - min(s)

        mbd.ids.aux <- POIFD(data.aux, type = 'MBD')[ids]
        epi.ids.aux <- POIFD(data.aux, type = 'MHRD')$mepi[ids]

        if(mbd.ids.aux[i] < 2*epi.ids.aux[i]*(1-epi.ids.aux[i])-limit){ #if now it is an outlier
          shape.out.further <- append(shape.out.further, i)

          mbd.final[i] <- mbd.ids.aux[i]
          epi.final[i] <- epi.ids.aux[i]
          d.final[i] <- mbd.ids.aux[i]-2*epi.ids.aux[i]*(1-epi.ids.aux[i])
        }
      }
    }
  }

  #plotting outliergram

  outgramParabolaOutliersPartially <- data.frame(id = rep(1, 100), x = seq(0,1,length.out = 100),  y = Ps)
  outgramParabolaEnhancedOutliersPartially <- data.frame(id = rep(1, 100), x = seq(0,1,length.out = 100), y = Ps - limit)
  outgramParabolaEnhancedOutliersPartially <- outgramParabolaEnhancedOutliersPartially[rowSums(outgramParabolaEnhancedOutliersPartially>0) ==3,]


  auxType <- rep("No-Outlier", N)
  auxType[shape.out.further] <- "Shape-Adjusted"
  auxType[shape.out] <- "Shape"

  outliergramDataOutliersPartially <- data.frame(
    id = c(1:N),
    group = factor(auxType),
    mepi = epi.final,
    mbd = mbd.final
  )

  outgramOutliersPartially <- outliergramDataOutliersPartially %>%
    ggplot2::ggplot(aes(x = .data$mepi, y = .data$mbd, group = .data$id, col = .data$group)) +
    ggplot2::geom_point(size = 1.5, aes(x = .data$mepi, y = .data$mbd, group = .data$id, color = .data$group)) +
    ggplot2::geom_line(data = outgramParabolaOutliersPartially, aes(.data$x, .data$y), color = "black")+
    ggplot2::geom_line(data = outgramParabolaEnhancedOutliersPartially, aes(.data$x, .data$y), color = "black", linetype = 2)+
    ggplot2::xlab("Partially Observed MEI") +
    ggplot2::ylab("Partially Observed MBD") +
    ggplot2::xlim(0,1)+
    ggplot2::ylim(0,0.5)+
    ggplot2::theme(legend.position = "bottom",
          legend.key = element_rect(colour = NA, fill = NA),
          legend.key.size = unit(0.75, "cm"),
          legend.title = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.line = element_line(colour = "black",
                                   size = rel(1))
    )

  return(list(outliergram = outgramOutliersPartially, shape = c(shape.out, shape.out.further)))

}
