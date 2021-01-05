#' Plot Partially Observed Functional Data
#'
#' Plot the sample of partially observed curves and the proportion of observed functions.
#'
#' @param data matrix p by n, being n the number of functions and p the number of grid points.

#' @return Plot of the partially observed functinal data and the proportion of observed functions
#' at each time point.
#'
#' @examples
#' #plot_PoFD(data)
#'
#' @import ggplot2
#' @import tibble
#' @importFrom magrittr "%>%"
#' @import patchwork
#' @importFrom reshape2 melt
#'
#' @export
#'
plot_PoFD <- function(data){

  N <- dim(data)[2]
  P <- dim(data)[1]

  if(is.null(rownames(data))){rownames(data) <- x <- c(1:P)}else{x <- as.numeric(rownames(data))}
  if(is.null(colnames(data))){colnames(data) <- c(1:N)}

  plotPoFD <- reshape2::melt(data, id='x') %>%
    ggplot2::ggplot(aes(x = .data$Var1, y = .data$value, col = as.factor(.data$Var2), group = as.factor(.data$Var2))) +
    ggplot2::geom_line(col = "black", size = 0.25, alpha = 1) +
    ggplot2::theme(legend.position = "none",
              legend.title = element_blank(),
              panel.background = element_blank(),
              plot.title = element_text(hjust = 0.5),
              axis.title =element_blank(),
              axis.line = element_line(colour = "black",
                                     size = rel(1)))+
    ggplot2::ggtitle("Partially Observed Funcional Data")


  w <- (N-rowSums(is.na(data)))/N
  w_dataFrame <- tibble(
    X1 = round(w, digits = 2),
    X2 = x,
    X3 =  c(x[2:P], x[P])
  )

  plotPropPoFD <- ggplot2::ggplot() +
    ggplot2::geom_rect(w_dataFrame, mapping = aes(xmin = .data$X2, xmax = .data$X3, ymin = 0, ymax = 1, fill = .data$X1, color = .data$X1)) +
    ggplot2::geom_line(data = w_dataFrame, mapping = aes(x = .data$X2, y = .data$X1), col= 'black') +
    ggplot2::scale_fill_gradient(low = "white", high = "black", limits = c(0, 1), aesthetics = c("color", "fill")) +
    ggplot2::labs(y=expression('q'[n]), x= "") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0)), breaks = c(0, 0.5, 1) , labels=c("0", "0.5", "1"))+
    ggplot2::theme(legend.position = "none",
                          legend.title = element_blank(),
                          panel.background = element_blank(),
                          plot.title = element_text(hjust = 0.5),
                          axis.line = element_line(colour = "black",
                                                   size = rel(1))
                    ) +
    ggplot2::ggtitle("Proportion of observed functional data")


  # Annotate if there is a common domain
  if(sum(w == 1) >0){
    plotPropPoFD <- plotPropPoFD +
      ggplot2::annotate("segment", x = x[which(w==1)[1]], xend = x[which(w==1)[sum(w==1)]], y = 0.5, yend = 0.5, colour = "white", size = 0.25, alpha=1, arrow=arrow(ends = "both", angle = 5, type = "closed"))+
      ggplot2::annotate(geom="text", x = (x[which(w==1)[sum(w==1)]] - x[which(w==1)[1]])/2 + x[which(w==1)[1]], y = 0.7, label="Common Domain", color="white", size = 2.5)
  }


  finalPlot <- plotPoFD / plotPropPoFD + patchwork::plot_layout(heights = c(0.85, 0.15))

  return(finalPlot)
}
