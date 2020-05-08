library(ggplot2)
library(grid)

my_theme <- function(){
  theme_bw(base_size=12) %+replace%
    theme(
      
    )
}


make_correlation_plot <- function(data,
                                  labels,
                                  corrs.scale,
                                  corrs.min,
                                  text.size,
                                  use_lines = FALSE,
                                  custom_theme = my_theme
){
  corrs_theme <- function(){
    my_theme() %+replace%
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle=45, hjust = 1, vjust = 1),
        axis.text.y = element_text(hjust = 1),
        panel.border = element_blank()
      )
  }
  
  highlight.scale <- 0.9 * corrs.scale
  highlight.stroke <- 0.1 * corrs.scale
  p_scale <- 0.05
  
  p <- length(data$corrs)
  signs <- rep(NA, p)
  for(i in 1:p){
    if (!is.na(data$corrs[i])){
      if (data$corrs[i] > 0){
        signs[i] = "purple"
      }else{
        signs[i] = "orange"
      }
    }
  }
  
  data$signs <- signs
  
  print(data$signs)
  
  myplot <- ggplot(data, aes(row, col)) + 
    geom_point(aes(size=abs(data$corrs), col=data$corrs)) + 
    geom_point(shape = 1, aes(size = abs(data$corrs)), col = data$signs) +
    geom_point(shape = 8, aes(size = data$fill_p * p_scale)) +
    scale_size(range = c(corrs.min, 1) * corrs.scale, guide=FALSE) +
    scale_color_gradient2(low="orange", mid="white", high="purple", limits=c(-1, 1), name="Correlation") +
    # guides(fill=guide_legend(title="Correlation")) +
    geom_point(aes(size=corr_ref), col="grey") +
    geom_text(aes(label=data$string_p, size=text.size * corrs.scale, angle=45)) + 
    # annotate("text", x = 2, y = 2, label = "test") +
    # annotate("point", shape = 8, x = 3, y = 3) +
    coord_fixed() +
    corrs_theme()
    # annotation_custom(grob = pointsGrob(x = 0.35,
    #                                     y = .99, 
    #                                     pch= 8, size = unit(0.45, "char"),
    #                                     default.units = "native", name = NULL,
    #                                     gp = gpar(), vp = NULL))
    
  
  if (use_lines){
    for (coord in line_locations) {
      myplot <- myplot + geom_segment(x = coord, y = 0.5, xend = coord, yend = p + 0.5)
      myplot <- myplot + geom_segment(y = coord, x = 0.5, yend = coord, xend = p + 0.5)
    }
  }
  
  
 
  myplot
}

