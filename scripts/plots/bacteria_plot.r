library(ggplot2)
options(device="windows")
source("correlation_plot.r")

data.path <- "bacteriaCorrelation.csv"
labels.path <- "bacteriaCorrelationLabels.csv"

data <- read.csv(data.path)
labels <- read.csv(labels.path)

data$row <- factor(data$row, levels = labels$labels)
data$col <- factor(data$col, levels = rev(labels$labels))

corrs.abs <- abs(data$corrs)

corrs.min <- min(corrs.abs, na.rm=TRUE)
corrs.max <- max(corrs.abs, na.rm=TRUE)
corrs.scale <- 12.5
text.size <- 0.006
highlight.scale <- 0.9 * corrs.scale
highlight.stroke <- 0.1 * corrs.scale

p <- length(labels$labels)
line_locations <- seq(0.5, p + 0.5, by = 1)





myplot <- make_correlation_plot(data, labels, corrs.scale, corrs.min, text.size, use_lines=FALSE)


myplot <- myplot + annotation_custom(grob = textGrob(label = "indicates > 0.99", gp = gpar(fontsize=9)), 
                                         ymin = 7.5, xmin = 1, xmax = Inf, ymax = Inf) +
                    annotation_custom(grob = pointsGrob(x = 0.35,
                          y = .99,
                          pch= 8, size = unit(0.45, "char"),
                          default.units = "native", name = NULL,
                          gp = gpar(), vp = NULL))
                                      
  
gt <- ggplot_gtable(ggplot_build(myplot))
gt$layout$clip[gt$layout$name == "panel"] <- "off"


pdf("bacteriaCorrs.pdf", 
    height = 4.5, width = 5.5)
grid.draw(gt)
dev.off()
