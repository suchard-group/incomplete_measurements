library(ggplot2)
options(device="windows")
source("correlation_plot.r")


data.path <- file.path("mammalsCorrelation.csv")
labels.path <- file.path("mammalsCorrelationLabels.csv")

data <- read.csv(data.path)
labels <- read.csv(labels.path)

data$row <- factor(data$row, levels = labels$labels)
data$col <- factor(data$col, levels = rev(labels$labels))

corrs.abs <- abs(data$corrs)

corrs.min <- min(corrs.abs, na.rm=TRUE)
corrs.max <- max(corrs.abs, na.rm=TRUE)
corrs.scale <- 12.5
text.size <- 0.007
highlight.scale <- 0.8 * corrs.scale
highlight.stroke <- 0.1 * corrs.scale

p <- length(labels$labels)
line_locations <- seq(0.5, p + 0.5, by = 1)

highlighted_index = c(7, 8)
highlight_color = "red"
highlight_width = 3



myplot <- make_correlation_plot(data, labels, corrs.scale, corrs.min, text.size, use_lines=FALSE)




# Code to override clipping

myplot <- myplot + annotation_custom(grob = textGrob(label = "indicates > 0.99", gp = gpar(fontsize=9)), 
                                ymin = 8.5, xmin = 1, xmax = Inf, ymax = Inf) +
          annotation_custom(grob = pointsGrob(x = 0.37,
                                    y = .993, 
                                    pch= 8, size = unit(0.45, "char"),
                                    default.units = "native", name = NULL,
                                    gp = gpar(), vp = NULL))
  
gt <- ggplot_gtable(ggplot_build(myplot))
gt$layout$clip[gt$layout$name == "panel"] <- "off"

pdf("mammalsCorrs.pdf", width = 6, height = 4.8)
grid.draw(gt)
dev.off()