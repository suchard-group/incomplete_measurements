this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(ggplot2)
library(grid)
library(gtable)
library(gridExtra)
library(RColorBrewer)
library(ggpubr)
library(magrittr)
library(dplyr)


my_theme <- function(){
  theme_bw(base_size=12) +
  theme(
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(hjust = 0.5),
      plot.title = element_text(hjust = 0.5)
    )
}

# modified from https://github.com/tidyverse/ggplot2/blob/master/R/stat-boxplot.r
stat_boxplot_custom <- function(mapping = NULL, data = NULL,
                                geom = "boxplot", position = "dodge",
                                ...,
                                qs = c(.025, .25, 0.5, 0.75, 0.975),
                                na.rm = FALSE,
                                show.legend = NA,
                                inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatBoxplotCustom,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      qs = qs,
      ...
    )
  )
}

StatBoxplotCustom <- ggproto("StatBoxplotCustom", Stat,
                             required_aes = c("x", "y"),
                             non_missing_aes = "weight",

                             setup_params = function(data, params) {
                               params$width <- ggplot2:::"%||%"(
                                 params$width, (resolution(data$x) * 0.75)
                               )

                               if (is.double(data$x) && !ggplot2:::has_groups(data) && any(data$x != data$x[1L])) {
                                 warning(
                                   "Continuous x aesthetic -- did you forget aes(group=...)?",
                                   call. = FALSE
                                 )
                               }

                               params
                             },

                             compute_group = function(data, scales, width = NULL, na.rm = FALSE, qs = c(.025, .25, 0.5, 0.75, 0.975)) {

                               if (!is.null(data$weight)) {
                                 mod <- quantreg::rq(y ~ 1, weights = weight, data = data, tau = qs)
                                 stats <- as.numeric(stats::coef(mod))
                               } else {
                                 stats <- as.numeric(stats::quantile(data$y, qs))
                               }
                               names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
                               iqr <- diff(stats[c(2, 4)])

                               outliers <- (data$y < stats[1]) | (data$y > stats[5])

                               if (length(unique(data$x)) > 1)
                                 width <- diff(range(data$x)) * 0.9

                               df <- as.data.frame(as.list(stats))
                               df$outliers <- list(data$y[outliers])

                               if (is.null(data$weight)) {
                                 n <- sum(!is.na(data$y))
                               } else {
                                 # Sum up weights for non-NA positions of y and weight
                                 n <- sum(data$weight[!is.na(data$y) & !is.na(data$weight)])
                               }

                               df$notchupper <- df$middle + 1.58 * iqr / sqrt(n)
                               df$notchlower <- df$middle - 1.58 * iqr / sqrt(n)

                               df$x <- if (is.factor(data$x)) data$x[1] else mean(range(data$x))
                               df$width <- width
                               df$relvarwidth <- sqrt(n)
                               df
                             }
)


plot.labels <- c("HIV", "Mammals", "Prokaryotes", "N", "logMSE", "Sparsity", "Diffusion Correlation", "Residual Correlation", "Heritability", "Traits", "Bias")
names(plot.labels) <- c("hivSim", "mammalsSim", "prokSim", "nTaxa", "logmse", "sparsity", "diffCorr", "resCorr", "her", "traits", "bias")


simBoxPlot <- function(data, run, xVar, yVar, shadeVar, title="", zero.line=FALSE) {
  data.run = data
  if (run != "") {
    data.run <- data[which(data$run == run),]
  }

  print(title)



  dodge.width <- 0.95

  p <- ggplot(data.run, aes(x = data.run[,xVar], y = data.run[,yVar], fill=data.run[,shadeVar]))
  if (zero.line) {
    p = p + geom_hline(yintercept=0, linetype="dashed", size=.5)
  }
  p = p +
    stat_boxplot_custom(lwd=0.25, outlier.size=0.75) +
    # stat_summary(fun.data=box_plot_quantiles, geom="boxplot", lwd=0.25, position=position_dodge(width=dodge.width)) +
    # geom_boxplot(outlier.size=0, lwd=0.75) +
    colFill.sparsity +
    # stat_summary(fun.y = box_plot_outliers, geom="point", size=1, position=position_jitterdodge(dodge.width=dodge.width), pch=21) +
    # geom_point(position=position_jitterdodge()) +
    # dummy.sparsity +
    # colCol.sparsity +
    ggtitle(title) +
    my_theme() +
    labs(fill = plot.labels[shadeVar], x = plot.labels[xVar], y = plot.labels[yVar])
}

box_plot_quantiles <- function(x) {
  r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

box_plot_outliers <- function(x) {
  r <- quantile(x, probs=c(0.025, 0.975))
  subset(x, x < r[1] | r[2] < x)
}

simLinePlot <- function(data, run, xVar, yVar, colVar, title="") {
  data.run = data
  if (run != "") {
    data.run <- data[which(data$run == run),]
  }

  p = ggplot() +
    # geom_point(aes(x = data.run[,xVar], y = data.run[,yVar], color=data.run[,colVar])) +
    stat_summary(aes(x = data.run[,xVar], y = data.run[,yVar], color=data.run[,colVar]), fun.data = "stat_sum_quantiles") +
    colFill.sparsity +
    ggtitle(title) +
    my_theme() +
    labs(color = plot.labels[colVar], x = plot.labels[xVar], y = plot.labels[yVar])

}

stat_sum_quantiles <- function(data) {
  y = mean(data)
  qs = quantile(data, c(0.05, 0.95))
  df = data.frame(y = y, ymin = qs[[1]], ymax = qs[[2]])
}

simCoveragePlot <- function(data, run, xVar, yVar, colVar) {
  data.run = data
  if (run != "") {
    data.run <- data[which(data$run == run),]
  }

  p = ggplot() +
    # geom_point(aes(x = data.run[,xVar], y = data.run[,yVar], color=data.run[,colVar])) +
    stat_summary(aes(x = data.run[,xVar], y = data.run[,yVar], color=data.run[,colVar]), fun.y = "mean", geom = "point") +
    labs(color = colVar, x = xVar, y = yVar) +
    ylim(0, 1)
}

gridBox <- function(datas, run, xVar, yVar, colVar, titles) {

  n = 4
  stopifnot(length(datas) == n)

  plots = vector("list", n + 1)

  for (i in 1:n) {
    plots[[i]] = simBoxPlot(datas[[i]], run, xVar, yVar, colVar, title=plot.labels[titles[i]])
  }
  plots[[1]] = plots[[1]] + theme(legend.direction = "horizontal")
  legend = get_legend(plots[[1]])
  for (i in 1:n) {
    plots[[i]] = plots[[i]] + theme(legend.position = "none")
  }

  plots[[n + 1]] = legend

  p = grid.arrange(grobs = plots,
                   nrow = 3, ncol= 2,
                   heights = c(1, 1, 0.1),
                   layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 5)),
                   top = textGrob(plot.labels[run], gp=gpar(fontsize=20)))

}

gridBigBox <- function(datas, run, xVar, yVarOne, yVarTwo, colVar, titles) {

  n = 4
  stopifnot(length(datas) == n)

  plots = vector("list", 2 * n + 1)

  for (i in 1:n) {
    plots[[i]] = simBoxPlot(datas[[i]], run, xVar, yVarOne, colVar, title=plot.labels[titles[i]])
  }

  for (i in (n + 1):(2 * n)) {
    plots[[i]] = simBoxPlot(datas[[i - n]], run, xVar, yVarTwo, colVar, title=plot.labels[titles[i - n]], zero.line = TRUE)
  }


  plots[[1]] = plots[[1]] + theme(legend.direction = "horizontal")
  legend = get_legend(plots[[1]])
  for (i in 1:(2 * n)) {
    plots[[i]] = plots[[i]] + theme(legend.position = "none")
  }

  plots[[2 * n + 1]] = legend

  p = grid.arrange(grobs = plots,
                   nrow = 5, ncol= 2,
                   heights = c(1, 1, 1, 1, 0.1),
                   layout_matrix = rbind(c(1, 5), c(2, 6), c(3, 7), c(4, 8), c(9, 9)),
                   top = textGrob(plot.labels[run], gp=gpar(fontsize=20)))

}

paperGrid <- function(data, runs, xVar, yVar, colVar, titles) {
  n = 3
  stopifnot(length(runs) == n)

  plots = vector("list", n + 1)

  for (i in 1:n) {
    plots[[i]] = simBoxPlot(data, runs[i], xVar, yVar, colVar, title=plot.labels[titles[i]])
  }
  # plots[[1]] = plots[[1]] + theme(legend.direction = "vertical")
  legend = get_legend(plots[[1]])
  for (i in 1:n) {
    plots[[i]] = plots[[i]] + theme(legend.position = "none")
  }

  plots[[n + 1]] = legend

  p = grid.arrange(grobs = plots,
                   nrow = 2, ncol= 2
                   # heights = c(1, 1, 0.1),
                   # layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 5)),
                   # top = textGrob(plot.labels[run], gp=gpar(fontsize=20))
                   )
}

gridLine <- function(datas, run, xVar, yVar, colVar, titles) {

  n = 4
  stopifnot(length(datas) == n)

  plots = vector("list", n + 1)

  for (i in 1:n) {
    plots[[i]] = simLinePlot(datas[[i]], run, xVar, yVar, colVar, title=plot.labels[titles[i]])
  }
  plots[[1]] = plots[[1]] + theme(legend.direction = "horizontal")
  legend = get_legend(plots[[1]])
  for (i in 1:n) {
    plots[[i]] = plots[[i]] + theme(legend.position = "none")
  }

  plots[[n + 1]] = legend

  p = grid.arrange(grobs = plots,
                   nrow = 3, ncol= 2,
                   heights = c(1, 1, 0.1),
                   layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 5)),
                   top = textGrob(plot.labels[run], gp=gpar(fontsize=20)))

}

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}




storage.dir <- file.path(this.dir, "storage")
traits.path <- file.path(storage.dir, "trait_simulation.csv")
matrix.path <- file.path(storage.dir, "matrix_simulation.csv")
coverage.path <- file.path(storage.dir, "coverage_simulation.csv")

traits <- read.csv(traits.path)
mats <- read.csv(matrix.path)
coverage <- read.csv(coverage.path)

traits <- traits[which(traits$isRandom == "true"),] #remove rows with no error (i.e. the diagonals of a correlation matrix)
mats <- mats[which(mats$isRandom == "true"),]
coverage <- coverage[which(coverage$isRandom == "true"),]




traits$sparsity <- factor(traits$sparsity)
mats$sparsity <- factor(mats$sparsity)
coverage$sparsity <- factor(coverage$sparsity)

stopifnot(levels(mats$sparsity)[1] == "0")
levels(mats$sparsity)[1] = "0.0"
stopifnot(levels(coverage$sparsity)[1] == "0")
levels(coverage$sparsity)[1] = "0.0"

colors.sparsity <- brewer.pal(length(levels(mats$sparsity)), "Set1")
names(colors.sparsity) <- levels(mats$sparsity)
colFill.sparsity <- scale_fill_manual(values=colors.sparsity)
colCol.sparsity <- scale_color_manual(values=colors.sparsity)
colors.dummy <- c("blue", "black", "black", "black")
dummy.sparsity <- scale_color_manual(values=colors.dummy)


traits$nTaxa <- factor(traits$nTaxa)
mats$nTaxa <- factor(mats$nTaxa)
coverage$nTaxa <- factor(coverage$nTaxa)

traits$logmse <- unlist(lapply(traits$mse, log))

mats$logmse <- unlist(lapply(mats$mse, log))

mats.diag <- mats[which(mats$component == "diagonal"),]
mats.offDiag <- mats[which(mats$component == "offDiagonal"),]

mats.diffCorr <- mats[which(mats$variable == "diffCorr"),]
mats.resCorr <- mats[which(mats$variable == "resCorr"),]
mats.her <- mats[which(mats$variable == "heritability"),]
mats.her.diag <- mats.her[which(mats.her$component == "diagonal"),]
mats.her.offdiag <- mats.her[which(mats.her$component == "offDiagonal"),]

coverage.diffCorr <- coverage[which(coverage$variable == "diffCorr"),]
coverage.resCorr <- coverage[which(coverage$variable == "resCorr"),]
coverage.her <- coverage[which(coverage$variable == "heritability"),]
coverage.her.diag <- coverage.her[which(coverage.her$component == "diagonal"),]
coverage.her.offdiag <- coverage.her[which(coverage.her$component == "offDiagonal"),]
coverage.trait <- coverage[which(coverage$variable == "trait"),]


#
#
# p <- simBoxPlot(mats.resCorr, "", "nTaxa", "logmse", "sparsity", "testTitle")
# # p
# #
# ps <- simLinePlot(mats.her.diag, "mammalsSim", "nObs", "logmse", "sparsity")
# # ps
# #
# pc <- simCoveragePlot(coverage.her.diag, "", "nObs", "coverage","run")
# # pc

mat_dats <- list(mats.diffCorr, mats.resCorr, mats.her.diag, traits)
titles <- c("diffCorr", "resCorr", "her", "traits")

# grid.hivMSE <- gridBox(mat_dats, "hivSim", "nTaxa", "logmse", "sparsity", titles)
# grid.mammalsMSE <- gridBox(mat_dats, "mammalsSim", "nTaxa", "logmse", "sparsity", titles)
# grid.prokMSE <- gridBox(mat_dats, "prokSim", "nTaxa", "logmse", "sparsity", titles)

plot.width <- 7.5
# plot.height <- 5
#
# ggsave("hivSimMSE.pdf", plot=grid.hivMSE, width=plot.width, height=plot.height, units="in")
# ggsave("mammalsSimMSE.pdf", plot=grid.mammalsMSE, width=plot.width, height=plot.height, units="in")
# ggsave("prokSimMSE.pdf", plot=grid.prokMSE, width=plot.width, height=plot.height, units="in")
#
# grid.hivBias <- gridBox(mat_dats, "hivSim", "nObs", "bias", "sparsity", titles)
# grid.mammalsBias <- gridBox(mat_dats, "mammalsSim", "nTaxa", "bias", "sparsity", titles)
# grid.prokBias <- gridBox(mat_dats, "prokSim", "nTaxa", "bias", "sparsity", titles)
#
# ggsave("hivSimBias.pdf", plot=grid.hivBias, width=plot.width, height=plot.height, units="in")
# ggsave("mammalsSimBias.pdf", plot=grid.mammalsBias, width=plot.width, height=plot.height, units="in")
# ggsave("prokSimBias.pdf", plot=grid.prokBias, width=plot.width, height=plot.height, units="in")


plot.height <- 8.5

grid.hivBoth <- gridBigBox(mat_dats, "hivSim", "nTaxa", "logmse", "bias", "sparsity", titles)
grid.mammalsBoth <- gridBigBox(mat_dats, "mammalsSim", "nTaxa", "logmse", "bias", "sparsity", titles)
grid.prokBoth <- gridBigBox(mat_dats, "prokSim", "nTaxa", "logmse", "bias", "sparsity", titles)

ggsave("hivSimBoth.pdf", plot=grid.hivBoth, width=plot.width, height=plot.height, units="in")
ggsave("mammalsSimBoth.pdf", plot=grid.mammalsBoth, width=plot.width, height=plot.height, units="in")
ggsave("prokSimBoth.pdf", plot=grid.prokBoth, width=plot.width, height=plot.height, units="in")

runs <- c("mammalsSim", "prokSim", "hivSim")
grid.paper <- paperGrid(mats.diffCorr, runs, "nTaxa", "logmse", "sparsity", runs)


plot.height <- 5
ggsave("diffCorSim.pdf", plot=grid.paper, width=plot.width, height=plot.height, units="in")


gc()
