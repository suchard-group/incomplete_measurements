this.dir <- dirname(parent.frame(2)$ofile)

library(ggplot2)
library(grid)
library(gtable)
library(gridExtra)
library(RColorBrewer)
library(ggpubr)




my_theme <- function(){
  theme_bw(base_size=12) +
  theme(
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(hjust = 0.5),
      plot.title = element_text(hjust = 0.5)
    )
}



plot.labels <- c("HIV", "Mammals", "Prokaryotes", "N", "logMSE", "Sparsity", "Diffusion Correlation", "Residual Correlation", "Heritability", "Traits", "Bias")
names(plot.labels) <- c("hivSim", "mammalsSim", "prokSim", "nTaxa", "logmse", "sparsity", "diffCorr", "resCorr", "her", "traits", "bias")


simBoxPlot <- function(data, run, xVar, yVar, shadeVar, title="") {
  
  data.run = data
  if (run != "") {
    data.run <- data[which(data$run == run),]
  }
  
  print(title)
  
  p = ggplot() + geom_boxplot(aes(x = data.run[,xVar], y = data.run[,yVar], fill=data.run[,shadeVar]), lwd=.25, outlier.size=0.25) +
    colFill.sparsity + 
    ggtitle(title) +
    my_theme() + 
    labs(fill = plot.labels[shadeVar], x = plot.labels[xVar], y = plot.labels[yVar])
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
    plots[[i]] = simBoxPlot(datas[[i - n]], run, xVar, yVarTwo, colVar, title=plot.labels[titles[i - n]])
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




p <- simBoxPlot(mats.resCorr, "", "nTaxa", "logmse", "sparsity", "testTitle")
# p
# 
ps <- simLinePlot(mats.her.diag, "mammalsSim", "nObs", "logmse", "sparsity")
# ps
# 
pc <- simCoveragePlot(coverage.her.diag, "", "nObs", "coverage","run")
# pc

mat_dats <- list(mats.diffCorr, mats.resCorr, mats.her.diag, traits)
titles <- c("diffCorr", "resCorr", "her", "traits")
grid.hivMSE <- gridBox(mat_dats, "hivSim", "nTaxa", "logmse", "sparsity", titles)
grid.mammalsMSE <- gridBox(mat_dats, "mammalsSim", "nTaxa", "logmse", "sparsity", titles)
grid.prokMSE <- gridBox(mat_dats, "prokSim", "nTaxa", "logmse", "sparsity", titles)

plot.width <- 7.5
plot.height <- 5

ggsave("hivSimMSE.pdf", plot=grid.hivMSE, width=plot.width, height=plot.height, units="in")
ggsave("mammalsSimMSE.pdf", plot=grid.mammalsMSE, width=plot.width, height=plot.height, units="in")
ggsave("prokSimMSE.pdf", plot=grid.prokMSE, width=plot.width, height=plot.height, units="in")

grid.hivBias <- gridBox(mat_dats, "hivSim", "nObs", "bias", "sparsity", titles)
grid.mammalsBias <- gridBox(mat_dats, "mammalsSim", "nTaxa", "bias", "sparsity", titles)
grid.prokBias <- gridBox(mat_dats, "prokSim", "nTaxa", "bias", "sparsity", titles)

ggsave("hivSimBias.pdf", plot=grid.hivBias, width=plot.width, height=plot.height, units="in")
ggsave("mammalsSimBias.pdf", plot=grid.mammalsBias, width=plot.width, height=plot.height, units="in")
ggsave("prokSimBias.pdf", plot=grid.prokBias, width=plot.width, height=plot.height, units="in")


plot.height <- 10

grid.hivBoth <- gridBigBox(mat_dats, "hivSim", "nTaxa", "logmse", "bias", "sparsity", titles)
grid.mammalsBoth <- gridBigBox(mat_dats, "mammalsSim", "nTaxa", "logmse", "bias", "sparsity", titles)
grid.prokBoth <- gridBigBox(mat_dats, "prokSim", "nTaxa", "logmse", "bias", "sparsity", titles)

ggsave("hivSimBoth.pdf", plot=grid.hivBoth, width=plot.width, height=plot.height, units="in")
ggsave("mammalsSimBoth.pdf", plot=grid.mammalsBoth, width=plot.width, height=plot.height, units="in")
ggsave("prokSimBoth.pdf", plot=grid.prokBoth, width=plot.width, height=plot.height, units="in")


gc()
