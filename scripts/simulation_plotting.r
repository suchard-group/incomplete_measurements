this.dir <- dirname(parent.frame(2)$ofile)

library(ggplot2)

simBoxPlot <- function(data, run, xVar, yVar, shadeVar) {
  
  data.run <- data[which(data$run == run),]
  p <- ggplot() + geom_boxplot(aes(x = data.run[,xVar], y = data.run[,yVar], fill=data.run[,shadeVar], group=data.run[,xVar])) +
    labs(fill = shadeVar, x = xVar, y = yVar) + scale_x_log10()
}

simScatterPlot <- function(data, run, xVar, yVar, colVar) {
  data.run <- data[which(data$run == run),]
  p <- ggplot() + geom_point(aes(x = data.run[,xVar], y = data.run[,yVar], color=data.run[,colVar]), position="jitter") +
    labs(color = colVar, x = xVar, y = yVar)
  
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


traits$nTaxa <- factor(traits$nTaxa)
mats$nTaxa <- factor(mats$nTaxa)
coverage$nTaxa <- factor(coverage$nTaxa)

# traits$nObs <- factor(traits$nObs)
# mats$nObs <- factor(mats$nObs)
# coverage$nObs <- factor(coverage$nObs)

traits$logmse <- unlist(lapply(traits$mse, log))

mats$logmse <- unlist(lapply(mats$mse, log))

mats.diag <- mats[which(mats$component == "diagonal"),]
mats.offDiag <- mats[which(mats$component == "offDiagonal"),]


p <- simBoxPlot(mats.diag, "mammalsSim", "nObs", "logmse", "sparsity")
p
# ps <- simScatterPlot(mats.offDiag, "mammalsSim", "nObs", "logmse", "sparsity")
# 
# ps