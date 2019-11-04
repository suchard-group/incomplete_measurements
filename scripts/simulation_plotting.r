this.dir <- dirname(parent.frame(2)$ofile)

library(ggplot2)

storage.dir <- file.path(this.dir, "storage")
traits.path <- file.path(storage.dir, "trait_simulation.csv")
matrix.path <- file.path(storage.dir, "matrix_simulation.csv")
coverage.path <- file.path(storage.dir, "coverage_simulation.csv")

traits <- read.csv(traits.path)
mats <- read.csv(matrix.path)
coverage <- read.csv(coverage.path)

mats <- mats[which(mats$mse != 0.0),] #remove rows with no error (i.e. the diagonals of a correlation matrix)

simBoxPlot <- function(data, run, xVar, yVar) {
  
  data.run <- data[which(data$run == run),]
  print(head(data.run))
  p <- ggplot() + geom_boxplot(aes(x = data.run[,xVar], y = data.run[,yVar]))

}
  


traits$sparsity <- factor(traits$sparsity)
traits$logmse <- unlist(lapply(traits$mse, log))
mats$logmse <- unlist(lapply(mats$mse, log))
mats$sparsity <- factor(mats$sparsity)

p <- simBoxPlot(mats, "hivSim", "sparsity", "logmse")
p
# mats.hiv <- mats[which(mats$run == "hivSim"),]
# 
# p <- ggplot() + geom_boxplot(aes(x = mats.hiv[,"sparsity"], y = mats.hiv[,"logmse"]))
# p
# 
# traits.ten <- traits[which(traits$sparsity == 0.1),]
# traits.fifty <- traits[which(traits$sparsity == 0.5),]

