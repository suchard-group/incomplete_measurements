library(PCMBase)
library(PCMBaseCpp)
library(ape)

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
data.dir <- file.path("..", "data")


PCMBaseTiming <- function(tree_path, data_path, Sigma_x_path, Sigmae_x_path, 
                          Nreps = 10000, 
                          fixed_tree = TRUE) {
  
  #Read tree and make PCMBase tree
  tree <- read.tree(tree_path)
  ptree <- PCMTree(tree)
  
  #Read data and format into matrix
  data <- read.csv(data_path)
  
  N <- nrow(data)
  P <- ncol(data) - 1
  taxa.data <- data$taxon
  taxa.tree <- tree$tip.label
  
  X <- matrix(NA, nrow=P, ncol=N) 
  
  # Reorder data
  for (i in 1:N) {
    ind <- match(taxa.tree[i], taxa.data)
    for (j in 1:P) {
      X[j, i] <- data[[j + 1]][ind]
    }
  }
  
  

  #Set PCMBase options to avoid errors from small eigenvalues / nearly singular matrices
  listOptions = list(PCMBase.Threshold.EV = 1e-10, PCMBase.Threshold.SV = 1e-10)
  do.call(options, listOptions)
  
  
  #Build model and 'mataI' objects from PCMBase
  model <- PCM(model="BM", k = P)
  
  
  metaI <- PCMInfo(
    X = X,
    tree = ptree,
    model = model)
  
  metaICpp <- PCMInfoCpp(X, ptree, model, metaI = metaI)
  
  
  #Load diffusion variance and residual variance matrices
  diff.var <- scan(Sigma_x_path)
  diff.var <- matrix(diff.var, nrow=P)
  model$Sigma_x[,,] <- rbind(diff.var)
  
  res.var <- matrix(scan(Sigmae_x_path), nrow = P)
  model$Sigmae_x[,,] <- rbind(res.var)
  
  #Likelihood computations
  t1 <- system.time(
    for (i in 1:Nreps) {
      if ((i %% 100) == 0) {
        print(sprintf("%d completed", i))
      }
      if (!fixed_tree) { #change tree branch lengths at random
        ptree$edge.length[sample(1:(2 * N - 2), 1)] = 100
        metaICpp <- PCMInfoCpp(X, ptree, model, metaI = metaI)
      }

      ll <- PCMLik(X = X, tree = ptree, model = model, metaI = metaICpp)
    }
  )
  
  t1
  
}