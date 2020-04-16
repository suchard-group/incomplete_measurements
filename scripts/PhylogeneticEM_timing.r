# this.dir <- dirname(parent.frame(2)$ofile)
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}
this.dir <- dirname(thisFile())

data.dir <- file.path("..", "data")



library(PhylogeneticEM)
library(ape)


PhylogeneticEMTiming <- function(tree_path, data_path, diffVar_path, resVar_path, 
                          Nreps = 10000, 
                          fixed_tree = TRUE) {
  
  #Read tree and make PCMBase tree
  tree <- read.tree(tree_path)

  
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
  
  
  
 
  
  
  
  #Load diffusion variance and residual variance matrices
  diff.var <- scan(diffVar_path)
  diff.var <- matrix(diff.var, nrow=P)

  # res.var <- matrix(scan(resVar_path), nrow = P)
  

  #Build model
  params <- params_BM(p = P, variance = diff.var, phylo = phy)
  
  phy <- reorder(tree, order = "postorder")
  Beta <- matrix(0, ncol = length(phy$edge.length), nrow = P)
  params$root.state$var.root <- as.matrix(params$root.state$var.root)
  params$variance <- as.matrix(params$variance) 
  ll_bare_tree_param <- function(data, phy, Beta, params) {
    PhylogeneticEM:::log_likelihood_BM(data, phy$edge, Beta,
                                       params$variance, phy$edge.length,
                                       params$root.state) }
  
  

  #Likelihood computations
  t1 <- system.time(
    for (i in 1:Nreps) {
      if ((i %% 100) == 0) {
        print(sprintf("%d completed", i))
      }
      if (!fixed_tree) { #change tree branch lengths at random
        stop("not implemented")
      }
      
      ll_bare_tree_param(X, phy, Beta, params)
    }
  )
  
  t1
  
}
