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

setwd(this.dir)

source("PhylogeneticEM_timing.r")

tree <- file.path(data.dir, "mammals_trimmed_newick.txt")
data <- file.path(data.dir, "mammals_log_data.csv")
diff_mat <- file.path("storage", "mammals_diff_mat.mat")

time <- PhylogeneticEMTiming(tree, data, diff_mat, fixed_tree = TRUE, Nreps=1000)

setwd(this.dir)

write(time, file=file.path("mammalsPCMComparisonTimer_phy.txt"))