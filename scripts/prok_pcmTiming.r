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

source("pcm_timing.r")

prok_tree <- file.path(data.dir, "prokaryotes_newick.txt")
prok_data <- file.path(data.dir, "prokaryotes_processed_data.csv")
prok_diff_mat <- file.path("storage", "prok_diff_mat.mat")
prok_res_mat <- file.path("storage", "prok_res_mat.mat")

prok_time <- PCMBaseTiming(prok_tree, prok_data, prok_diff_mat, prok_res_mat, fixed_tree = TRUE, Nreps=1000)

setwd(this.dir)
# write(prok_time, file=file.path("..", "logs", "PCMBase_timing", "prokPCMComparisonTimer_pcm.txt"))
write(prok_time, file=file.path("prokPCMComparisonTimer_pcm.txt"))