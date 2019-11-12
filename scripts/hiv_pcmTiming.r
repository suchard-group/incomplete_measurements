# this.dir <- dirname(parent.frame(2)$ofile)
# this.dir <- getSrcDirectory(function(x) {x})
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

print(this.dir)

setwd(this.dir)

source("pcm_timing.r")

hiv_tree <- file.path(data.dir, "hiv_newick.txt")
hiv_data <- file.path(data.dir, "hiv_processed_data.csv")
hiv_diff_mat <- file.path("storage", "hiv_diff_mat.mat")
hiv_res_mat <- file.path("storage", "hiv_res_mat.mat")

hiv_time <- PCMBaseTiming(hiv_tree, hiv_data, hiv_diff_mat, hiv_res_mat, fixed_tree = TRUE, Nreps=1000)

setwd(this.dir)
# write(hiv_time, file=file.path("..", "logs", "PCMBase_timing", "hivPCMComparisonTimer_pcm.txt"))
write(hiv_time, file=file.path("hivPCMComparisonTimer_pcm.txt"))