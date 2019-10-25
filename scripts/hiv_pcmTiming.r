this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

source("pcm_timing.r")

mammals_tree <- file.path(data.dir, "hiv_newick.txt")
mammals_data <- file.path(data.dir, "hiv_processed_data.csv")
mammals_diff_mat <- file.path("storage", "hiv_diff_mat.mat")
mammals_res_mat <- file.path("storage", "hiv_res_mat.mat")

mammals_time <- PCMBaseTiming(mammals_tree, mammals_data, mammals_diff_mat, mammals_res_mat, fixed_tree = TRUE, Nreps=1000)

setwd(this.dir)
write(mammals_time, file=file.path("..", "logs", "PCMBase_timing", "hivPCMComparisonTimer_pcm.txt"))