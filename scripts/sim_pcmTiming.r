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

source(file.path(this.dir, "pcm_timing.r"))

args = commandArgs(trailingOnly = TRUE)

stopifnot(length(args) == 2)

filename = args[1]
storage.dir = args[2]

tree <- normalizePath(file.path(storage.dir, sprintf("%s_newick.txt", filename)))
print(tree)
data <- normalizePath(file.path(storage.dir, sprintf("%s_data.csv", filename)))

diff_mat <- normalizePath(file.path(storage.dir, sprintf("%s_diffVar.mat", filename)))
res_mat <- normalizePath(file.path(storage.dir, sprintf("%s_resVar.mat", filename)))

time <- PCMBaseTiming(tree, data, diff_mat, res_mat, fixed_tree = TRUE, Nreps=1000)

options(digits=10)
write(time, sprintf("%sTimer_pcm.txt", filename))

