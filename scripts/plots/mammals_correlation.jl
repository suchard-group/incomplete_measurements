include(joinpath(@__DIR__, "correlation_prep.jl"))

mammals_logpath = joinpath(log_directory, "mammals.log")
mammals_outfile = joinpath(plot_directory, "mammalsCorrelation.csv")
mammals_labels_outfile = joinpath(plot_directory, "mammalsCorrelationLabels.csv")

mammals_header = "correlation.inverted.diffusion.precision.diffusion.precision"

mammals_labels = ["body mass", "age at first birth", "gestation time",
                    "litter size", "litters per year", "neonatal body mass",
                    "weaning age", "reproductive lifespan"]

mammals_neworder = [1, 3, 7, 6, 2, 8, 4, 5]

print("Processing mammals log file ... ")
mammals_df = make_corrDF(mammals_logpath, mammals_outfile,
    mammals_labels_outfile, mammals_labels, mammals_header, mammals_neworder)
println("done")
