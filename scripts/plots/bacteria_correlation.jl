include(joinpath(@__DIR__, "correlation_prep.jl"))

bacteria_logpath = joinpath(log_directory, "prokaryotes.log")
bacteria_outfile = joinpath(plot_directory, "bacteriaCorrelation.csv")
bacteria_labels_outfile = joinpath(plot_directory, "bacteriaCorrelationLabels.csv")

bacteria_header = "diffCorr.traits.precision.col"

bacteria_labels = ["cell diameter", "optimal pH", "optimal temperature",
    "cell length", "genome length", "GC percent", "CDS length"]
bacteria_neworder = [5, 7, 3, 6, 1, 4, 2]



print("Processing prokaryotes log file ... ")
bacteria_df = make_corrDF(bacteria_logpath, bacteria_outfile,
    bacteria_labels_outfile, bacteria_labels, bacteria_header,
    bacteria_neworder)
println("done")
