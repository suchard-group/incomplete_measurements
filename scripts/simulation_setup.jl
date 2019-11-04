using Revise

using LinearAlgebra, LightXML, DataFrames, CSV
using XMLConstructor2, DiffusionSimulation2, Trees2, ReadLogs, MyFunctions #personal packages

function sim_xml(xml_dir::String,
                newick::String,
                Σ::AbstractArray{Float64, 2}, #diffusion variance
                Γ::AbstractArray{Float64, 2}, #residual variance
                sparsity::Float64,
                base_name::String;
                standardize_tree::Bool = true,
                rep::Int = 0)

    tree = Trees2.parse_newick(newick)
    if standardize_tree
        max_height = maximum([Trees2.distance_to_root(tree, i) for i = 1:tree.n_tips])
        tree.edge_lengths ./= max_height
    end
    taxa = tree.tip_labels
    p = size(Σ, 1)

    tdm = DiffusionSimulation2.TreeDiffusionModel(tree, Σ, Γ, zeros(p))
    data = DiffusionSimulation2.simulate_data(tdm, taxa)
    n, p = size(data)

    mis_data = copy(data)
    make_sparse!(mis_data, sparsity)

    bx = XMLConstructor2.make_validation_MBD_XML(mis_data, data, taxa, newick,
                                                chain_length = 20_000)
    bx.mcmc_el.file_logEvery = 50
    bx.mcmc_el.screen_logEvery = 100

    perc_missing = Int(round(sparsity * 100))
    filename = "$(base_name)_$(n)Taxa_$(p)Traits_$(perc_missing)percMissing"

    if rep > 0
        filename = "$(filename)_$(rep)"
    end

    xml_name = "$filename.xml"
    for file in readdir(xml_dir)
        if file == xml_name
            error("File \"$xml_name\" already exists. Please delete it to confirm you want to overwrite it and its associate data")
        end
    end

    bx.mcmc_el.filename = filename
    xdoc = XMLConstructor2.make_xml(bx)
    save_file(xdoc, joinpath(xml_dir, xml_name))
    bx = nothing
    free(xdoc)

    current_dir = pwd()
    cd(joinpath(@__DIR__, "storage", "simulation"))

    store_data("$(filename)_traitData.csv", taxa, data)
    write("$(filename)_newick.txt", newick)

    param_string = "diffusion_variance $(join(Σ, ' '))"
    param_string = "$param_string\nresidual_variance $(join(Γ, ' '))"
    write("$(filename)_parameters.txt", param_string)

    cd(current_dir)

end


function sim_from_example(xml_dir::String,
                            log_path::String, newick_path::String,
                            diff_string::String, res_string::String,
                            sparsity::Float64, base_name::String)

    cols, data = ReadLogs.get_log(log_path)
    Σ_hat = ReadLogs.make_meanmatrix(data, cols, diff_string)
    Γ_hat = ReadLogs.make_meanmatrix(data, cols, res_string)

    make_symmetric!(Σ_hat)
    make_symmetric!(Γ_hat)


    @assert isposdef(Σ_hat)
    @assert isposdef(Γ_hat)
    @assert size(Σ_hat) == size(Γ_hat)

    p = size(Σ_hat, 1)

    newick = read(newick_path, String)

    sim_xml(xml_dir, newick, Σ_hat, Γ_hat, sparsity, base_name, standardize_tree = true)
end


function make_sparse!(X::Matrix{Float64}, sparsity::Float64)

    n, p = size(X)
    p_buffer = zeros(p)

    for i = 1:n

        all_missing = true
        p_buffer .= X[i, :]

        for j = 1:p
            if rand() < sparsity
                X[i, j] = NaN
            else
                all_missing = false
            end
        end

        if all_missing #don't want a taxon with no observations
            ind = rand(1:p)
            X[i, ind] = p_buffer[ind]
        end
    end
end

function store_data(path::String, taxa::Vector{String}, data::Matrix{Float64})
    df = data_to_df(taxa, data)
    CSV.write(path, df)
end

function data_to_df(taxa::Vector{String}, data::Matrix{Float64})
    col_names = ["trait_$i" for i = 1:size(data, 2)]
    return data_to_df(taxa, data, col_names)
end

function data_to_df(taxa::Vector{String}, data::Matrix{Float64}, col_names::Vector{String})
    n, p = size(data)
    @assert length(taxa) == n
    @assert length(col_names) == p
    df = DataFrame()
    df.taxon = taxa

    mis_inds = MyFunctions.findnans(data)
    mis_data = convert(Matrix{Union{Float64, Missing}}, data)
    mis_data[mis_inds] .= missing

    for i = 1:p
        df[!, Symbol(col_names[i])] = mis_data[:, i]
    end

    return df
end





Ns = [100, 1000, 10000]
Ps = [3, 10, 20]
miss_perc = [0.0, 0.25, 0.5, 0.75]
heritability = ["high", "low"]

xml_dir = joinpath(@__DIR__, "..", "xml", "simulation_study")

#simulations on real data sets
log_dir = joinpath(@__DIR__, "..", "logs")
data_dir = joinpath(@__DIR__, "..", "data")
diff_start = "inverse.diffusion.precision.diffusion.precision"
res_start = "inverse.residualPrecision.residualPrecision"

mammals_log_path = joinpath(log_dir, "mammals.log")
mammals_newick_path = joinpath(data_dir, "mammals_trimmed_newick.txt")

# sim_from_example(xml_dir, mammals_log_path, mammals_newick_path, diff_start, res_start, 0.5, "mammalsSim")

hiv_log_path = joinpath(log_dir, "hiv.log")
hiv_newick_path = joinpath(data_dir, "hiv_newick.txt")
# sim_from_example(xml_dir, hiv_log_path, hiv_newick_path, diff_start, res_start, 0.1, "hivSim")




println("done")
#test
# import Likelihoods, Trees, Random
# seed = 666
# Random.seed!(seed)
#
#
# n = 5
# p = 2
# taxa = ["taxa$i" for i = 1:n]
#
# Σ = randn(p, p)
# Σ = Σ * Σ'
# Γ = randn(p, p)
# Γ = Γ * Γ'
#
# tree = Trees2.rtree(taxa, seed)
# tree.edge_lengths[rand(1:(2 * n - 2))] = 0.0
# newick = Trees2.make_newick(tree)
#
# tdm2 = DiffusionSimulation2.TreeDiffusionModel(tree, Σ, Γ, zeros(p))
# tdm = Likelihoods.TreeDataModel(zeros(n, p), taxa, newick, standardize_tree = false)
# tdm.Σ = Hermitian(Σ)
# tdm.Γ = Hermitian(Γ)
# tdm.pss = Inf
#  V = Likelihoods.data_distribution(tdm)[2].Σ.mat
#
# cov = zeros(n * p, n * p)
# reps = 10_000_000
# for i = 1:reps
#     d = vec(DiffusionSimulation2.simulate_data(tdm2, taxa))
#     LinearAlgebra.BLAS.syr!('U', 1.0, d, cov)
# end
#
# make_symmetric!(cov)
# cov ./= reps
# display(V - cov)
# display(maximum(abs.(V - cov)))