using LinearAlgebra, LightXML, DataFrames, CSV
using BeastUtils.XMLConstructor, BeastUtils.DiffusionSimulation,
        BeastUtils.RTrees, BeastUtils.Logs, BeastUtils.MatrixUtils,
        BeastUtils.DataStorage

function sim_xml(xml_dir::String,
                newick::String,
                Σ::AbstractArray{Float64, 2}, #diffusion variance
                Γ::AbstractArray{Float64, 2}, #residual variance
                sparsity::Float64,
                base_name::String;
                standardize_tree::Bool = true,
                rep::Int = 0)

    tree = RTrees.parse_newick(newick)
    if standardize_tree
        max_height = maximum([RTrees.distance_to_root(tree, i) for i = 1:tree.n_tips])
        tree.edge_lengths ./= max_height
    end
    taxa = tree.tip_labels
    p = size(Σ, 1)

    tdm = DiffusionSimulation.TreeDiffusionModel(tree, Σ, Γ, zeros(p))
    data = DiffusionSimulation.simulate_data(tdm, taxa)
    n, p = size(data)

    mis_data = copy(data)
    DataStorage.make_sparse!(mis_data, sparsity)

    bx = XMLConstructor.make_validation_MBD_XML(mis_data, data, taxa, newick,
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
            println("WARNING: File \"$xml_name\" already exists. File was NOT overwritten. Delete the file to create new xml.")
            return filename
        end
    end

    bx.mcmc_el.filename = filename
    xdoc = XMLConstructor.make_xml(bx)
    save_file(xdoc, joinpath(xml_dir, xml_name))
    bx = nothing
    free(xdoc)

    current_dir = pwd()
    cd(joinpath(@__DIR__, "storage", "simulation"))

    DataStorage.store_data("$(filename)_traitData.csv", taxa, data)
    write("$(filename)_newick.txt", newick)

    param_string = "diffusion_variance $(join(Σ, ' '))"
    param_string = "$param_string\nresidual_variance $(join(Γ, ' '))"
    write("$(filename)_parameters.txt", param_string)

    cd(current_dir)

    println("$(xml_name) completed")

    return filename
end


function sim_from_example(xml_dir::String,
                            log_path::String, newick_path::String,
                            diff_string::String, res_string::String,
                            sparsity::Vector{Float64},
                            ns::Vector{Int},
                            base_name::String;
                            repeats::Int = 0)

    cols, data = Logs.get_log(log_path)
    Σ_hat = Logs.make_meanmatrix(data, cols, diff_string)
    Γ_hat = Logs.make_meanmatrix(data, cols, res_string)

    MatrixUtils.make_symmetric!(Σ_hat)
    MatrixUtils.make_symmetric!(Γ_hat)


    @assert isposdef(Σ_hat)
    @assert isposdef(Γ_hat)
    @assert size(Σ_hat) == size(Γ_hat)

    p = size(Σ_hat, 1)

    newick = read(newick_path, String)
    tree = RTrees.parse_newick(newick)

    for s in sparsity
        for n in ns

            trimmed_tree = RTrees.trim_to_n(tree, n)
            trimmed_newick = RTrees.make_newick(trimmed_tree)


            if repeats == 0

                sim_xml(xml_dir, trimmed_newick, Σ_hat, Γ_hat, s,
                        base_name,
                        standardize_tree = true)

            else
                @assert repeats > 0
                for i = 1:repeats
                    sim_xml(xml_dir, trimmed_newick, Σ_hat, Γ_hat, s,
                            base_name,
                            standardize_tree = true, rep = i)
                end
            end
        end
    end
end



reps = 10

xml_dir = joinpath(@__DIR__, "..", "xml", "simulation_study")

### Simulations on real data sets



log_dir = joinpath(@__DIR__, "..", "logs")
data_dir = joinpath(@__DIR__, "..", "data")
diff_start = "inverse.diffusion.precision.diffusion.precision"
res_start = "inverse.residualPrecision.residualPrecision"


#mammals
mammals_log_path = joinpath(log_dir, "mammals.log")
mammals_newick_path = joinpath(data_dir, "mammals_trimmed_newick.txt")
mammals_mis_percs = [0.0, 0.25, 0.5, 0.75]
mammals_ns = [100, 500, 1000, 3649]

sim_from_example(xml_dir, mammals_log_path, mammals_newick_path,
                diff_start, res_start, mammals_mis_percs, mammals_ns,
                "mammalsSim", repeats = reps)


#hiv
hiv_log_path = joinpath(log_dir, "hiv.log")
hiv_newick_path = joinpath(data_dir, "hiv_newick.txt")
hiv_mis_percs = [0.0, 0.25, 0.5]
hiv_ns = [100, 500, 1000, 1536]

sim_from_example(xml_dir, hiv_log_path, hiv_newick_path,
                diff_start, res_start, hiv_mis_percs, hiv_ns,
                "hivSim", repeats = reps)


#prokaryotes
prok_log_path = joinpath(log_dir, "prokaryotes.log")
prok_newick_path = joinpath(data_dir, "prokaryotes_newick.txt")
prok_diff_start = "diffVar"
prok_res_start = "resVar"
prok_mis_percs = mammals_mis_percs
prok_ns = [100, 500, 705]

sim_from_example(xml_dir, prok_log_path, prok_newick_path,
                prok_diff_start, prok_res_start, prok_mis_percs, prok_ns,
                "prokSim", repeats = reps)


println("done")
