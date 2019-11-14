using Revise

cd(@__DIR__)

using LinearAlgebra, CSV, LightXML, Distributions
using ReadLogs, MyFunctions, XMLConstructor2, Trees2, DiffusionSimulation2, DataStorage #Personal modules

### Set up diffusion and residual variance based on real data

function setup_variances(log_path::String,
        Σ_outpath::String,
        Γ_outpath::String,
        Σ_start::String,
        Γ_start::String)

    cols, data = get_log(log_path)

    Σ_hat = make_meanmatrix(data, cols, Σ_start)
    Γ_hat = make_meanmatrix(data, cols, Γ_start)

    make_symmetric!(Σ_hat)
    make_symmetric!(Γ_hat)

    @assert isposdef(Σ_hat)
    @assert isposdef(Γ_hat)

    p = size(Σ_hat, 1)
    @assert size(Γ_hat) == (p, p)

    Σ_string = join(Σ_hat, ' ')
    Γ_string = join(Γ_hat, ' ')

    write(Σ_outpath, Σ_string)
    write(Γ_outpath, Γ_string)

    return Σ_hat, Γ_hat
end

function setup_xml(data_path::String, newick_path::String, xml_dir::String,
                    filename::String, Σ::Matrix{Float64}, Γ::Matrix{Float64})

    df = CSV.read(data_path)
    taxa, data = XMLConstructor2.df_to_matrix(df)
    newick = read(newick_path, String)

    return setup_xml(taxa, data, newick, Σ, Γ, xml_dir, filename)
end

function setup_xml(taxa::Vector{String}, data::Matrix{Float64}, newick::String,
                    Σ::Matrix{Float64}, Γ::Matrix{Float64},
                    xml_dir::String, filename::String)

    bx = XMLConstructor2.make_timing_MBD_XML(data, taxa, newick, chain_length = 1_000)

    bx.mcmc_el.filename = filename
    bx.mcmc_el.screen_logEvery = 100
    bx.timer_el.filename = "$(filename)Timer_beast.txt"
    bx.MBD_el.precision = inv(Σ)
    bx.extension_el.precision = inv(Γ)

    save_file(XMLConstructor2.make_xml(bx), joinpath(xml_dir, "$filename.xml"))

    return bx
end

# Universal variables

log_dir = joinpath("..", "logs")
data_dir = joinpath("..", "data")
xml_dir = joinpath("..", "xml", "PCMBase_comparison")


### Real data examples

Σ_start = "inverse.diffusion.precision"
Γ_start = "inverse.residualPrecision"

# Mammals example

mammals_filename = "mammalsPCMComparison"

mammals_log_path = joinpath(log_dir, "mammals.log")
mammals_data_path = joinpath(data_dir, "mammals_log_data.csv")
mammals_newick_path = joinpath(data_dir, "mammals_trimmed_newick.txt")

mammals_Σ_outpath = joinpath(@__DIR__, "storage", "mammals_diff_mat.mat")
mammals_Γ_outpath = joinpath(@__DIR__, "storage", "mammals_res_mat.mat")


mammals_Σ, mammals_Γ = setup_variances(mammals_log_path,
                                        mammals_Σ_outpath,
                                        mammals_Γ_outpath,
                                        Σ_start, Γ_start)

setup_xml(mammals_data_path, mammals_newick_path, xml_dir,
            mammals_filename, mammals_Σ, mammals_Γ)

# HIV example
hiv_filename = "hivPCMComparison"
hiv_log_path = joinpath(log_dir, "hiv.log")
hiv_data_path = joinpath(data_dir, "hiv_processed_data.csv")
hiv_newick_path = joinpath(data_dir, "hiv_newick.txt")

hiv_Σ_outpath = joinpath(@__DIR__, "storage", "hiv_diff_mat.mat")
hiv_Γ_outpath = joinpath(@__DIR__, "storage", "hiv_res_mat.mat")

hiv_Σ, hiv_Γ = setup_variances(hiv_log_path,
                                        hiv_Σ_outpath,
                                        hiv_Γ_outpath,
                                        Σ_start, Γ_start)

setup_xml(hiv_data_path, hiv_newick_path, xml_dir,
            hiv_filename, hiv_Σ, hiv_Γ)

# Prokaryote example
prok_filename = "prokPCMComparison"
prok_log_path = joinpath(log_dir, "prokaryotes.log")
prok_data_path = joinpath(data_dir, "prokaryotes_processed_data.csv")
prok_newick_path = joinpath(data_dir, "prokaryotes_newick.txt")

prok_Σ_outpath = joinpath(@__DIR__, "storage", "prok_diff_mat.mat")
prok_Γ_outpath = joinpath(@__DIR__, "storage", "prok_res_mat.mat")

prok_Σ, prok_Γ = setup_variances(prok_log_path,
                                        prok_Σ_outpath,
                                        prok_Γ_outpath,
                                        "diffVar", "resVar")

setup_xml(prok_data_path, prok_newick_path, xml_dir,
            prok_filename, prok_Σ, prok_Γ)


### Simulated examples
import Random
Random.seed!(666)


storage_dir = joinpath(@__DIR__, "storage", "PCMBase_comparison")

sim_name = "sim"
ns = [100, 1000, 10000]
ps = [2, 10, 20]
reps = 1
sparsity = 0.5

nn = length(ns)
np = length(ps)

filenames = Vector{String}(undef, nn * np * reps)
file_ind = 1

for i = 1:nn

    n = ns[i]

    #Setup taxa names
    taxa = ["taxon_$l" for l = 1:n]

    for j = 1:np

        p = ps[j]

        println("n = $n, p = $p")


        #Setup prior distribution for Σ and Γ
        W = Wishart(p, Matrix(Diagonal(ones(p))))

        for k = 1:reps

            #Draw diffusion and residual variances from prior distribution
            Σ = inv(rand(W))
            Γ = inv(rand(W))

            MyFunctions.make_symmetric!(Σ)
            MyFunctions.make_symmetric!(Γ)

            #Random tree
            tree = Trees2.rtree(taxa)
            newick = Trees2.make_newick(tree)

            #Simulate data
            tdm = DiffusionSimulation2.TreeDiffusionModel(tree, Σ, Γ, zeros(p))
            data = DiffusionSimulation2.simulate_data(tdm, taxa)
            DataStorage.make_sparse!(data, sparsity)

            #setup xml
            filename = "$(sim_name)_$(n)Taxa_$(p)Traits_$k"
            setup_xml(taxa, data, newick, Σ, Γ, xml_dir, filename)

            #save data to storage directory so R can access it

            cd(storage_dir)
            DataStorage.store_parameters(
                    ["$(filename)_diffVar.mat", "$(filename)_resVar.mat"],
                    [Σ, Γ])

            DataStorage.store_data("$(filename)_data.csv", taxa, data)

            write("$(filename)_newick.txt", newick)

            cd(@__DIR__)

            local local_ind = file_ind
            filenames[local_ind] = filename
            global file_ind += 1

        end
    end
end

write("sim_timing_files.txt", join(filenames, '\n'))
