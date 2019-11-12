using Revise

cd(@__DIR__)

using LinearAlgebra, CSV, LightXML
using ReadLogs, MyFunctions, XMLConstructor2 #Personal modules

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

function setup_xml(data_path::String, newick_path::String, xml_path::String,
                    filename::String, Σ::Matrix{Float64}, Γ::Matrix{Float64})

    df = CSV.read(data_path)
    taxa, data = XMLConstructor2.df_to_matrix(df)
    newick = read(newick_path, String)

    bx = XMLConstructor2.make_timing_MBD_XML(data, taxa, newick, chain_length = 1_000)

    bx.mcmc_el.filename = filename
    bx.mcmc_el.screen_logEvery = 100
    bx.timer_el.filename = "$(filename)Timer_beast.txt"
    bx.MBD_el.precision = inv(Σ)
    bx.extension_el.precision = inv(Γ)
    save_file(XMLConstructor2.make_xml(bx), xml_path)
    return bx
end

# Universal variables

log_dir = joinpath("..", "logs")
data_dir = joinpath("..", "data")
xml_dir = joinpath("..", "xml", "PCMBase_comparison")

Σ_start = "inverse.diffusion.precision"
Γ_start = "inverse.residualPrecision"

# Mammals example

mammals_filename = "mammalsPCMComparison"

mammals_log_path = joinpath(log_dir, "mammals.log")
mammals_data_path = joinpath(data_dir, "mammals_log_data.csv")
mammals_newick_path = joinpath(data_dir, "mammals_trimmed_newick.txt")
mammals_xml_path = joinpath(xml_dir, "$(mammals_filename).xml")

mammals_Σ_outpath = joinpath(@__DIR__, "storage", "mammals_diff_mat.mat")
mammals_Γ_outpath = joinpath(@__DIR__, "storage", "mammals_res_mat.mat")


mammals_Σ, mammals_Γ = setup_variances(mammals_log_path,
                                        mammals_Σ_outpath,
                                        mammals_Γ_outpath,
                                        Σ_start, Γ_start)

setup_xml(mammals_data_path, mammals_newick_path, mammals_xml_path,
            mammals_filename, mammals_Σ, mammals_Γ)

hiv_filename = "hivPCMComparison"
hiv_log_path = joinpath(log_dir, "hiv.log")
hiv_data_path = joinpath(data_dir, "hiv_processed_data.csv")
hiv_newick_path = joinpath(data_dir, "hiv_newick.txt")
hiv_xml_path = joinpath(xml_dir, "$(hiv_filename).xml")

hiv_Σ_outpath = joinpath(@__DIR__, "storage", "hiv_diff_mat.mat")
hiv_Γ_outpath = joinpath(@__DIR__, "storage", "hiv_res_mat.mat")

hiv_Σ, hiv_Γ = setup_variances(hiv_log_path,
                                        hiv_Σ_outpath,
                                        hiv_Γ_outpath,
                                        Σ_start, Γ_start)

setup_xml(hiv_data_path, hiv_newick_path, hiv_xml_path,
            hiv_filename, hiv_Σ, hiv_Γ)

prok_filename = "prokPCMComparison"
prok_log_path = joinpath(log_dir, "prokaryotes.log")
prok_data_path = joinpath(data_dir, "prokaryotes_processed_data.csv")
prok_newick_path = joinpath(data_dir, "prokaryotes_newick.txt")
prok_xml_path = joinpath(xml_dir, "$(prok_filename).xml")

prok_Σ_outpath = joinpath(@__DIR__, "storage", "prok_diff_mat.mat")
prok_Γ_outpath = joinpath(@__DIR__, "storage", "prok_res_mat.mat")

prok_Σ, prok_Γ = setup_variances(prok_log_path,
                                        prok_Σ_outpath,
                                        prok_Γ_outpath,
                                        "diffVar", "resVar")

setup_xml(prok_data_path, prok_newick_path, prok_xml_path,
            prok_filename, prok_Σ, prok_Γ)
