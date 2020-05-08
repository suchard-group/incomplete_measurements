using Revise

using CSV, DataFrames
using BeastUtils.Logs, BeastUtils.RTrees, BeastUtils.MatrixUtils
using ReadLogs, Trees2, MyFunctions #my packages

const VALUE_DIR = joinpath(@__DIR__, "storage", "simulation")
const LOG_DIR = joinpath(@__DIR__, "..", "logs", "simulation_study")

struct SimStats
    bias::Vector{Float64}
    mse::Vector{Float64}
    coverage::Float64
end

struct MatrixSimStats
    diagonal::SimStats
    off_diagonal::SimStats
end

struct SimVersion
    filename::String
    version::String
    N::Int #number of taxa
    P::Int #number of traits
    sparsity::Float64
    repeat::Int
end

struct SimAnalysis
    sv::SimVersion
    traits::SimStats
    diff_corr::MatrixSimStats
    res_corr::MatrixSimStats
    heritability::MatrixSimStats
end

function parse_filename(filename::String)
    base_name = split(filename, '.')[1]
    elements = split(base_name, '_')
    version = elements[1]

    @assert endswith(elements[2], "Taxa")
    N = parse(Int, elements[2][1:(end - length("Taxa"))])

    @assert endswith(elements[3], "Traits")
    P = parse(Int, elements[3][1:(end - length("Traits"))])

    @assert endswith(elements[4], "percMissing")
    sparsity = parse(Int, elements[4][1:(end - length("percMissing"))]) / 100

    repeat = 0
    if length(elements) > 4
        @assert length(elements) == 5
        repeat = parse(Int, elements[5])
    end

    return SimVersion(base_name, version, N, P, sparsity, repeat)
end

function parse_parameters(param_string::String)
    diff_line, res_line = split(param_string, '\n')

    diff_start = "diffusion_variance"
    res_start = "residual_variance"
    @assert startswith(diff_line, diff_start)
    @assert startswith(res_line, res_start)

    diff_vec = split(diff_line, ' ')[2:end]
    res_vec = split(res_line, ' ')[2:end]

    Σ_vec = parse.(Float64, diff_vec)
    Γ_vec = parse.(Float64, res_vec)

    q = length(Σ_vec)
    @assert length(Γ_vec) == q
    p = Int(round(sqrt(q)))
    @assert p * p == q

    Σ = reshape(Σ_vec, p, p)
    Γ = reshape(Γ_vec, p, p)
    return Σ, Γ
end

function df_to_data(df::DataFrame)
    nms = names(df)
    @assert nms[1] == :taxon

    taxa = Vector{String}(df[!, 1])

    p = length(nms) - 1
    n = length(taxa)

    data = fill(NaN, n, p)

    for j = 1:p
        for i = 1:n
            if !ismissing(df[i, j + 1])
                data[i, j] = df[i, j + 1]
            end
        end
    end

    return taxa, data
end

function make_heritability(tree::RTrees.PhyloTree, Σ::AbstractArray{Float64, 2},
                        Γ::AbstractArray{Float64, 2})


    diag_sum, all_sum = RTrees.tree_variance_sums(tree, standardize_tree = true)

    n = tree.n_tips
    cσ = 1.0 / n * diag_sum - 1.0 / (n^2) * all_sum
    cγ = (n - 1) / n

    return make_heritability(cσ, cγ, Σ, Γ)
end

function make_heritability(cσ::Float64, cγ::Float64,
                        Σ::AbstractArray{Float64, 2},
                        Γ::AbstractArray{Float64, 2})
    p = size(Σ, 1)
    @assert size(Γ) == size(Σ) == (p, p)

    H = zeros(p, p)
    for i = 1:p
        ti = cσ * Σ[i, i] + cγ * Γ[i, i]
        H[i, i] = cσ * Σ[i, i] / ti

        for j = (i + 1):p
            tj = cσ * Σ[j, j] + cγ * Γ[j, j]

            H[i, j] = cσ * Σ[i, j] / sqrt(ti * tj)
            H[j, i] = H[i, j]
        end
    end
    return H
end

function compute_stats(x::Vector{Float64}, value::Float64)
    sum_x = 0.0
    sum_square_x = 0.0

    n = length(x)
    for i = 1:n
        sum_x += x[i]
        sum_square_x += x[i]^2
    end

    mean_x = sum_x / n
    var_x = sum_square_x / n - mean_x^2

    bias = mean_x - value
    mse = var_x + bias^2

    hpd_interval = Logs.HPD_intervals(x, conf = 0.95) #95% hpd intervals

    covered = 0
    if hpd_interval[1] <= value <= hpd_interval[2]
        covered = 1
    end

    return bias, mse, covered
end

function compute_trait_stats(taxa::Vector{String}, data::Matrix{Float64},
                        log_cols::Vector{String}, log_data::Matrix{Float64},
                        log_range::UnitRange{Int}; col_start::String = "tip.traits.")

    n, p = size(data)
    n_missing = length(log_range)
    m = size(log_data, 1)

    buffer = zeros(m)

    bias = zeros(n_missing)
    mse = zeros(n_missing)
    n_covered = 0

    ind = 0
    for i in log_range
        ind += 1

        col_name = log_cols[i]
        col_name = col_name[(length(col_start) + 1):end]
        taxon, dim_string = split(col_name, '.')

        taxon_ind = findfirst(x -> x == taxon, taxa)
        dim = parse(Int, dim_string)

        true_value = data[taxon_ind, dim]

        buffer .= @view log_data[:, i]

        trait_bias, trait_mse, trait_covered = compute_stats(buffer, true_value)

        bias[ind] = trait_bias
        mse[ind] = trait_mse
        n_covered += trait_covered
    end

    return SimStats(bias, mse, n_covered / n_missing)
end

function compute_matrix_stats(X::AbstractArray{Float64, 2}, log_data::Matrix{Float64},
                            log_range::UnitRange{Int}) #assumes a symmetric matrix

    p = size(X, 1)
    @assert size(X, 2) == p
    @assert length(log_range) == p^2

    m = size(log_data, 1)
    buffer = zeros(m)

    diag_bias = zeros(p)
    diag_mse = zeros(p)
    diag_covered = 0

    p_off = div(p * (p - 1), 2)

    off_diag_bias = zeros(p_off)
    off_diag_mse = zeros(p_off)
    off_diag_covered = 0

    ind = 0
    off_ind = 0
    for i in log_range
        ind += 1

        row = div(ind - 1, p) + 1
        col = ind - (row - 1) * p

        if row > col
            continue #don't want to do double work on off-diagonals
        end

        buffer .= @view log_data[:, i]

        bias, mse, covered = compute_stats(buffer, X[row, col])

        if row == col
            diag_bias[row] = bias
            diag_mse[row] = mse
            diag_covered += covered

        else
            off_ind += 1
            off_diag_bias[off_ind] = bias
            off_diag_mse[off_ind] = mse
            off_diag_covered += covered
        end
    end

    diag_stats = SimStats(diag_bias, diag_mse, diag_covered / p)
    off_diag_stats = SimStats(off_diag_bias, off_diag_mse, off_diag_covered / p_off)

    return MatrixSimStats(diag_stats, off_diag_stats)
end

function sim_analysis(sv::SimVersion)
    filename = sv.filename
    println("starting $filename")

    current_dir = pwd()

    #Read in 'true' parameters used for simulation
    cd(VALUE_DIR)

    newick = read("$(filename)_newick.txt", String)
    params = read("$(filename)_parameters.txt", String)
    data_df = CSV.read("$(filename)_traitData.csv")

    tree = RTrees.parse_newick(newick)

    Σ, Γ = parse_parameters(params)
    Σcorr = MatrixUtils.cov2corr(Σ)
    Γcorr = MatrixUtils.cov2corr(Γ)

    H = make_heritability(tree, Σ, Γ) #TODO: check this

    taxa, data = df_to_data(data_df)

    n, p = size(data)
    @assert size(Σ) == size(Γ) == (p, p)

    #Read in log file and split it into components

    print("loading log file ... ")

    cd(LOG_DIR)
    cols, log_data = Logs.get_log("$(sv.filename).log")
    cd(current_dir)

    println("done")

    trait_range = 1:0
    if sv.sparsity > 0.0
        f, l = Logs.find_cols(cols, "tip.traits")
        trait_range = f:l
    end

    f, l = Logs.find_cols(cols, "correlation.inverted.diffusion.precision")
    Σcorr_range = f:l

    f, l = Logs.find_cols(cols, "correlation.inverted.residualPrecision")
    Γcorr_range = f:l

    f, l = Logs.find_cols(cols, "varianceProportionStatistic")
    H_range = f:l

    #Compute statistics for traits
    print("computing trait statistics ... ")
    trait_stats = compute_trait_stats(taxa, data, cols, log_data, trait_range)
    println("done")

    #Compute statistics for matrices
    print("computing matrix statistics ... ")
    Σcorr_stats = compute_matrix_stats(Σcorr, log_data, Σcorr_range)
    Γcorr_stats = compute_matrix_stats(Γcorr, log_data, Γcorr_range)
    H_stats = compute_matrix_stats(H, log_data, H_range)
    println("done")
    println("")

    return SimAnalysis(sv, trait_stats, Σcorr_stats, Γcorr_stats, H_stats)
end

const TRAIT = "trait"
const DIAGONAL = "diagonal"
const OFF_DIAGONAL = "offDiagonal"
const BIAS = "bias"
const MSE = "mse"
const DIFFUSION_CORRELATION = "diffCorr"
const RESIDUAL_CORRELATION = "resCorr"
const HERITABILITY = "heritability"
const NONE = "none"

function make_dataframe(sas::Array{SimAnalysis},
                        trait_path::String,
                        matrix_path::String,
                        coverage_path::String)
    n = length(sas)

    #compute size of data frame
    n_missing = [length(sas[i].traits.bias) for i = 1:n]
    dim_mats = [length(sas[i].diff_corr.diagonal.bias) for i = 1:n]

    # for each trait, there are 3 * p diagonal and 3 * p * (p - 1) off-diagonal entries
    n_matrix_slots = [3 * (p + div(p * (p - 1), 2)) for p in dim_mats]

    N_trait = sum(n_missing)
    N_matrix = sum(n_matrix_slots)
    n_coverage = 7
    N_coverage = n * n_coverage


    trait_df = DataFrame(run = Vector{String}(undef, N_trait),
                        rep = Vector{Int}(undef, N_trait),
                        nTaxa = Vector{Int}(undef, N_trait),
                        nTraits = Vector{Int}(undef, N_trait),
                        nObs = Vector{Float64}(undef, N_trait),
                        sparsity = Vector{Float64}(undef, N_trait),
                        bias = Vector{Float64}(undef, N_trait),
                        mse = Vector{Float64}(undef, N_trait),
                        isRandom = fill(true, N_trait))

    matrix_df = DataFrame(run = Vector{String}(undef, N_matrix),
                        rep = Vector{Int}(undef, N_matrix),
                        nTaxa = Vector{Int}(undef, N_matrix),
                        nTraits = Vector{Int}(undef, N_matrix),
                        nObs = Vector{Float64}(undef, N_matrix),
                        sparsity = Vector{Float64}(undef, N_matrix),
                        variable = Vector{String}(undef, N_matrix),
                        component = Vector{String}(undef, N_matrix),
                        bias = Vector{Float64}(undef, N_matrix),
                        mse = Vector{Float64}(undef, N_matrix),
                        isRandom = Vector{Bool}(undef, N_matrix))

    coverage_df = DataFrame(run = Vector{String}(undef, N_coverage),
                        rep = Vector{Int}(undef, N_coverage),
                        nTaxa = Vector{Int}(undef, N_coverage),
                        nTraits = Vector{Int}(undef, N_coverage),
                        nObs = Vector{Float64}(undef, N_coverage),
                        sparsity = Vector{Float64}(undef, N_coverage),
                        variable = Vector{String}(undef, N_coverage),
                        component = Vector{String}(undef, N_coverage),
                        coverage = Vector{Float64}(undef, N_coverage),
                        isRandom = Vector{Bool}(undef, N_coverage))

    # set up for matrix statistics
    variables = [DIFFUSION_CORRELATION, RESIDUAL_CORRELATION, HERITABILITY]
    var_dict = Dict(DIFFUSION_CORRELATION => :diff_corr,
                    RESIDUAL_CORRELATION => :res_corr,
                    HERITABILITY => :heritability)



    trait_ind = 1
    matrix_ind = 1

    for i = 1:n
        coverage_ind = (i - 1) * n_coverage + 1

        sa = sas[i]
        sv = sa.sv

        # n_obs = sv.N * sv.P - n_missing[i]
        n_obs = sv.N * (1 - sv.sparsity) * sv.P / (1 - sv.sparsity ^ sv.P) #expected number of observations

        trait_range = trait_ind:(trait_ind + n_missing[i] - 1)
        mat_range = matrix_ind:(matrix_ind + n_matrix_slots[i] - 1)
        cov_range = coverage_ind:(coverage_ind + n_coverage - 1)

        trait_df.run[trait_range] .= sv.version
        matrix_df.run[mat_range] .= sv.version
        coverage_df.run[cov_range] .= sv.version

        trait_df.rep[trait_range] .= sv.repeat
        matrix_df.rep[mat_range] .= sv.repeat
        coverage_df.rep[cov_range] .= sv.repeat

        trait_df.nTaxa[trait_range] .= sv.N
        matrix_df.nTaxa[mat_range] .= sv.N
        coverage_df.nTaxa[cov_range] .= sv.N

        trait_df.nTraits[trait_range] .= sv.P
        matrix_df.nTraits[mat_range] .= sv.P
        coverage_df.nTraits[cov_range] .= sv.P

        trait_df.sparsity[trait_range] .= sv.sparsity
        matrix_df.sparsity[mat_range] .= sv.sparsity
        coverage_df.sparsity[cov_range] .= sv.sparsity

        trait_df.nObs[trait_range] .= n_obs
        matrix_df.nObs[mat_range] .= n_obs
        coverage_df.nObs[cov_range] .= n_obs


        #missing trait values

        trait_df.bias[trait_range] .= sa.traits.bias
        trait_df.mse[trait_range] .= sa.traits.mse

        coverage_df.variable[coverage_ind] = TRAIT
        coverage_df.component[coverage_ind] = NONE
        coverage_df.coverage[coverage_ind] = sa.traits.coverage

        no_traits = length(sa.traits.bias) == 0
        coverage_df.isRandom[coverage_ind] = !no_traits




        # matrix stats setup
        p = dim_mats[i] #number diagonal entries
        q = div(p * (p - 1), 2) #number off-diagonal entries
        t = p + q #total number of entries


        n_var = length(variables)
        for j = 1:n_var
            s = (matrix_ind + (j - 1) * t) #start index for this variable
            matrix_df.variable[s:(s + t - 1)] .= variables[j]


            mat_stats = getfield(sa, var_dict[variables[j]])

            #diagonal
            diag_range = s:(s + p - 1)
            matrix_df.component[diag_range] .= DIAGONAL
            diag_stats = mat_stats.diagonal


            matrix_df.bias[diag_range] .= diag_stats.bias
            matrix_df.mse[diag_range] .= diag_stats.mse

            corr_diag = (variables[j] == "diffCorr" || variables[j] == "resCorr")
            matrix_df.isRandom[diag_range] .= !corr_diag

            diag_cov_ind = coverage_ind + 2 * (j - 1) + 1
            coverage_df.variable[diag_cov_ind] = variables[j]
            coverage_df.component[diag_cov_ind] = DIAGONAL
            coverage_df.coverage[diag_cov_ind] = diag_stats.coverage
            coverage_df.isRandom[diag_cov_ind] = !corr_diag

            #off-diagonal
            off_diag_range = (s + p):(s + t - 1)
            matrix_df.component[off_diag_range] .= OFF_DIAGONAL
            od_stats = mat_stats.off_diagonal

            matrix_df.bias[off_diag_range] .= od_stats.bias
            matrix_df.mse[off_diag_range] .= od_stats.mse
            matrix_df.isRandom[off_diag_range] .= true

            off_diag_cov_ind = coverage_ind + 2 * (j - 1) + 2
            coverage_df.variable[off_diag_cov_ind] = variables[j]
            coverage_df.component[off_diag_cov_ind] = OFF_DIAGONAL
            coverage_df.coverage[off_diag_cov_ind] = od_stats.coverage
            coverage_df.isRandom[off_diag_cov_ind] = true


        end

        trait_ind += n_missing[i]
        matrix_ind += n_matrix_slots[i]

    end

    CSV.write(trait_path, trait_df)
    CSV.write(matrix_path, matrix_df)
    CSV.write(coverage_path, coverage_df)
end



files = readdir(LOG_DIR)
deleteat!(files, findfirst(x -> x == ".gitignore", files))

n = length(files)
sas = Vector{SimAnalysis}(undef, n)



for i = 1:n
    sim_version = parse_filename(files[i])
    sas[i] = sim_analysis(sim_version)
    println("$i of $n completed\n")
end

storage_dir = joinpath(@__DIR__, "storage")
trait_path = joinpath(storage_dir, "trait_simulation.csv")
matrix_path = joinpath(storage_dir, "matrix_simulation.csv")
coverage_path = joinpath(storage_dir, "coverage_simulation.csv")

df = make_dataframe(sas, trait_path, matrix_path, coverage_path)
