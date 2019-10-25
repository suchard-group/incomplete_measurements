using Revise

using CSV, DataFrames
using ReadLogs, Trees2, MyFunctions #my packages

const VALUE_DIR = joinpath(@__DIR__, "storage")
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

function make_heritability(tree::Trees2.PhyloTree, Σ::AbstractArray{Float64, 2},
                        Γ::AbstractArray{Float64, 2})


    diag_sum, all_sum = Trees2.tree_variance_sums(tree, standardize_tree = true)

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

    hpd_interval = ReadLogs.HPD_intervals(x, conf = 0.95) #95% hpd intervals

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

    current_dir = pwd()

    #Read in 'true' parameters used for simulation
    cd(VALUE_DIR)

    newick = read("$(filename)_newick.txt", String)
    params = read("$(filename)_parameters.txt", String)
    data_df = CSV.read("$(filename)_traitData.csv")

    tree = Trees2.parse_newick(newick)

    Σ, Γ = parse_parameters(params)
    Σcorr = MyFunctions.cov2corr(Σ)
    Γcorr = MyFunctions.cov2corr(Γ)

    H = make_heritability(tree, Σ, Γ) #TODO: check this

    display(Σcorr)
    display(Γcorr)
    display(H)

    taxa, data = df_to_data(data_df)

    n, p = size(data)
    @assert size(Σ) == size(Γ) == (p, p)

    #Read in log file and split it into components

    print("loading log file ... ")

    cd(LOG_DIR)
    cols, log_data = ReadLogs.get_log("$(sv.filename).log")
    cd(current_dir)

    println("done")

    f, l = ReadLogs.find_cols(cols, "tip.traits")
    trait_range = f:l

    f, l = ReadLogs.find_cols(cols, "correlation.inverted.diffusion.precision")
    Σcorr_range = f:l

    f, l = ReadLogs.find_cols(cols, "correlation.inverted.residualPrecision")
    Γcorr_range = f:l

    f, l = ReadLogs.find_cols(cols, "varianceProportionStatistic")
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

    return SimAnalysis(sv, trait_stats, Σcorr_stats, Γcorr_stats, H_stats)

end

files = readdir(LOG_DIR)
n = length(files)
sas = Vector{SimAnalysis}(undef, n)

x = 0
for i = 1:n #TODO: change back to 1
    sim_version = parse_filename(files[i])
    sas[i] = sim_analysis(sim_version)
end


#testing
using Plots

function display_analysis(sa::SimAnalysis)
    @show sa.traits.coverage
    @show sa.diff_corr.off_diagonal.coverage
    @show sa.res_corr.off_diagonal.coverage
    @show sa.heritability.diagonal.coverage
    @show sa.heritability.off_diagonal.coverage
end

sa = sas[1]
display_analysis(sa)

using Distributions

p = 8
W = Wishart(p, Matrix(Diagonal(ones(p))))


newick = read(joinpath(VALUE_DIR, "$(sa.sv.filename)_newick.txt"), String)
tree = Trees2.parse_newick(newick)

diag_sum, all_sum = Trees2.tree_variance_sums(tree, standardize_tree = true)

n = tree.n_tips
cσ = 1.0 / n * diag_sum - 1.0 / (n^2) * all_sum
cγ = (n - 1) / n



# n = 100000
# d = zeros(n, p)
#
# for i = 1:n
#     Σ = (rand(W))
#     Γ = (rand(W))
#
#     H = make_heritability(cσ, cγ, Σ, Γ)
#
#     d[i, :] .= diag(H)
# end

H = make_heritability(tree, Diagonal(ones(p)), Diagonal(ones(p)))
