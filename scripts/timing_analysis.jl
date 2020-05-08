using Statistics, DataFrames, CSV
using BeastUtils.Logs, BeastUtils.MatrixUtils, BeastUtils.ESS


const HOUR_MULT = Dict("days" => 24.0, "hours" => 1.0, "minutes" => 1 / 60, "seconds" => 1 / 3600)

function get_time(path::String)
    lines = readlines(path)
    try
        time_line = lines[end - 1]
        time_line = split(time_line)
        val = time_line[1]
        val = parse(Float64, val)
        unit = time_line[2]
        hrs = val * HOUR_MULT[unit]
        return hrs
    catch
        last_line = lines[end]
        split_line = split(last_line)
        states = split_line[1]
        states = parse(Float64, states)
        time_per_mil = parse(Float64, split_line[end - 2])

        unit = split(split_line[end - 1], "/")[1]
        @assert split(split_line[end - 1], "/")[2] == "million"
        hrs = states / 1e6 * time_per_mil * HOUR_MULT[unit]
        return hrs
    end
end

function ess_per_hour(log_path::String, time_path::String,
        relevant_cols::Vector{Int})
    hrs = get_time(time_path)
    cols, data = get_log(log_path)
    esses = ess(data, relevant_cols)
    med = median(esses)
    mn = minimum(esses)
    return med / hrs, mn / hrs
end


function compare_ess_per_hour(log_names::Vector{String},
        time_names::Vector{String}, log_dir::String,
        time_dir::String, cols::Vector{Int})


    p = length(log_names)
    @assert length(time_names) == p
    mins = zeros(p)
    meds = zeros(p)
    for i = 1:p
        time_path = joinpath(time_dir, time_names[i])
        log_path = joinpath(log_dir, log_names[i])
        med, mn = ess_per_hour(log_path, time_path, cols)
        mins[i] = mn
        meds[i] = med
    end

    return mins, meds
end

function ess_stats(data_path::String, timing_path::String, col_start::String)
    cols = get_cols(data_path)
    f, l = find_cols(cols, col_start)
    return ess_stats(data_path, timing_path, collect(f:l))
end

function ess_stats(data_path::String, timing_path::String, rel_cols::Vector{Int})
    cols, data = get_log(data_path)
    return ess_stats(data, timing_path, rel_cols)
end

function ess_stats(data::Matrix{Float64}, timing_path::String,
        rel_cols::Vector{Int})

    esses = ess(data, rel_cols)
    n = size(data, 1)
    states = data[n, 1]
    time = get_time(timing_path)
    return esses, states, time
end

function get_states(data::AbstractArray{Float64, 2})
    return data[end, 1]
end

function average_ess_stats(data_dir::String, data_files::Array{String},
        timing_dir::String, timing_files::Array{String}, rel_cols::Vector{Int})

    data_path = joinpath(data_dir, data_files[1])
    timing_path = joinpath(timing_dir, timing_files[1])
    esses, states, time = ess_stats(data_path, timing_path, rel_cols)
    n = length(data_files)
    @assert length(timing_files) == n

    for i = 2:n
        data_path = joinpath(data_dir, data_files[i])
        timing_path = joinpath(timing_dir, timing_files[i])
        new_esses, new_states, new_time =
            ess_stats(data_path, timing_path, rel_cols)
        esses .+= new_esses
        states += new_states
        time += new_time
    end
    return esses / n, states / n, time / n
end

function get_stats(log_dir::String, timing_dir::String,
        filenames::Array{Array{String, N}, M},
        param_names::Vector{Vector{String}};
        corrections::Vector{Float64} = ones(length(filenames))) where {N, M}

    p = length(filenames)
    time_names = [["$x.txt" for x in y] for y in filenames]
    log_names = [["$x.log" for x in y] for y in filenames]

    min_ess = zeros(p)
    med_ess = zeros(p)
    states = zeros(p)
    hrs = zeros(p)

    for i = 1:p
        @show log_names[i]
        cols = get_cols(joinpath(log_dir, log_names[i][1]))
        # list_cols(cols)
        rel_cols = Vector{Int}(undef, 0)
        @show log_names[i][1]
        for j = 1:length(param_names[i])
            first, last = find_cols(cols, param_names[i][j])
            rel_cols = [rel_cols; matrix_inds(first, last)]
        end
        @show cols[rel_cols]
        esses, state, hr = average_ess_stats(log_dir, log_names[i],
            timing_dir, time_names[i], rel_cols)
        min_ess[i] = minimum(esses)
        med_ess[i] = median(esses)
        states[i] = state
        hrs[i] = hr
    end
    min_stats = hcat((min_ess ./ hrs),
        (states ./ hrs ./corrections),
        (corrections .* min_ess ./ states))
    med_stats = hcat((med_ess ./ hrs),
        (states ./ hrs ./corrections),
        (corrections .* med_ess ./ states))
    return min_stats, med_stats
end


function fold_change(stats::Array{Float64}, prs::Matrix{Int})
    @assert size(prs, 2) == 2
    @assert 2 * size(prs, 1) == size(stats, 1)
    n = size(prs, 1)
    p = size(stats, 2)
    folds = zeros(n, p)
    for i = 1:n
        for j = 1:p
            folds[i, j] = stats[prs[i, 1], j] / stats[prs[i, 2], j]
        end
    end
    return folds
end


log_dir = joinpath(@__DIR__, "..", "logs", "timing")
timing_dir = log_dir

# Corrections for sampling frequency
mammals_corrs = [1.0, 501/1]
hiv_corrs = [1.0, 11/1]

mammals_corrs_rm = [1.0, 52/2]
hiv_corrs_rm = [1.0, 52/2]


corrs = [mammals_corrs; hiv_corrs; mammals_corrs_rm; hiv_corrs_rm]


repeats = 10

filenames = ["mammals_diff", "mammals_newTiming", "hiv_diff", "hiv_newTiming",
    "mammals_res", "mammals_newTiming_rm", "hiv_res", "hiv_newTiming_rm"]

prs = [1 2; 3 4; 5 6; 7 8]

pnames = [["variance.precisionMatrix"], ["variance.precisionMatrix"],
        ["variance"], ["variance"],
        ["DiffVariance.precisionMatrix", "SampVariance.rm_precision"],
            ["variance.precisionMatrix", "residualVariance.residualPrecision"],
        ["diffVariance", "resVariance"], ["variance", "rm_variance"]]

r_filenames = [["$(x)_r$i" for i = 1:repeats] for x in filenames]


min_stats, med_stats = get_stats(log_dir, timing_dir, r_filenames, pnames,
    corrections = corrs)

min_folds = fold_change(min_stats, prs)
med_folds = fold_change(med_stats, prs)

df = DataFrame()
df.data = ["mammals", "mammals", "HIV", "HIV", "mammals", "mammals", "HIV", "HIV"]
df.model = ["diff", "diff", "diff", "diff", "diff and res", "diff and res", "diff and res", "diff and res"]
df.method = ["analytic", "sampling", "analytic", "sampling", "analytic", "sampling", "analytic", "sampling"]
df.minESSpHr = min_stats[:, 1]
df.medESSpHr = med_stats[:, 1]
df.minESSpSample = min_stats[:, 3]
df.medESSpSample = med_stats[:, 3]
df.SamplepHr = min_stats[:, 2]

CSV.write(joinpath(@__DIR__, "storage", "timing_summary.csv"), df)

### Code below used to write variables to latex-readable format
#
# meths = ["An", "Sam", "An", "Sam", "An", "Sam", "An", "Sam"]
# datasets = ["Mam", "Mam", "HIV", "HIV", "Mam", "Mam", "HIV", "HIV"]
# rms = ["", "", "", "", "RM", "RM", "RM", "RM"]
#
# var_names = ["EssPHr", "StPHr", "ESSPSt"]
# p = length(var_names)
# @assert p == size(min_stats, 2) == size(med_stats, 2)
# n = length(filenames)
# min_names = Matrix{String}(undef, n, p)
# med_names = Matrix{String}(undef, n, p)
#
# for i = 1:n
#     for j = 1:p
#         add_on = "$(var_names[j])$(meths[i])$(datasets[i])$(rms[i])"
#         min_names[i, j] = "min$add_on"
#         med_names[i, j] = "med$add_on"
#     end
# end
#
# m = size(prs, 1)
# @assert 2 * m == n
# min_fold_names = Matrix{String}(undef, m, p)
# med_fold_names = Matrix{String}(undef, m, p)
# fold_data = datasets[prs[:, 1]]
# fold_rms = rms[prs[:, 1]]
#
# for i = 1:m
#     for j = 1:p
#         add_on = "$(var_names[j])$(fold_data[i])$(fold_rms[i])Fold"
#         min_fold_names[i, j] = "min$add_on"
#         med_fold_names[i, j] = "med$add_on"
#     end
# end
#
# storage_path = joinpath(@__DIR__, "storage", "time_table.tex")
#
# insert_defs(storage_path, [min_names; med_names; min_fold_names; med_fold_names],
#     [min_stats; med_stats; min_folds; med_folds], overwrite = true, round_to = 2,
#     sigfigs = true, cut_zeros = true, add_commas = true)
