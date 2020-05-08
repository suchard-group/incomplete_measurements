using BeastUtils.Logs
using Statistics, DataFrames, CSV


function compute_means(log_path::String; header::String = "validation")
    cols, data = get_log(log_path)
    validation_cols = Vector{Int}(undef, 0)
    for i = 1:length(cols)
        col = cols[i]
        if startswith(col, header)
            push!(validation_cols, i)
        end
    end

    validation_data = data[:, validation_cols]

    return mean(validation_data, dims = 1)
end

function process_means(mses::Matrix{Float64}, out_path::String,
        dims::Vector{String},
        proc::Vector{String};
        log_scale::Bool = false)

    n, p = size(mses)
    @assert length(dims) == p && length(proc) == p
    vdata = vec(mses)
    if log_scale
        vdata = log.(vdata)
    end
    vdims = repeat(dims, inner = n)
    vproc = repeat(proc, inner = n)
    df = DataFrame()
    df[!, :MSE] = vdata
    df[!, :Dimension] = vdims
    df[!, :Process] = vproc
    CSV.write(out_path, df)
end


function analyze_results(template_names::Vector{String}, n_reps::Int,
                        log_dir::String, csv_path::String,
                        dimensionality::Vector{String}, processes::Vector{String})

    n = length(template_names)
    mses = Matrix{Float64}(undef, 0, n)
    ind = 1
    for i = 1:n_reps
        m = 0
        @show i

        for j = 1:n
            log_name = "$(template_names[j])_$i.log"
            log_path = joinpath(log_dir, log_name)
            means = compute_means(log_path)
            m = length(means)
            if j == 1
                mses = [mses; zeros(m, n)]
            end

            mses[ind:(ind + m - 1), j] .= vec(means)
        end
        ind += m
    end

    process_means(mses, csv_path, dimensionality, processes, log_scale = true)
    return mses
end

log_dir = joinpath(@__DIR__, "..", "logs", "hiv_prediction")


original_names = ["hiv_mse_d", "hiv_mse_dr", "hiv_mse_du", "hiv_mse_dru"]
dims = ["Bivariate", "Bivariate", "Univariate", "Univariate"]

procs = ["Diffusion Only", "Diffusion and\nResidual",
            "Diffusion Only", "Diffusion and\nResidual"]

plot_csv_path = joinpath(@__DIR__, "plots", "hiv_prediction.csv")


p = 50

analyze_results(original_names, p, log_dir, plot_csv_path, dims, procs)
