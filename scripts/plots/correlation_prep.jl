using DataFrames, CSV, LinearAlgebra
using BeastUtils.Logs

const DATA_LABELS = ["corrs", "p", "signed_p", "corr_ref", "string_p", "fill_p"]
const P_FILL = "*"
const TRIM_P = true


function matrix_to_df(ind_labels::Vector{String}, data_labels::Vector{String},
        x::Matrix{}...)
    (m, n) = size(x[1])
    @assert m == n && m == length(ind_labels) &&
        length(data_labels) == length(x)
    for i = 1:length(x)
        a, b = size(x[i])
        @assert a == b && n == a
    end

    label_rows = repeat(ind_labels, outer=n)
    label_cols = repeat(ind_labels, inner=n)
    df = DataFrame(row = label_rows, col = label_cols)
    for i = 1:length(x)
        df[!, Symbol(data_labels[i])] = vec(x[i])
    end
    return df
end

function make_signed_pmat(corrs::Matrix{Float64}, p_vals::Matrix{Float64})
    n = size(corrs, 1)
    signed_p = zeros(n, n)
    for i = 1:n
        for j = 1:n
            if corrs[i, j] > 0.0
                signed_p[i, j] = p_vals[i, j]
            else
                signed_p[i, j] = -p_vals[i, j]
            end
        end
    end
    return signed_p
end

function make_corr_ref(p::Int)
    T = Union{Float64, Missing}
    corr_ref = Matrix{T}(undef, p, p)
    fill!(corr_ref, missing)
    for i = 1:p
        corr_ref[i, i] = 1.0
    end
    return corr_ref
end

function ind_convert(x::Int, p::Int)
    i = div(x - 1, p) + 1
    j = x - (i - 1) * p
    return [i, j]
end

function ind_convert(x::Vector{Int}, p::Int)
    return (x[1] - 1) * p + x[2]
end

function compute_sign_prob(corrmat, data)
    p = size(corrmat, 1)
    n, q = size(data)
    @assert q == p^2
    sign_mat = zeros(Bool, p, p)
    for i = 1:p
        for j = 1:p
            if corrmat[i, j] > 0.0
                sign_mat[i, j] = true
            end
        end
    end

    sign_vec = vec(sign_mat)
    counts = zeros(Int, q)
    for i = 1:n
        for j = 1:q
            test = data[i, j] > 0.0
            if test == sign_vec[j]
                counts[j] += 1
            end
        end
    end
    ps = counts ./ n
    pmat = reshape(ps, p, p)
    return pmat
end

function remove_data!(X::Matrix{}; upper::Bool=false,
        diagonal::Bool=false, lower::Bool=false)
    p = size(X, 1)
    if lower
        for i = 1:(p - 1)
            for j = (i + 1):p
                X[i, j] = missing
            end
        end
    end
    if upper
        for i = 2:p
            for j = 1:(i - 1)
                X[i, j] = missing
            end
        end
    end
    if diagonal
        for i = 1:p
            X[i, i] = missing
        end
    end
end

function remove_data!(X::Array{Union{Missing, T}}, val::T) where T <: Any
    n = length(X)
    for i = 1:n
        if !ismissing(X[i]) && X[i] == val
            X[i] = missing
        end
    end
end

function make_string_pmat(X::Matrix{T}, n::Int,
        sigfigs::Int) where T <: Union{Missing, Float64}
    p, q = size(X)
    S = Matrix{String}(undef, p, q)
    s = ""
    lower_bound = string((n - 1) / n)[1:(sigfigs + 2)]
    lb = parse(Float64, lower_bound)
    lower_bound = lower_bound[2:end]
    # lower_bound = join(["> ", lower_bound])
    lower_bound = P_FILL
    for i = 1:p
        for j = 1:q
            if round(X[i, j], digits = sigfigs) > lb
                S[i, j] = lower_bound
            else
                S[i, j] = string(X[i, j])[2:(sigfigs + 2)]
            end
        end
    end
    return S
end

function make_corrDF(logpath::String,
        outfile::String,
        labels_outfile::String,
        labels::Vector{String},
        header::String,
        new_order::Vector{Int};
        trim_p::Bool = TRIM_P,
        burnin::Float64 = 0.1)
    cols, data = get_log(logpath, header, burnin = burnin)
    n, q = size(data)
    p = Int(sqrt(q))
    corrmat = make_meanmatrix(data, 1, q)
    pmat = compute_sign_prob(corrmat, data)
    signed_pmat = make_signed_pmat(corrmat, pmat)
    corr_ref = make_corr_ref(p)
    string_pmat = make_string_pmat(pmat, n, 2)

    T = Union{Float64, Missing}
    S = Union{String, Missing}
    corrmat = convert(Matrix{T},corrmat[new_order, new_order])
    pmat = convert(Matrix{T}, pmat[new_order, new_order])
    signed_pmat = convert(Matrix{T}, signed_pmat[new_order, new_order])
    string_pmat = convert(Matrix{S}, string_pmat[new_order, new_order])

    remove_data!(corrmat, upper=true, diagonal=true)
    remove_data!(pmat, lower = true, diagonal=true)
    remove_data!(signed_pmat, lower = true, diagonal=true)
    remove_data!(string_pmat, lower=true, diagonal=true)

    #Starting asterisk stuff
    tsk_mat = Matrix{T}(undef, p, p)

    if trim_p
        for i = 1:q
            if !ismissing(string_pmat[i]) && string_pmat[i] == P_FILL
                tsk_mat[i] = 1.0
            end
        end
        remove_data!(string_pmat, P_FILL)
    end


    df = matrix_to_df(labels[new_order], DATA_LABELS,
        corrmat,
        pmat,
        signed_pmat,
        corr_ref,
        string_pmat,
        tsk_mat)
    CSV.write(outfile, df)
    CSV.write(labels_outfile, DataFrame(labels = labels[new_order]))
    return df
end

log_directory = joinpath(@__DIR__, "..", "..", "logs")
plot_directory = @__DIR__
