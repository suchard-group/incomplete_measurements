using Statistics, DataFrames, CSV

function pcm_time(path::String)
    s = read(path, String)
    string_times = split(s, ' ')
    return parse(Float64, string_times[3])
end

function beast_time(path::String)
    s = read(path, String)
    s = rstrip(s)
    string_time, unit = split(s, ' ')
    tm = parse(Float64, string_time)
    unit_mult = Dict("seconds" => 1, "minutes" => 6, "hours" => 3600)
    return tm * unit_mult[unit]
end

function get_time(path::String, method::String)
    if method == "beast"
        return beast_time(path)
    elseif method == "pcm"
        return pcm_time(path)
    else
        error("The second argument must be \"beast\" or \"pcm\".")
    end
end

const Nreps = 1_000
const Nruns = 6

time_dir = joinpath(@__DIR__, "..", "logs", "PCMBase_timing")
cd(time_dir)

sets = ["hiv", "prok", "mammals"]

nTraits = [3, 7, 8]
nTaxa = [1536, 705, 3649]

meths = ["pcm", "beast"]

n = length(sets)
p = length(meths)

paths = Array{String, 3}(undef, n, p, Nruns)

for i = 1:n
    for j = 1:p
        for k = 1:Nruns
            paths[i, j, k] = "$(sets[i])PCMComparisonTimer_$(meths[j])_$k.txt"
        end
    end
end

times = zeros(n, p, Nruns)
for i = 1:n
    for j = 1:p
        for k = 1:Nruns
            times[i, j, k] = get_time(paths[i, j, k], meths[j])
        end
    end
end

times = times / Nreps
speed_ups = times[:, 1, :] ./ times[:, 2, :]

mins = zeros(n)
maxs = zeros(n)
means = zeros(n)

for i = 1:n
    x = speed_ups[i, :]
    means[i] = mean(x)
    mins[i] = minimum(x)
    maxs[i] = maximum(x)
end

df = DataFrame(run = sets, N = nTaxa, P = nTraits, mean = means, min = mins, max = maxs)
CSV.write(joinpath(@__DIR__, "storage", "speed_results.csv"), df)
