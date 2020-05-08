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
    unit_mult = Dict("seconds" => 1, "minutes" => 60, "hours" => 3600)
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
const Nruns = 10

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

times = Nreps ./ times
speed_ups = times[:, 2, :] ./ times[:, 1, :]

min_speedup = zeros(n)
max_speedup = zeros(n)
mean_speedup = zeros(n)

min_beast = zeros(n)
max_beast = zeros(n)
mean_beast = zeros(n)

min_pcm = zeros(n)
max_pcm = zeros(n)
mean_pcm = zeros(n)

for i = 1:n
    su = speed_ups[i, :]
    mean_speedup[i] = mean(su)
    min_speedup[i] = minimum(su)
    max_speedup[i] = maximum(su)

    beast = times[i, 2, :]
    mean_beast[i] = mean(beast)
    min_beast[i] = minimum(beast)
    max_beast[i] = maximum(beast)

    pcm = times[i, 1, :]
    mean_pcm[i] = mean(pcm)
    min_pcm[i] = minimum(pcm)
    max_pcm[i] = maximum(pcm)

end

df = DataFrame(run = sets, N = nTaxa, P = nTraits,
            meanSpeedup = mean_speedup,
            minSpeedup = min_speedup,
            maxSpeedup = max_speedup,
            meanBeast = mean_beast,
            minBeast = min_beast,
            maxBeast = max_beast,
            meanPCM = mean_pcm,
            minPCM = min_pcm,
            maxPCM = max_pcm)


### Sim speed study


cd(@__DIR__)
filenames = readlines(joinpath("storage", "PCMBase_comparison",
                                "sim_timing_files.txt"))

n = length(filenames)

Ns = zeros(Int, n)
Ps = zeros(Int, n)
Rs = zeros(Int, n)
b_times = zeros(n)
p_times = zeros(n)

for i = 1:n
    N = parse(Int, match(r"_(\d+)Taxa", filenames[i])[1])
    Ns[i] = N

    P = parse(Int, match(r"_(\d+)Traits", filenames[i])[1])
    Ps[i] = P

    rep = parse(Int, match(r"Traits_(\d+)", filenames[i])[1])
    Rs[i] = rep

    # if N == 10000 && P == 20
    #     println("!!!!!!!!!!!!!!!!!!!REMOVE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    #     continue
    # end

    pcm_path = joinpath(time_dir, "$(filenames[i])Timer_pcm.txt")
    beast_path = joinpath(time_dir, "$(filenames[i])Timer_beast.txt")
    b_times[i] = get_time(beast_path, "beast")
    p_times[i] = get_time(pcm_path, "pcm")


end

b_rates = Nreps ./ b_times
p_rates = Nreps ./ p_times

sus = b_rates ./ p_rates

NPs = [(Ns[i], Ps[i]) for i = 1:n]

unique_NPs = unique(NPs)

n = length(unique_NPs)
b_maxRates = zeros(n)
b_minRates = zeros(n)
b_meanRates = zeros(n)

p_maxRates = zeros(n)
p_minRates = zeros(n)
p_meanRates = zeros(n)

s_max = zeros(n)
s_min = zeros(n)
s_mean = zeros(n)



for i = 1:n

    inds = findall(x -> x == unique_NPs[i], NPs)

    b_reps = b_rates[inds]
    p_reps = p_rates[inds]
    su_reps = sus[inds]

    b_minRates[i] = minimum(b_reps)
    b_maxRates[i] = maximum(b_reps)
    b_meanRates[i] = mean(b_reps)

    p_minRates[i] = minimum(p_reps)
    p_maxRates[i] = maximum(p_reps)
    p_meanRates[i] = mean(p_reps)

    s_min[i] = minimum(su_reps)
    s_max[i] = maximum(su_reps)
    s_mean[i] = mean(su_reps)

end

df_sim = DataFrame(run = fill("sim", n),
                N = [unique_NPs[i][1] for i = 1:n],
                P = [unique_NPs[i][2] for i = 1:n],
                meanSpeedup = s_mean,
                minSpeedup = s_min,
                maxSpeedup = s_max,
                meanBeast = b_meanRates,
                minBeast = b_minRates,
                maxBeast = b_maxRates,
                meanPCM = p_meanRates,
                minPCM = p_minRates,
                maxPCM = p_maxRates
                )

append!(df, df_sim)

CSV.write(joinpath(@__DIR__, "storage", "speed_results.csv"), df)
