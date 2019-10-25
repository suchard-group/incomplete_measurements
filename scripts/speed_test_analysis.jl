using Revise

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

const Nreps = 1_000

time_dir = joinpath(@__DIR__, "..", "logs", "PCMBase_timing")
cd(time_dir)


mammals_beast_path = "mammalsPCMComparisonTimer_beast.txt"
mammals_pcm_path = "mammalsPCMComparisonTimer_pcm.txt"
hiv_beast_path = "hivPCMComparisonTimer_beast.txt"
hiv_pcm_path = "hivPCMComparisonTimer_pcm.txt"
#TODO: add prokaryotes

#all times are in seconds
m_pcm_time = pcm_time(mammals_pcm_path)
m_beast_time = beast_time(mammals_beast_path)
h_pcm_time = pcm_time(hiv_pcm_path)
h_beast_time = beast_time(hiv_beast_path)


data_sets = ["mammals", "mammals", "hiv", "hiv"]
software = ["beast", "pcm", "beast", "pcm"]
times = [m_beast_time, m_pcm_time, h_beast_time, h_pcm_time]
reps_per_sec = Nreps ./ times

mammals_increase = times[2] / times[1]
hiv_increase = times[4] / times[3]

@show mammals_increase
@show hiv_increase
