using LightXML
using BeastUtils.XMLTools

const REPEATS = 1
const OPERATORNAME = "newLatentLiabilityGibbsOperator"

function format_xml(filenames::Array{String},
        xml_dir::String,
        chain_lengths::Array{Int},
        log_everys::Array{Int},
        weights::Array{Int};
        use_weights::Bool = true
        )

    n = length(filenames)
    for i = 1:n
        xml_filename = "$(filenames[i]).xml"
        xdoc = parse_file(joinpath(xml_dir, xml_filename))
        replace_chainLength!(xdoc, chain_lengths[i])
        replace_logEvery!(xdoc, log_everys[i])
        if use_weights
            set_operator_weight!(xdoc, OPERATORNAME, weights[i])
        end
        for j = 1:REPEATS
            new_name = "$(filenames[i])_r$j"
            new_logname = "$new_name.log"
            new_xmlname = "$new_name.xml"
            replace_filelog!(xdoc, new_logname)
            save_file(xdoc, joinpath(xml_dir, new_xmlname))
        end
        free(xdoc)
    end
end

slow_filenames = ["hiv_newTiming", "hiv_newTiming_rm", "mammals_newTiming",
    "mammals_newTiming_rm"]

fast_filenames = ["hiv_diff", "hiv_res", "mammals_diff", "mammals_res"]

slow_weights = [10, 50, 500, 50]
slow_logEverys = [100, 1000, 1000, 5000]
slow_chainLengths = [1000000, 50000000, 25000000, 50000000]
# slow_chainLengths = [1000, 1000, 1000, 1000]

fast_logEverys = [1, 25, 1, 50]
fast_chainLengths = [100000, 100000, 100000, 100000]
# fast_chainLengths = [1000, 1000, 1000, 1000]


n = 4
m = 2



filenames = [["$(slow_filenames[i])_r$j", "$(fast_filenames[i])_r$j"]
        for i = 1:n, j = 1:REPEATS]


#scramble the order
for i = 1:(n * REPEATS)
    nms = copy(filenames[i])
    x = rand([0, 1])
    dims = [1, 2]
    if x == 1
        dims = [2, 1]
    end
    filenames[i] = nms[dims]
end


xml_dir = joinpath(@__DIR__, "..", "xml", "timing")


format_xml(slow_filenames, xml_dir, slow_chainLengths, slow_logEverys,
    slow_weights, use_weights = true)

format_xml(fast_filenames, xml_dir, fast_chainLengths, fast_logEverys, [0],
    use_weights = false)


### Below code used to setting up files to submit to cluster

# batch_names = ["hiv_diff", "hiv_res", "mammals_diff", "mammals_res"]
#
# batch_names = ["$(x)_r$j.txt" for x in batch_names, j = 1:REPEATS]
#
#
# using BatchSetup
# setup_sh(Directories.working_batch, filenames, batch_names,
#     source_dir = "$(BatchSetup.SOURCEDIR)/timing", tmp_dir = true)
