using BeastUtils.DataStorage, BeastUtils.MatrixUtils, BeastUtils.XMLTools

using Statistics, LightXML, StatsBase

function setup_xml(format_xmlpath::String,
    output_xmlpath::String,
    full_data::Matrix{Float64},
    removed_data::Matrix{Float64},
    taxa::Vector{String},
    dates::Vector{Float64},
    percs::Vector{Float64};
    trait_name::String = "traits")

    xdoc = parse_file(format_xmlpath)
    setup_xml(xdoc, output_xmlpath, replace_data,
        keep_dims, taxa, dates, percs, trait_name = trait_name)
    free(xdoc)

end

function setup_xml(xdoc::XMLDocument,
    output_xmlpath::String,
    full_data::Matrix{Float64},
    removed_data::Matrix{Float64},
    taxa::Vector{String},
    dates::Vector{Float64},
    percs::Vector{Float64};
    trait_name::String = "traits")


    trait_names = [trait_name, "$(trait_name)True"]
    replace_all_data!(xdoc, taxa, [removed_data, full_data], trait_names, dates = dates)
    # setup_trait_validation!(xdoc)
    save_file(xdoc, output_xmlpath)
end


function cull_data(data::Matrix{Float64}, percs::Vector{Float64})

    @assert minimum(percs) >= 0.0 && maximum(percs) < 1.0

    n, p = size(data)

    n_remove = zeros(Int, p)
    for i = 1:p
        m = Int(round(percs[i] * n))
        n_remove[i] = m
    end

    new_data = copy(data)

    for i = 1:p
        removed = sample(1:n, n_remove[i], replace = false)

        new_data[removed, i] .= NaN

    end
    return new_data
end



function make_xml(data_path::String,
                    dates_path::String,
                    newick_path::String,
                    remove_percs::Vector{Float64},
                    template_names::Vector{String},
                    xml_path::String,
                    n_reps::Int)


    taxa, data = DataStorage.csv_to_data(data_path)
    taxa2, dates = DataStorage.csv_to_data(dates_path)
    dates = vec(dates)

    @assert taxa == taxa2

    MatrixUtils.standardize_data!(data)

    newick = read(newick_path, String)


    template_paths = [joinpath(xml_path, "$i.xml") for i in template_names]
    n = length(template_names)
    new_names = Matrix{String}(undef, n_reps, n)
    for j = 1:n_reps

        culled_data = cull_data(data, remove_percs)

        for i = 1:n

            culled_data_i = culled_data[:, keep_dims[i]]
            data_i = data[:, keep_dims[i]]
            percs_i = remove_percs[keep_dims[i]]

            xdoc = parse_file(template_paths[i])
            replace_chainLength!(xdoc, CHAINLENGTH)
            replace_logEvery!(xdoc, LOGEVERY)

            replace_newick!(xdoc, newick)




            new_name = "$(template_names[i])_$j"
            @show new_name
            new_names[j, i] = new_name
            replace_filelog!(xdoc, "$new_name.log")

            new_path = joinpath(xml_path, "$new_name.xml")

            setup_xml(xdoc, new_path, data_i, culled_data_i,
                taxa, dates, percs_i, trait_name = "trait")

            free(xdoc)

        end

    end


end


CHAINLENGTH = 100000
LOGEVERY = 100


original_names = ["hiv_mse_d", "hiv_mse_dr", "hiv_mse_du", "hiv_mse_dru"]

keep_dims = [[1, 2], [1, 2], [2], [2]]


data_dir = joinpath(@__DIR__, "..", "data")
hiv_data_path = joinpath(data_dir, "hiv_processed_data.csv")
hiv_dates_path = joinpath(data_dir, "hiv_dates.csv")
hiv_newick_path = joinpath(data_dir, "hiv_newick.txt")

xml_dir = joinpath(@__DIR__, "..", "xml", "hiv_prediction")


remove_percs = [0.0, 0.5, 0.0]

p = 20

make_xml(hiv_data_path, hiv_dates_path, hiv_newick_path, remove_percs, original_names, xml_dir, p)

# Code below is for creating scrip to run xml on a cluster

# setup_sh(Directories.working_batch, new_names,
#     run_time = "01:00:00",
#     source_dir = "$(BatchSetup.SOURCEDIR)/many_traits/hiv_mse",
#     email = false)
