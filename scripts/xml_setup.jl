using Revise
using CSV, LightXML, DataFrames
using XMLConstructor2 #personal module


data_dir = joinpath(@__DIR__, "..", "data")
xml_dir = joinpath(@__DIR__, "..", "xml")

function make_xml(data_path::String, newick_path::String, xml_path::String,
                filename::String; dates_path::String = "")

    df = CSV.read(data_path)

    use_dates = false
    if length(dates_path) > 0
        dates_df = CSV.read(dates_path)
        @assert dates_df[!, :taxon] == df[!, :taxon]
        use_dates = true
    end

    newick = read(newick_path, String)

    taxa, data = XMLConstructor2.df_to_matrix(df)


    bx = XMLConstructor2.make_MBD_XML(data, taxa, newick, chain_length = 100_000)
    if use_dates
        XMLConstructor2.use_dates!(bx)
        bx.data_el.node_times = dates_df[!, :date]
    end
    bx.mcmc_el.screen_logEvery = 100
    bx.mcmc_el.file_logEvery = 10
    bx.mcmc_el.filename = filename
    XMLConstructor2.add_MBD_loggables!(bx)


    xdoc = XMLConstructor2.make_xml(bx)

    save_file(xdoc, xml_path)
    free(xdoc)
end

mammals_data_path = joinpath(data_dir, "mammals_log_data.csv")
mammals_newick_path = joinpath(data_dir, "mammals_trimmed_newick.txt")
mammals_xml_path = joinpath(xml_dir, "mammals2.xml") #TODO: change to just 'mammals.xml'

hiv_data_path = joinpath(data_dir, "hiv_processed_data.csv")
hiv_dates_path = joinpath(data_dir, "hiv_dates.csv")
hiv_newick_path = joinpath(data_dir, "hiv_newick.txt")
hiv_xml_path = joinpath(xml_dir, "hiv2.xml") #TODO: change to just 'hiv.xml'

make_xml(mammals_data_path, mammals_newick_path, mammals_xml_path, "mammals")
make_xml(hiv_data_path, hiv_newick_path, hiv_xml_path, "hiv", dates_path = hiv_dates_path)
