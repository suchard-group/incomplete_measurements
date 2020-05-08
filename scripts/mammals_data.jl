cd(@__DIR__)
push!(LOAD_PATH, pwd())

using DelimitedFiles, DataFrames, CSV

import BeastUtils.RTrees
trees = RTrees


const MISSING_VAL = -999.0
const PSS = 0.001


### Setting up directoreisand
root_dir = joinpath(pwd(), "..")
data_dir = joinpath(root_dir, "data")

newick_path = joinpath(data_dir, "mammals_newick.txt")
full_newick = read(newick_path, String)
tree = trees.parse_newick(full_newick)

data_path = joinpath(data_dir, "PanTHERIA_1-0_WR05_Aug2008.txt")


raw_data = readdlm(data_path)

### Need to combine 'Genus species' to 'Genus_species'
comb_cols = [5, 6]
n, p = size(raw_data)
n = n - 1 #Omit header row
p = p - 1 #Extra column due to the fact that the binomial name was split into two columns

all_data = Matrix{Any}(undef, n, p)
all_data[:, 1:4] .= raw_data[2:end, 1:4]
for i = 1:n
    all_data[i, 5] = join(raw_data[i + 1, comb_cols], '_')
end
all_data[:, 6:end] .= raw_data[2:end, 7:end]
header = raw_data[1, 1:p]

binomial_col = 5 # all_data column corresponding to binomial name
@show header[binomial_col]
cols = [7, 11, 16, 21, 22, 24, 33, 23, 11] # columns corresponding to values of interest
col_names = ["body_mass",
            "age_at_first_birth",
            "gestation_length",
            "litter_size",
            "litters_per_year",
            "neonate_body_mass",
            "weaning_age",
            "max_longevity",
            "age_at_first_birth"]

display([header[cols] col_names])


### Converting the missing values (-999.0) to 'missing' before working with the data

trimmed_data = all_data[:, cols]
for i = 1:length(trimmed_data)
    if trimmed_data[i] == MISSING_VAL
        trimmed_data[i] = missing
    end
end

col_names = ["body_mass", # reassigning col_names to reflect traits used in final analysis
            "age_at_first_birth",
            "gestation_length",
            "litter_size",
            "litters_per_year",
            "neonate_body_mass",
            "weaning_age",
            "reproductive_lifespan"]


taxa = all_data[:, binomial_col]
p = length(cols) - 1

m_mult = 365 / 12 # Average number of days in a month

data = Matrix{Union{Float64, Missing}}(undef, n, p)
data[:, 1:(p - 1)] .= trimmed_data[:, 1:(p - 1)]
data[:, p] .= trimmed_data[:, p] * m_mult .- trimmed_data[:, p + 1] #Longevity is in months, but age at first birth is in days

### Remove taxa from tree and data list for various reasons

# Remove rows with all missing values

keep_inds = zeros(Int, n)
current_ind = 0
for i = 1:n
    keep = false
    for j = 1:p
        if !ismissing(data[i, j])
            keep = true
            break
        end
    end
    if keep
        global current_ind += 1
        keep_inds[current_ind] = i
    end
end

keep_inds = keep_inds[1:current_ind]

# Need to manually remove one additional taxon (it has a reproductive_lifespan of 0.0)
rem_ind = findfirst(x -> x == "Monodelphis_dimidiata", taxa)
deleteat!(keep_inds, findfirst(x -> x == rem_ind, keep_inds))


taxa = taxa[keep_inds]
data = data[keep_inds, :]


# Remove taxa from the data matrix not represented in the tree
not_on_tree = setdiff(taxa, tree.tip_labels)
keep_inds = findall(x -> !(taxa[x] in not_on_tree), 1:length(taxa))

taxa = taxa[keep_inds]
data = data[keep_inds, :]
log_data = log10.(data)


# Remove taxa from the tree that aren't present in the data
not_in_data = setdiff(tree.tip_labels, taxa)
for taxon in not_in_data
    trees.trim_tree!(tree, taxon)
end

### Save newick and data

# newick
newick = trees.make_newick(tree)
write(joinpath(data_dir, "mammals_trimmed_newick.txt"), newick)


# trait data
df = DataFrame()

df[!, :taxon] = taxa
for i = 1:p
    df[!, Symbol(col_names[i])] = log_data[:, i]
end

CSV.write(joinpath(data_dir, "mammals_log_data.csv"), df)


### Additional newicks for PCMBase comparison

# Standardize tree and add prior sample size
max_distance = maximum([trees.distance_to_root(tree, i) for i = 1:tree.n_tips])
tree.edge_lengths ./= max_distance
pcm_newick_noRoot = trees.make_newick(tree)

trees.add_root!(tree, 1 / PSS)
pcm_newick = trees.make_newick(tree)

# Save newick and data


write(joinpath(data_dir, "mammals_trimmed_pcm_newick.txt"), pcm_newick)
write(joinpath(data_dir, "mammals_trimmed_pcm_newick_noRoot.txt"), pcm_newick_noRoot)




### Checking that the trimming process worked. This doesn't check topology (I did that by hand on smaller examples), but it does make

old_tree = trees.parse_newick(full_newick)
new_tree = trees.parse_newick(newick)


tol = 1e-10
for i in 1:tree.n_tips
    taxon = tree.tip_labels[i]
    d1 = trees.distance_to_root(old_tree, taxon)
    d2 = trees.distance_to_root(new_tree, taxon)
    if abs(d1 - d2) > tol
        error("The distances for taxon $taxon do not match.")
    end
end
