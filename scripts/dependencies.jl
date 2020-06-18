required_packages = [
    "CSV",
    "DataFrames",
    "DelimitedFiles",
    "Distributions",
    "LightXML",
    "LinearAlgebra",
    "Statistics",
    "StatsBase"
]


import Pkg

for pkg in required_packages
    Pkg.add(pkg)
end

# installs BeastUtils package which stores lots of utilities for working with BEAST
Pkg.add(Pkg.PackageSpec(url="https://github.com/gabehassler/BeastUtils.jl.git", rev="9b5b58e0155daf15148464bb1dfada7a124e5043"))
