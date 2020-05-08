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
Pkg.add(PackageSpec(url="https://github.com/gabehassler/BeastUtils.jl.git"))
