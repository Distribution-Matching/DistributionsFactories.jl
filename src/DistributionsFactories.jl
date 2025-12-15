module DistributionsFactories

using Distributions
using SpecialFunctions
using Roots


include("exists_unique_dist_from_mean_var.jl")
include("dist_from_mean_var.jl")

export dist_from_mean_var,
        exists_unique_dist_from_mean_var

end # end the module