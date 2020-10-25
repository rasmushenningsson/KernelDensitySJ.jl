module KernelDensitySJ

using Statistics
using DataStructures
using Roots

export bwsj

include("sums.jl")
include("phi.jl")
include("bwsj.jl")
include("smooth.jl")
include("reference.jl")

end
