module KernelDensitySJ

using Statistics
using DataStructures
using Roots

export bwsj

include("sums.jl")
include("phi.jl")
include("reference.jl")
include("bwsj.jl")

end
