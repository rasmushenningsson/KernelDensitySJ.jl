module KernelDensitySJ

using Statistics
using DataStructures
using Roots

export bwsj

include("sums.jl")
include("approximations.jl")
include("sj.jl")

end
