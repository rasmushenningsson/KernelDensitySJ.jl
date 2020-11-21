module KernelDensitySJ

using Statistics
using DataStructures
using Roots

export bwsj, smooth, GaussianKernelSmoother

include("sums.jl")
include("phi.jl")
include("bwsj.jl")
include("smooth.jl")

end
