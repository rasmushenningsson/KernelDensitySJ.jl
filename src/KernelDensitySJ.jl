module KernelDensitySJ

using Statistics
using DataStructures
using Roots

export bwsj, gaussiansmoothing, GaussianKernelSmoother

include("sums.jl")
include("phi.jl")
include("bwsj.jl")
include("smooth.jl")

end
