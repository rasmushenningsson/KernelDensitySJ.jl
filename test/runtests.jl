using KernelDensitySJ
using Distributions
using StableRNGs
using Test

import KernelDensitySJ: ϕ4,ϕ6,ϕ4bounds,ϕ6bounds

include("reference.jl")

@testset "Bandwidth" begin
	include("test_bw.jl")
end
@testset "Smoothing" begin
	include("test_smooth.jl")
end
