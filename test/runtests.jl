using KernelDensitySJ
using Distributions
using StableRNGs
using Test

# Ground truth was generated using the bw.SJ() R function, with nb=10000 and tol=1e-3.
# Some differences are expected because of how the double sums are approximated.

@testset "randn" begin
	rng = StableRNG(1234)
	# @info "A"
	X = randn(rng, 20)
	@test bwsj(X)≈0.6189722 rtol=0.1

	# @info "B"
	X = randn(rng, 100)
	@test bwsj(X)≈0.4314236 rtol=0.1

	# @info "C"
	X = randn(rng, 1000)
	@test bwsj(X)≈0.2823561 rtol=0.1
end

@testset "quantile" begin
	N = 100
	α = 2.0.^(-4:3)
	X = quantile.(Normal(), range(1/2N,stop=1-1/2N,length=N))

	@test bwsj(sign.(X).*abs.(X).^α[1])≈0.05301544  rtol=0.1
	@test bwsj(sign.(X).*abs.(X).^α[2])≈0.07877113  rtol=0.1
	@test bwsj(sign.(X).*abs.(X).^α[3])≈0.1279672   rtol=0.1
	@test bwsj(sign.(X).*abs.(X).^α[4])≈0.2092015   rtol=0.1
	@test bwsj(sign.(X).*abs.(X).^α[5])≈0.4747508   rtol=0.1
	@test bwsj(sign.(X).*abs.(X).^α[6])≈0.1394881   rtol=0.1
	@test bwsj(sign.(X).*abs.(X).^α[7])≈0.03132362  rtol=0.1
	@test bwsj(sign.(X).*abs.(X).^α[8])≈0.004705298 rtol=0.1 # NB: run with nb=10000000 to get accurate answer
end
