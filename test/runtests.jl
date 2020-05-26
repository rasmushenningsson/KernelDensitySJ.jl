using KernelDensitySJ
using Distributions
using StableRNGs
using Test

# Ground truth was generated using the bw.SJ() R function, with nb=10000 and tol=1e-3.
# Some differences are expected because of how the double sums are approximated.

@testset "randn" begin
	rng = StableRNG(1234)
	X = randn(rng, 20)
	@test bwsj(X)≈0.6189722 rtol=0.1

	X = randn(rng, 100)
	@test bwsj(X)≈0.4314236 rtol=0.1

	X = randn(rng, 1000)
	@test bwsj(X)≈0.2823561 rtol=0.1
end

@testset "quantile" begin
	N = 100
	expected = [0.05301544, 0.1279672, 0.4747508, 0.03132362, 5.242246e-05, 3.903804e-13, 4.232559e-44]
	for (i,α) in enumerate(4.0.^(-2:4))
		# NB, the laste three tests do not yet pass due to precision issues
		X = quantile.(Normal(), range(1/2N,stop=1-1/2N,length=N))
		X = sign.(X) .* abs.(X).^α
		@test bwsj(X)≈expected[i] rtol=0.1
	end
end
