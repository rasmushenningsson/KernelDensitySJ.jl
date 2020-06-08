using KernelDensitySJ
using Distributions
using StableRNGs
using Test

# Ground truth was generated using the bw.SJ() R function, with nb=10000 and tol=1e-3 unless further specified.
# Some differences are expected because of how the double sums are approximated.

@testset "range" begin
	# R: nb=10000000, tol=1e-6
	@test bwsj(range(0,stop=1,length= 96); rtol=1e-3)≈0.1109613 rtol=1e-2
	@test bwsj(range(0,stop=1,length= 97); rtol=1e-3)≈0.1105281 rtol=1e-2
	@test bwsj(range(0,stop=1,length= 98); rtol=1e-3)≈0.1101011 rtol=1e-2
	@test bwsj(range(0,stop=1,length= 99); rtol=1e-3)≈0.1096801 rtol=1e-2
	@test bwsj(range(0,stop=1,length=100); rtol=1e-3)≈0.1092651 rtol=1e-2
	@test bwsj(range(0,stop=1,length=101); rtol=1e-3)≈0.1088558 rtol=1e-2
	@test bwsj(range(0,stop=1,length=102); rtol=1e-3)≈0.1084521 rtol=1e-2
	@test bwsj(range(0,stop=1,length=103); rtol=1e-3)≈0.1080540 rtol=1e-2
	@test bwsj(range(0,stop=1,length=104); rtol=1e-3)≈0.1076612 rtol=1e-2
	@test bwsj(range(0,stop=1,length=105); rtol=1e-3)≈0.1072737 rtol=1e-2
end

@testset "scale" begin
	@test bwsj(range(0,stop=1e3,  length=101))≈0.1088518*1e3   rtol=0.1
	@test bwsj(range(0,stop=1e-3, length=101))≈0.1088518*1e-3  rtol=0.1
	@test bwsj(range(0,stop=1e24, length=101))≈0.1088518*1e24  rtol=0.1
	@test bwsj(range(0,stop=1e-24,length=101))≈0.1088518*1e-24 rtol=0.1
end

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

@testset "outliers" begin
	# ground truth generated with bwsj_reference(; rtol=1e-24)
	rng = StableRNG(2017)
	X = rand(rng, 100)

	Y1 = 1e9 .+ randn(rng,1)
	@test bwsj(vcat(X,Y1))≈0.11706501880745836 rtol=0.1

	Y2 = 1e9 .+ randn(rng,10)
	@test bwsj(vcat(X,Y2))≈0.1197213803664104 rtol=0.1

	Y3 = 1e9.*sign.(randn(rng,10)) .+ randn(rng,10)
	@test bwsj(vcat(X,Y3))≈0.12247926574875453 rtol=0.1
end
