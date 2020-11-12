@testset "$f" for f in (gaussiansmoothing_reference,gaussiansmoothing)
	@testset "twopoints" begin
		x = Float64[0,1]
		y = Float64[0,1]
		xeval = Float64[0]
		@test f(x,y,1/√2,xeval) ≈ [1/(ℯ+1)]
		@test f(x,y,1,xeval) ≈ [1/(√ℯ+1)]
		@test f(x,y,√2,xeval) ≈ [1/(exp(1/4)+1)]

		@test f(x.+1,y,1,xeval.+1) ≈ [1/(√ℯ+1)] # translation
		@test f(x.*5,y,5,xeval) ≈ [1/(√ℯ+1)] # scaling
		@test f(-x,y,1,xeval) ≈ [1/(√ℯ+1)] # mirror

		@test f(x,y.+1,1,xeval) ≈ [1+1/(√ℯ+1)] # translation
		@test f(x,3y,1,xeval) ≈ [3/(√ℯ+1)] # scaling
		@test f(x,-y,1,xeval) ≈ -[1/(√ℯ+1)] # mirror

		# extreme bandwidths
		@test f(x,y,1e-5,[0.5]) ≈ [0.5]
		@test f(x,y,1e5,[0.5]) ≈ [0.5]
		@test f(x,y,1e-100,[0.5]) ≈ [0.5]
		@test f(x,y,1e100,[0.5]) ≈ [0.5]

		@test f(x,y,1e-5,[0.6]) ≈ [1.0]
		@test f(x,y,1e5,[0.6]) ≈ [0.5]
		@test f(x,y,1e-100,[0.6]) ≈ [1.0]
		@test f(x,y,1e100,[0.6]) ≈ [0.5]
	end

	@testset "fewpoints" begin
		x = [-4, 3.5, 10, 12]
		y = [1.0, -2.8, 8, 13]
		xeval = [0.0]
		bandwidths = 10.0.^(-3:3)
		ground_truth = [-2.79999999999999982236, -2.79999999999999982236, -2.79999999999999982236, -2.29473588724824868121, 3.20241290611236881292, 4.78316896715214792770, 4.79983163801801130898]
		@testset "bandwidth=$bw" for (bw,gt) in zip(bandwidths,ground_truth)
			@test f(x,y,bw,xeval) ≈ [gt] rtol=1e-3
			@test f(x.*-3.1.-7.8,y,bw.*3.1,xeval.-7.8) ≈ [gt] rtol=1e-3
		end
	end
end


@testset "npoints=N" for N in [99,100,101,995,1000,1005,9997,10000,10003]
	rng = StableRNG(42)
	x = randn(rng, N)
	y = randn(rng, N)
	xeval = randn(rng,10)

	bandwidths = 10.0.^ [0,-1,1,-2,2,-3,3,-5,5,-10,10,-100,100]
	@testset "bandwidth=$bw" for bw in bandwidths
		gt = gaussiansmoothing_reference(x,y,bw,xeval)
		@test gaussiansmoothing(x,y,bw,xeval;rtol=1e-3) ≈ gt rtol=1e-3
		@test gaussiansmoothing(x,y,bw,xeval;rtol=1e-6) ≈ gt rtol=1e-6
		@test gaussiansmoothing(x,y,bw,xeval;rtol=1e-9) ≈ gt rtol=1e-9
	end
end
