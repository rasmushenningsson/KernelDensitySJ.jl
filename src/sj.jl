# Implementation of the "bw.sj" method, from
# Sheather and Jones, "A Reliable Data-based Bandwidth Selection Method for Kernel Density Estimation"
# Notation below follows roughly the notation in the paper.


ϕ(x) = 1/√(2π) * exp(-x^2/2)

ϕ4(x) = (x2=x*x; (x2*x2 - 6x2 + 3) * ϕ(x))
# ϕ5(x) = (x2=x*x; x4=x2*x2; -x*(x4 - 10x2 + 15) * ϕ(x)) # For reference: used to find local extrema of ϕ4
ϕ6(x) = (x2=x*x; x4=x2*x2; (x4*x2 - 15x4 + 45x2 - 15) * ϕ(x))
# ϕ7(x) = (x2=x*x; x4=x2*x2; x6=x4*x2; x*(x6 - 21x4 + 105x2 - 105)*ϕ(x))                 # For reference: used to find local extrema of ϕ6
# ϕ8(x) = (x2=x*x; x4=x2*x2; x6=x4*x2; x8=x4*x4; (x8 - 28x6 + 210x4 - 420x2 + 105)*ϕ(x)) # For reference: used to find intervals where ϕ6 is convex or concave.

ϕ(::Val{4},x) = ϕ4(x)
ϕ(::Val{6},x) = ϕ6(x)


function ϕ4extrema(a,b)
	# (x-√(5+√10))(x+√(5+√10))(x-√(5-√10))(x+√(5-√10))x
	# -√(5+√BigFloat(10)), √(5-√BigFloat(10)), 0, √(5-√BigFloat(10)), √(5+√BigFloat(10))
	# xExt = (-2.8569700138728056, 1.355626179974266, 0.0, 1.355626179974266, 2.8569700138728056)
	# yExt = (0.13912184973492348, -0.7399861849949221, 1.1968268412042982, -0.7399861849949221, 0.13912184973492348)
	# we only need to consider x>0 since the function is even
	xExt = (1.355626179974266, 2.8569700138728056)
	yExt = (-0.7399861849949221, 0.13912184973492348)

	# the implementation assumes 0<=a<=b
	if b<=xExt[1] || a>=xExt[2]
		ϕ4(b),ϕ4(a) # decreasing function in these intervals
	elseif a>=xExt[1]
		if b<=xExt[2]
			ϕ4(a),ϕ4(b) # increasing function in this interval
		else
			min(ϕ4(a),ϕ4(b)),yExt[2] # the function covers a local maximum
		end
	elseif b<=xExt[2]
		yExt[1],max(ϕ4(a),ϕ4(b)) # the function covers the global minimum
	else
		yExt[1], max(ϕ4(a),yExt[2]) # the function covers the global minimum and a local maximum
	end
end


function ϕ6extrema(a,b)
	# ϕ6 has global minimum at 0, global maximum at x≈1.15, two more local extrema and ϕ6(x)→0 as x→∞.
	xExt = (1.154405394739968, 2.366759410734541, 3.750439717725742)
	yExt = (4.2406124820663305, -1.4018233232694841, 0.15246889807941977)
	# the implementation assumes 0<=a<=b


	if xExt[3]<=a || (xExt[1]<=a && b<=xExt[2])
		ϕ6(b),ϕ6(a) # decreasing function in these intervals
	elseif b<=xExt[1] || (xExt[2]<=a && b<=xExt[3])
		ϕ6(a),ϕ6(b) # increasing function in these intervals
	elseif b<=xExt[2]
		min(ϕ6(a),ϕ6(b)),yExt[1] # covers global maximum
	elseif xExt[1]<=a && b<=xExt[3]
		yExt[2],max(ϕ6(a),ϕ6(b)) # covers local minimum
	elseif xExt[2]<=a
		min(ϕ6(a),ϕ6(b)),yExt[3] # covers second local maximum
	elseif xExt[1]<=a
		yExt[2],max(ϕ6(a),yExt[3])  # covers the local minimum and the second local maximum
	# elseif b<=xExt[3] # same as last case
	# 	min(ϕ6(a),yExt[2]),yExt[1] # covers global maximum and local minimum
	else
		min(ϕ6(a),yExt[2]),yExt[1] # covers global maximum and local minimum
	end
end



ϕextrema(::Val{4},a,b) = ϕ4extrema(a,b)
ϕextrema(::Val{6},a,b) = ϕ6extrema(a,b)



#---------------------------------------------------------------------------------------------------

struct Block
	depth::Int
	k1::Int
	k2::Int
	lb::Float64
	ub::Float64
end
Base.isless(b1::Block,b2::Block) = isless(b1.ub-b1.lb, b2.ub-b2.lb)


function ϕsum(::Val{N}, ::Type{T}, k1, k2, α, X, leafSize) where {N,T}
	n = length(X)

	i1 = (k1-1)*leafSize+1
	i2 = min(k1*leafSize, n)

	s = zero(T)
	if k1==k2
		@inbounds for i in i1:i2
			for j in i+1:i2 # upper triangular part
				s += ϕ(Val{N}(), (X[j]-X[i])/α)
			end
		end
	else
		j1 = (k2-1)*leafSize+1
		j2 = min(k2*leafSize, n)

		@inbounds for i in i1:i2
			for j in j1:j2
				s += ϕ(Val{N}(), (X[j]-X[i])/α)
			end
		end
	end
	s # exact
end


function bounds(::Val{N}, ::Type{T}, d, k1, k2, α, X, tree::SumTree)::Tuple{T,T} where {N,T}
	n = length(X)
	intervalSize = tree.leafSize * 2^(length(tree.intervalSums)-d)

	i1 = (k1-1)*intervalSize+1
	i2 = min(k1*intervalSize, n)

	if k1==k2 # block on the diagonal
		npoints = div((i2-i1+1)*(i2-i1),2)
		x = (tree.incSums[d][k1] - tree.decSums[d][k1])/(npoints*α)
		ϕbounds(Val{N}(), 0.0, (X[i2]-X[i1])/α, x, npoints)
	else
		j1 = (k2-1)*intervalSize+1
		j2 = min(k2*intervalSize, n)
		ni = intervalSize
		nj = (j2-j1+1)
		npoints = ni*nj
		meanI = tree.intervalSums[d][k1]/(ni*α)
		meanJ = tree.intervalSums[d][k2]/(nj*α)
		ϕbounds(Val{N}(), (X[j1]-X[i2])/α, (X[j2]-X[i1])/α, meanJ-meanI, npoints)
	end
end


function ssumstep!(heap::BinaryMaxHeap{Block}, ::Val{N}, ::Type{T},α,X,tree::SumTree)::Tuple{T,T} where {N,T}
	block = pop!(heap)

	# remove existing bounds
	lb = -block.lb
	ub = -block.ub

	if block.depth >= depth(tree) # Fallback to exact sum
		s = ϕsum(Val{N}(), T, block.k1, block.k2, α, X, tree.leafSize)
		lb += s
		ub += s
	else # create child blocks
		d = block.depth+1
		k1 = 2block.k1-1 # first child in first direction
		k2 = 2block.k2-1 # first child in second direction
		nbrBlocks = length(tree.intervalSums[d])

		lb11,ub11 = bounds(Val{N}(), T, d, k1, k2, α, X, tree)
		lb += lb11
		ub += ub11
		lb11!=ub11 && push!(heap, Block(d, k1, k2, lb11, ub11))

		if k1<k2<=nbrBlocks # only care about upper triangular part
			lb21,ub21 = bounds(Val{N}(), T, d, k1+1, k2, α, X, tree)
			lb += lb21
			ub += ub21
			lb21!=ub21 && push!(heap, Block(d, k1+1, k2, lb21, ub21))
		end

		if k2<nbrBlocks
			lb12,ub12 = bounds(Val{N}(), T, d, k1, k2+1, α, X, tree)
			lb += lb12
			ub += ub12
			lb12!=ub12 && push!(heap, Block(d, k1, k2+1, lb12, ub12))

			lb22,ub22 = bounds(Val{N}(), T, d, k1+1, k2+1, α, X, tree)
			lb += lb22
			ub += ub22
			lb22!=ub22 && push!(heap, Block(d, k1+1, k2+1, lb22, ub22))
		end
	end
	lb,ub
end


# Let s = ∑ᵢ∑ⱼϕ⁽ᴺ⁾(α⁻¹(Xᵢ-Xⱼ)).
# returns lower and upper bounds for s, terminating when the tolerance is good enough
function sumapprox(::Val{N},::Type{T},α,X,tree::SumTree;rtol)::Tuple{T,T} where {N,T}
	heap = BinaryMaxHeap{Block}()
	lb,ub = bounds(Val{N}(),T,1,1,1,α,X,tree)
	push!(heap,Block(1,1,1,lb,ub))
	while !isempty(heap) && !isapprox(lb,ub;rtol=rtol)
		lb,ub = (lb,ub) .+ ssumstep!(heap,Val{N}(),T,α,X,tree)
	end
	lb,ub
end

# Let s = ∑ᵢ∑ⱼϕⁱᵛ(α⁻¹(Xᵢ-Xⱼ)).
# returns sign(s-C)
function ssign(::Type{T},α,X,tree::SumTree,C)::Int where T
	heap = BinaryMaxHeap{Block}()
	lb,ub = bounds(Val{4}(),T,1,1,1,α,X,tree)
	push!(heap,Block(1,1,1,lb,ub))
	while !isempty(heap) && lb<=C<=ub
		lb,ub = (lb,ub) .+ ssumstep!(heap,Val{4}(),T,α,X,tree)
	end
	lb > C && return 1
	ub < C && return -1
	@warn "Numerical issues for upper/lower bounds, ub=$ub, lb=$lb, C=$C."
	return 0
end


# Computes the sign of the objective function by gradually approximating using lower/upper bounds
function objectivesign(h, X::AbstractArray{T}, tree::SumTree, α2Constant) where T
	n = length(X)
	α2 = α2Constant*h^(5/7) # 1.357[SD(a)/TD(b)]^(1/7) * h^(5/7)

	C = (n-1)*α2Constant^5*h^(-10/7)/(2*√π)
	C -= n*1.1968268412042982 # n*ϕ4(0) - get rid of the diagonal entries
	C/=2 # because of symmetry, we can sum over j>i (upper triangular part), effectively halving the sum

	ssign(promote_type(Float64,T), α2, X, tree, C)
end

function SD_bounded(α, X::AbstractArray{T}, tree::SumTree; rtol=0.05) where T
	# α⁻⁵∑ᵢ∑ⱼϕⁱᵛ(α⁻¹(Xᵢ-Xⱼ))
	n = length(X)
	lb,ub = sumapprox(Val{4}(), promote_type(Float64,T), α, X, tree; rtol=rtol)
	s = (lb+ub) + n*ϕ4(0) # 2*(lb+ub)/2 + diagonal entries
	s/α^5
end

function TD_bounded(b, X::AbstractArray{T}, tree::SumTree; rtol=0.05) where T
	# b⁻⁷∑ᵢ∑ⱼϕᵛⁱ(b⁻¹(Xᵢ-Xⱼ))
	n = length(X)
	lb,ub = sumapprox(Val{6}(), promote_type(Float64,T), b, X, tree; rtol=rtol)
	s = (lb+ub) + n*ϕ6(0) # 2*(lb+ub)/2 + diagonal entries
	-s/b^7
end


function _bwsj(X; rtol=0.1)
	n = length(X)
	tree = SumTree(X,10)

	q25,q75 = quantile(X,(0.25,0.75))
	λ = min(q75-q25, std(X)*1.349) # Givivng the same scale estimate as e.g. the KernSmooth R package.
	a = 0.920λ*n^(-1/7)
	b = 0.912λ*n^(-1/9)
	# α2Constant = 1.357*(SD(a,X)/TD(b,X))^(1/7)
	α2Constant = 1.357*(SD_bounded(a,X,tree;rtol=0.5rtol)/TD_bounded(b,X,tree;rtol=0.5rtol))^(1/7)


	# Decide if we need to deal with bounds better.
	# But note that execution is very fast for bad guesses for h.

	# The Bandwidth scales roughly with n^(-1/5)
	lower = λ*n^(-1/5)*1e-9
	upper = λ*n^(-1/5)*1e9

	find_zero(h->objectivesign(h,X,tree,α2Constant), (lower,upper), Bisection(), xrtol=rtol)
end


bwsj(X; kwargs...) = _bwsj(issorted(X) ? X : sort(X); kwargs...)



#---------------------------------------------------------------------------------------------------


# Unoptimized version for reference
function SD_reference(α, X::AbstractArray{T}) where T
	# (n(n-1))⁻¹α⁻⁵∑ᵢ∑ⱼϕⁱᵛ(α⁻¹(Xᵢ-Xⱼ))
	n = length(X)
	C = 1/(n*(n-1)*α^5)
	Y = X./α
	s = zero(T)
	@inbounds for i=1:n
		for j=i+1:n
			s += 2*ϕ4(Y[i]-Y[j])
		end
	end
	s += n*ϕ4(0)
	C * s
end

function TD_reference(b, X::AbstractArray{T}) where T
	# -(n(n-1))⁻¹b⁻⁷∑ᵢ∑ⱼϕᵛⁱ(b⁻¹(Xᵢ-Xⱼ))
	n = length(X)
	C = -1/(n*(n-1)*b^7)
	Y = X./b
	s = zero(T)
	@inbounds for i=1:n
		for j=i+1:n
			s += 2*ϕ6(Y[i]-Y[j])
		end
	end
	s += n*ϕ6(0)
	C * s
end


function objective_reference(h, X, α2Constant)
	# [ R(K) / (nσ⁴ₖSD(α₂(h))) ]^(1/5) - h = 0
	# R(K) = 1/(2√π)
	# σ⁴ₖ = 1
	n = length(X)
	α2 = α2Constant*h^(5/7) # 1.357[SD(a)/TD(b)]^(1/7) * h^(5/7)
	# denom = (n*2*√π) * SD(α2,X)
	# denom^(-1/5) - h
	1/(2*√(π)*n*h^5) - SD_reference(α2,X)
end


function _bwsj_reference(X; rtol=0.1)
	n = length(X)
	q25,q75 = quantile(X,(0.25,0.75))
	λ = min(q75-q25, std(X)*1.349) # Givivng the same scale estimate as e.g. the KernSmooth R package.
	a = 0.920λ*n^(-1/7)
	b = 0.912λ*n^(-1/9)
	α2Constant = 1.357*(SD_reference(a,X)/TD_reference(b,X))^(1/7)

	# Figure out lower and upper bounds for search.
	# The Bandwidth scales roughly with n^(-1/5). TODO: Use this fact!

	lower = 0.1*λ
	for i=1:11
		i==11 && error("Failed to find lower bound")
		objective(lower,X,α2Constant)>0 && break
		lower/=10
	end
	upper = 1*λ
	for i=1:11
		i==11 && error("Failed to find upper bound")
		objective(upper,X,α2Constant)<0 && break
		upper*=10
	end

	find_zero(h->objective_reference(h,X,α2Constant), (lower,upper), Bisection(), xrtol=rtol)
end


bwsj_reference(X; kwargs...) = _bwsj_reference(issorted(X) ? X : sort(X); kwargs...)

