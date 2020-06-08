# Implementation of the "bwsj" method, from
# Sheather and Jones, "A Reliable Data-based Bandwidth Selection Method for Kernel Density Estimation"
# Notation below follows roughly the notation in the paper.



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
		npoints==0 && return zero(T),zero(T) # degenerate case
		x = tree.diagonalSums[d][k1]/(npoints*α)
		npoints.*ϕbounds(Val{N}(), 0.0, (X[i2]-X[i1])/α, x)
	else
		j1 = (k2-1)*intervalSize+1
		j2 = min(k2*intervalSize, n)
		ni = intervalSize
		nj = (j2-j1+1)
		npoints = ni*nj
		meanI = tree.intervalSums[d][k1]/(ni*α)
		meanJ = tree.intervalSums[d][k2]/(nj*α)
		npoints.*ϕbounds(Val{N}(), (X[j1]-X[i2])/α, (X[j2]-X[i1])/α, meanJ-meanI)
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


# Let s = ∑ᵢ∑ⱼϕ⁽ᴺ⁾(α⁻¹(Xᵢ-Xⱼ)), where the sum is taken over i<j.
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

# Let s = ∑ᵢ∑ⱼϕⁱᵛ(α⁻¹(Xᵢ-Xⱼ)), where the sum is taken over i<j.
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
	C -= n*ϕ4(0) # get rid of the diagonal entries
	C/=2 # because of symmetry, we can sum over j>i (upper triangular part), effectively halving the sum

	ssign(promote_type(Float64,T), α2, X, tree, C)
end

function SD(α, X::AbstractArray{T}, tree::SumTree; rtol=0.05) where T
	# α⁻⁵∑ᵢ∑ⱼϕⁱᵛ(α⁻¹(Xᵢ-Xⱼ))
	lb,ub = sumapprox(Val{4}(), promote_type(Float64,T), α, X, tree; rtol=rtol)
	s = (lb+ub) + length(X)*ϕ4(0) # 2*(lb+ub)/2 + diagonal entries
	s/α^5
end

function TD(b, X::AbstractArray{T}, tree::SumTree; rtol=0.05) where T
	# b⁻⁷∑ᵢ∑ⱼϕᵛⁱ(b⁻¹(Xᵢ-Xⱼ))
	lb,ub = sumapprox(Val{6}(), promote_type(Float64,T), b, X, tree; rtol=rtol)
	s = (lb+ub) + length(X)*ϕ6(0) # 2*(lb+ub)/2 + diagonal entries
	-s/b^7
end

"""
	bwsj(X; rtol=0.1, leafsize=10, lower=λ*n^(-1/5)*1e-9, upper=λ*n^(-1/5)*1e9)

Bandwidth selection for kernel density estimation using Gaussian kernels.
Based on the paper Sheather, S. J., & Jones, M. C. (1991). A reliable data‐based bandwidth selection method for kernel density estimation. Journal of the Royal Statistical Society: Series B (Methodological), 53(3), 683-690.

Solves the equation [ R(K) / (nσ⁴ₖS_D(α₂(h))) ]^(1/5) - h = 0 by root finding, bounding the objective such that each iteration can be terminated when the sign of the objective is known.

* `rtol` - the relative tolerance for the bandwidth.
* `leafsize` - Below this size, exact computations are used instead of bounds. Does not impact the accuracy of the result, only the execution time.
* `lower` - Bandwidth lower bound. There is usally no need to change the bounds.
* `upper` - Bandwidth upper bound.
"""
function bwsj(X; rtol=0.1, leafsize=10, lower=nothing, upper=nothing)
	issorted(X) || (X=sort(X))
	n = length(X)
	tree = SumTree(X,leafsize)

	q25,q75 = quantile(X,(0.25,0.75))
	λ = min(q75-q25, std(X)*1.349) # Givivng the same scale estimate as e.g. the KernSmooth R package.
	a = 0.920λ*n^(-1/7)
	b = 0.912λ*n^(-1/9)
	α2Constant = 1.357*(SD(a,X,tree;rtol=0.5rtol)/TD(b,X,tree;rtol=0.5rtol))^(1/7)

	# Decide if we need to deal with bounds better.
	# But note that execution is very fast for bad guesses for h.

	# The Bandwidth scales roughly with n^(-1/5)
	isnothing(lower) && (lower = λ*n^(-1/5)*1e-9)
	isnothing(upper) && (upper = λ*n^(-1/5)*1e9)

	find_zero(h->objectivesign(h,X,tree,α2Constant), (lower,upper), Bisection(), xrtol=rtol)
end

