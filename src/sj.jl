# initial tests for implementation of the "bw.sj" method, from 
# Sheather and Jones, "A Reliable Data-based Bandwidth Selection Method for Kernel Density Estimation"
# Notation below follows roughly the notation in the paper.



ϕ(x) = 1/√(2π) * exp(-x^2/2)

ϕ4(x) = (x2=x*x; (x2*x2 - 6x2 + 3) * ϕ(x))
ϕ5(x) = (x2=x*x; x4=x2*x2; -x*(x4 - 10x2 + 15) * ϕ(x))
ϕ6(x) = (x2=x*x; x4=x2*x2; (x4*x2 - 15x4 + 45x2 - 15) * ϕ(x))
ϕ7(x) = (x2=x*x; x4=x2*x2; x6=x4*x2; x*(x6 - 21x4 + 105x2 - 105)*ϕ(x))


ϕ(::Val{4},x) = ϕ4(x)
ϕ(::Val{5},x) = ϕ5(x)
ϕ(::Val{6},x) = ϕ6(x)
ϕ(::Val{7},x) = ϕ7(x)


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
			min(ϕ4(a),ϕ4(b)),yExt[2] # the function covers one local extrema (maxima)
		end
	elseif b<=xExt[2]
		yExt[1],max(ϕ4(a),ϕ4(b)) # the function covers one local extrema (minima)
	else
		yExt[1], max(ϕ4(a),yExt[2]) # the function covers two local extrema
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
		yExt[2],max(ϕ6(a),ϕ6(b)) # covers local minima
	elseif xExt[2]<=a
		min(ϕ6(a),ϕ6(b)),yExt[3] # covers second local maxima
	elseif xExt[1]<=a
		yExt[2],max(ϕ6(a),yExt[3])  # covers the local minima and the second local maxima
	# elseif b<=xExt[3] # same as last case
	# 	min(ϕ6(a),yExt[2]),yExt[1] # covers global maximum and local minima
	else
		min(ϕ6(a),yExt[2]),yExt[1] # covers global maxima and local minima
	end
end



ϕextrema(::Val{4},a,b) = ϕ4extrema(a,b)
ϕextrema(::Val{6},a,b) = ϕ6extrema(a,b)


#---------------------------------------------------------------------------------------------------

function leafsums(::Type{T}, X, leafSize::Int) where T
	n = length(X)
	nbrRanges = div(n-1,leafSize)+1
	intervalSums = Vector{T}(undef, nbrRanges)
	incSums      = Vector{T}(undef, nbrRanges)
	decSums      = Vector{T}(undef, nbrRanges)

	for k in 1:nbrRanges
		offs = (k-1)*leafSize
		nk = k<nbrRanges ? leafSize : mod1(n,leafSize)

		s = zero(T)
		sInc = zero(T)
		sDec = zero(T)
		for i in 1:nk
			x = X[i+offs]
			s    += x
			sInc += (i-1)*x
			sDec += (nk-i)*x
		end
		intervalSums[k] = s
		incSums[k] = sInc
		decSums[k] = sDec
	end
	intervalSums,incSums,decSums
end
leafsums(X::AbstractArray{T}, leafSize::Int) where T = leafsums(promote_type(T,Float64),X,leafSize)

function aggregatesums(intervalSums::Vector{T},incSums::Vector{T},decSums::Vector{T},n::Int,intervalSize::Int) where T
	# needs more parameters, total length and current interval size

	nbrRanges = length(intervalSums)
	@assert nbrRanges==length(incSums)==length(decSums)
	nbrRanges2 = div(nbrRanges-1,2)+1
	intervalSums2 = Vector{T}(undef, nbrRanges2)
	incSums2      = Vector{T}(undef, nbrRanges2)
	decSums2      = Vector{T}(undef, nbrRanges2)

	for k in 1:nbrRanges2-1
		intervalSums2[k] = intervalSums[2k-1] + intervalSums[2k]
		incSums2[k] = incSums[2k-1] + incSums[2k] + intervalSize*intervalSums[2k]
		decSums2[k] = decSums[2k-1] + decSums[2k] + intervalSize*intervalSums[2k-1]
	end
	# special case for final one or two (possibly partial) subintervals

	k = nbrRanges2
	if isodd(nbrRanges)
		intervalSums2[k] = intervalSums[2k-1]
		incSums2[k] = incSums[2k-1]
		decSums2[k] = decSums[2k-1]
	else
		intervalSums2[k] = intervalSums[2k-1] + intervalSums[2k]
		incSums2[k] = incSums[2k-1] + incSums[2k] + intervalSize*intervalSums[2k]
		decSums2[k] = decSums[2k-1] + decSums[2k] + mod1(n,intervalSize)*intervalSums[2k-1]
	end

	intervalSums2, incSums2, decSums2
end


struct SumTree{T}
	leafSize::Int
	intervalSums::Vector{Vector{T}}
	incSums::Vector{Vector{T}}
	decSums::Vector{Vector{T}}
end
function SumTree(X, leafSize)
	intervalSums,incSums,decSums = leafsums(X,leafSize)
	allIntervalSums,allIncSums,allDecSums = [intervalSums],[incSums],[decSums]

	intervalSize = leafSize
	while length(intervalSums)>1
		intervalSums,incSums,decSums = aggregatesums(intervalSums,incSums,decSums,length(X),intervalSize)
		push!(allIntervalSums,intervalSums)
		push!(allIncSums,incSums)
		push!(allDecSums,decSums)
		intervalSize *= 2
	end

	SumTree(leafSize, reverse(allIntervalSums), reverse(allIncSums), reverse(allDecSums))
end
depth(t::SumTree) = length(t.intervalSums)



#---------------------------------------------------------------------------------------------------
# New algorithm v2

struct Block2
	depth::Int
	k1::Int
	k2::Int
	lb::Float64
	ub::Float64
end
Base.isless(b1::Block2,b2::Block2) = isless(b1.ub-b1.lb, b2.ub-b2.lb)


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


function bounds2(::Val{N}, ::Type{T}, d, k1, k2, α, X, tree::SumTree)::Tuple{T,T} where {N,T}
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


function ssumstep!(heap::BinaryMaxHeap{Block2}, ::Val{N}, ::Type{T},α,X,tree::SumTree)::Tuple{T,T} where {N,T}
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

		lb11,ub11 = bounds2(Val{N}(), T, d, k1, k2, α, X, tree)
		lb += lb11
		ub += ub11
		lb11!=ub11 && push!(heap, Block2(d, k1, k2, lb11, ub11))

		if k1<k2<=nbrBlocks # only care about upper triangular part
			lb21,ub21 = bounds2(Val{N}(), T, d, k1+1, k2, α, X, tree)
			lb += lb21
			ub += ub21
			lb21!=ub21 && push!(heap, Block2(d, k1+1, k2, lb21, ub21))
		end

		if k2<nbrBlocks
			lb12,ub12 = bounds2(Val{N}(), T, d, k1, k2+1, α, X, tree)
			lb += lb12
			ub += ub12
			lb12!=ub12 && push!(heap, Block2(d, k1, k2+1, lb12, ub12))

			lb22,ub22 = bounds2(Val{N}(), T, d, k1+1, k2+1, α, X, tree)
			lb += lb22
			ub += ub22
			lb22!=ub22 && push!(heap, Block2(d, k1+1, k2+1, lb22, ub22))
		end
	end
	lb,ub
end


# Let s = ∑ᵢ∑ⱼϕⁱᵛ(α⁻¹(Xᵢ-Xⱼ)).
# returns lower and upper bounds for s, terminating when the tolerance is good enough
function ssumapprox(::Type{T},α,X,tree::SumTree;rtol)::Tuple{T,T} where T
	heap = BinaryMaxHeap{Block2}()
	lb,ub = bounds2(Val{4}(),T,1,1,1,α,X,tree)
	push!(heap,Block2(1,1,1,lb,ub))
	while !isempty(heap) && !isapprox(lb,ub;rtol=rtol)
		lb,ub = (lb,ub) .+ ssumstep!(heap,Val{4}(),T,α,X,tree)
	end
	lb,ub
end

# Let s = ∑ᵢ∑ⱼϕⁱᵛ(α⁻¹(Xᵢ-Xⱼ)).
# returns sign(s-C)
function ssign2(::Type{T},α,X,tree::SumTree,C)::Int where T
	heap = BinaryMaxHeap{Block2}()
	lb,ub = bounds2(Val{4}(),T,1,1,1,α,X,tree)
	push!(heap,Block2(1,1,1,lb,ub))
	while !isempty(heap) && lb<=C<=ub
		lb,ub = (lb,ub) .+ ssumstep!(heap,Val{4}(),T,α,X,tree)
	end
	lb > C && return 1
	ub < C && return -1
	@warn "Numerical issues for upper/lower bounds, ub=$ub, lb=$lb, C=$C."
	return 0
end


# Computes the sign of the objective function by gradually approximating using lower/upper bounds
function objectivesign2(h, X::AbstractArray{T}, tree::SumTree, α2Constant) where T
	n = length(X)
	α2 = α2Constant*h^(5/7) # 1.357[SD(a)/TD(b)]^(1/7) * h^(5/7)

	C = (n-1)*α2Constant^5*h^(-10/7)/(2*√π)
	C -= n*1.1968268412042982 # n*ϕ4(0) - get rid of the diagonal entries
	C/=2 # because of symmetry, we can sum over j>i (upper triangular part), effectively halving the sum

	ssign2(promote_type(Float64,T), α2, X, tree, C)
end

function SD_bounded2(α, X::AbstractArray{T}, tree::SumTree; rtol=0.05) where T
	# (n(n-1))⁻¹α⁻⁵∑ᵢ∑ⱼϕⁱᵛ(α⁻¹(Xᵢ-Xⱼ))
	n = length(X)
	lb,ub = ssumapprox(promote_type(Float64,T), α, X, tree; rtol=rtol)
	s = (lb+ub) + n*1.1968268412042982 # 2*(lb+ub)/2 + n*ϕ4(0) - diagonal entries
	s/(α^5*n*(n-1)) # TODO: get rid of n*(n-1) factor from this and TD
end


function _bwsj_bounded2(X; rtol=0.1)
	n = length(X)
	tree = SumTree(X,10)

	q25,q75 = quantile(X,(0.25,0.75))
	λ = min(q75-q25, std(X)*1.349) # Givivng the same scale estimate as e.g. the KernSmooth R package.
	a = 0.920λ*n^(-1/7)
	b = 0.912λ*n^(-1/9)
	# α2Constant = 1.357*(SD(a,X)/TD(b,X))^(1/7)
	α2Constant = 1.357*(SD_bounded2(a,X,tree)/TD(b,X))^(1/7)


	# Decide if we need to deal with bounds better.
	# But note that execution is very fast for bad guesses for h.

	# The Bandwidth scales roughly with n^(-1/5)
	lower = λ*n^(-1/5)*1e-9
	upper = λ*n^(-1/5)*1e9

	find_zero(h->objectivesign2(h,X,tree,α2Constant), (lower,upper), Bisection(), xrtol=rtol)
end


bwsj_bounded2(X; kwargs...) = _bwsj_bounded2(issorted(X) ? X : sort(X); kwargs...)





#---------------------------------------------------------------------------------------------------
# New algorithm v1


struct Block
	i1::Int
	i2::Int
	j1::Int
	j2::Int
	lb::Float64
	ub::Float64
end
Base.isless(b1::Block,b2::Block) = isless(b1.ub-b1.lb, b2.ub-b2.lb)


function bounds(::Type{T}, i1, i2, j1, j2, α, X)::Tuple{T,T} where T
	if (i2-i1)<=20 || (j2-j1)<=20 # fallback - loop over the data points
		# @info "Fallback for block $block"
		s = zero(T)
		if i1==j1
			@inbounds for i in i1:i2
				for j in max(i+1,j1):j2 # upper triangular part
					s += ϕ4((X[i]-X[j])/α)
				end
			end
		else
			@inbounds for i in i1:i2
				for j in j1:j2
					s += ϕ4((X[i]-X[j])/α)
				end
			end
		end
		s,s # exact bounds
	elseif i1==j1 # block on the diagonal
		@assert i2==j2
		npoints = div((i2-i1+1)*(i2-i1),2)
		# npoints.*(-0.7399861849949221, 1.1968268412042982) # TODO: improve this
		npoints.*ϕ4extrema(0.0, (X[i2]-X[i1])/α) # TODO: improve this.
	else
		npoints = (i2-i1+1)*(j2-j1+1)
		npoints.*ϕ4extrema((X[j1]-X[i2])/α, (X[j2]-X[i1])/α)
	end
end



# Let s = ∑ᵢ∑ⱼϕⁱᵛ(α⁻¹(Xᵢ-Xⱼ)).
# returns sign(s-C)
function ssign(::Type{T},α,X,C)::Int where T
	n = length(X)
	C -= n*1.1968268412042982 # n*ϕ4(0) - get rid of the diagonal entries
	C/=2 # because of symmetry, we can sum over j>i (upper triangular part), effectively halving the sum

	heap = BinaryMaxHeap{Block}()
	lb,ub = bounds(T, 1, n, 1, n, α, X)
	push!(heap,Block(1,n,1,n,lb,ub))

	while true
		lb > C && return 1
		ub < C && return -1
		isempty(heap) && break
		block = pop!(heap)

		# remove existing bounds
		lb -= block.lb
		ub -= block.ub

		# create child blocks
		midI = div(block.i1+block.i2,2)
		midJ = div(block.j1+block.j2,2)

		lb11,ub11 = bounds(T, block.i1, midI, block.j1, midJ, α, X)
		lb += lb11
		ub += ub11
		lb11!=ub11 && push!(heap, Block(block.i1, midI, block.j1, midJ, lb11, ub11))

		lb12,ub12 = bounds(T, block.i1, midI, midJ+1, block.j2, α, X)
		lb += lb12
		ub += ub12
		lb12!=ub12 && push!(heap, Block(block.i1, midI, midJ+1, block.j2, lb12, ub12))

		if block.j1>block.i1 # only care about upper triangular part
			lb21,ub21 = bounds(T, midI+1, block.i2, block.j1, midJ, α, X)
			lb += lb21
			ub += ub21
			lb21!=ub21 && push!(heap, Block(midI+1, block.i2, block.j1, midJ, lb21, ub21))
		end

		lb22,ub22 = bounds(T, midI+1, block.i2, midJ+1, block.j2, α, X)
		lb += lb22
		ub += ub22
		lb22!=ub22 && push!(heap, Block(midI+1, block.i2, midJ+1, block.j2, lb22, ub22))
	end

	@warn "Numerical issues for upper/lower bounds, ub=$ub, lb=$lb, C=$C."
	return 0
end




# Computes the sign of the objective function by gradually approximating using lower/upper bounds
function objectivesign(h, X::AbstractArray{T}, α2Constant) where T
	n = length(X)
	α2 = α2Constant*h^(5/7) # 1.357[SD(a)/TD(b)]^(1/7) * h^(5/7)

	C = (n-1)*α2Constant^5*h^(-10/7)/(2*√π)
	ssign(promote_type(Float64,T), α2, X, C)
end


function _bwsj_bounded(X; rtol=0.1)
	n = length(X)
	q25,q75 = quantile(X,(0.25,0.75))
	λ = min(q75-q25, std(X)*1.349) # Givivng the same scale estimate as e.g. the KernSmooth R package.
	a = 0.920λ*n^(-1/7)
	b = 0.912λ*n^(-1/9)
	α2Constant = 1.357*(SD(a,X)/TD(b,X))^(1/7)

	# Decide if we need to deal with bounds better.
	# But note that execution is very fast for bad guesses for h.
	lower = λ*n^(-1/5)*1e-9
	upper = λ*n^(-1/5)*1e9

	find_zero(h->objectivesign(h,X,α2Constant), (lower,upper), Bisection(), xrtol=rtol)
end


bwsj_bounded(X; kwargs...) = _bwsj_bounded(issorted(X) ? X : sort(X); kwargs...)



#---------------------------------------------------------------------------------------------------



function SDrec(::Type{T},α,X,i1,i2,j1,j2)::T where T
	j1<i1 && return zero(T) # early out, skip lower triangular part


	if (i2-i1)<=5 || (j2-j1)<=5 # fallback - loop over the data points
		s = zero(T)

		if i1==j1
			@inbounds for i in i1:i2
				for j in max(i+1,j1):j2 # upper triangular part
					s += ϕ4((X[i]-X[j])/α)
				end
			end
		else
			@inbounds for i in i1:i2
				for j in j1:j2
					s += ϕ4((X[i]-X[j])/α)
				end
			end
		end
		return s
	end

	if i1!=j1 # only consider approximation if it's not the same interval
		npoints = (i2-i1+1)*(j2-j1+1)
		lb,ub = npoints.*ϕ4extrema((X[j1]-X[i2])/α, (X[j2]-X[i1])/α) # proper bounds

		if ub-lb < 1 # 20 # 0.1 # 0.01 # what is a good bound?
			return (lb+ub)/2 # approximation
		end
	end

	# split and recurse
	midI = div(i1+i2,2)
	midJ = div(j1+j2,2)
	SDrec(T,α,X,i1,midI,j1,midJ) + SDrec(T,α,X,i1,midI,midJ+1,j2) + SDrec(T,α,X,midI+1,i2,j1,midJ) + SDrec(T,α,X,midI+1,i2,midJ+1,j2)
end


# assumes X is sorted
function SDapprox(α, X::AbstractArray{T}) where T
	# @show α
	n = length(X)
	s = 2*SDrec(promote_type(T,Float64),α,X,1,n,1,n)
	s += n*ϕ4(0) # diagonal

	# @info s
	@assert s>0 "s>0 failed: $s"

	s/(n*(n-1)*α^5)
end



# unoptimized version
function SD(α, X::AbstractArray{T}) where T
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

function TD(b, X::AbstractArray{T}) where T
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


function objective(h, X, α2Constant)
	# [ R(K) / (nσ⁴ₖSD(α₂(h))) ]^(1/5) - h = 0
	# R(K) = 1/(2√π)
	# σ⁴ₖ = 1
	n = length(X)
	α2 = α2Constant*h^(5/7) # 1.357[SD(a)/TD(b)]^(1/7) * h^(5/7)
	# denom = (n*2*√π) * SD(α2,X)
	# denom^(-1/5) - h
	1/(2*√(π)*n*h^5) - SD(α2,X)
end


function _bwsj(X; rtol=0.1)
	n = length(X)
	q25,q75 = quantile(X,(0.25,0.75))
	λ = min(q75-q25, std(X)*1.349) # Givivng the same scale estimate as e.g. the KernSmooth R package.
	a = 0.920λ*n^(-1/7)
	b = 0.912λ*n^(-1/9)
	α2Constant = 1.357*(SD(a,X)/TD(b,X))^(1/7)


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

	find_zero(h->objective(h,X,α2Constant), (lower,upper), Bisection(), xrtol=rtol)
end


bwsj(X; kwargs...) = _bwsj(issorted(X) ? X : sort(X); kwargs...)
# bwsj(X; kwargs...) = _bwsj_bounded2(issorted(X) ? X : sort(X); kwargs...)

