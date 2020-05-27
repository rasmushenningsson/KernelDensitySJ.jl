# initial tests for implementation of the "bw.sj" method, from 
# Sheather and Jones, "A Reliable Data-based Bandwidth Selection Method for Kernel Density Estimation"
# Notation below follows roughly the notation in the paper.



ϕ(x) = 1/√(2π) * exp(-x^2/2)

# ϕ4(x) = (x^4 - 6x^2 + 3) * ϕ(x) # slow
# ϕ6(x) = (x^6 - 15x^4 + 45x^2 - 15) * ϕ(x) # slow

ϕ4(x) = (x2=x*x; (x2*x2 - 6x2 + 3) * ϕ(x))
ϕ6(x) = (x2=x*x; x4=x2*x2; (x4*x2 - 15x4 + 45x2 - 15) * ϕ(x))


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



#---------------------------------------------------------------------------------------------------
# New algorithm test


struct Block
	i1::Int
	i2::Int
	j1::Int
	j2::Int
end

struct Bounds
	lb::Float64
	ub::Float64
end
function Base.isless(b1::Bounds,b2::Bounds)
	r1 = b1.ub-b1.lb
	r2 = b2.ub-b2.lb
	r1!=r2 ? isless(r2,r1) : isless(b2.lb,b1.lb) # NB: this order to sort by max instead of mean
end
isexact(b::Bounds) = b.lb===b.ub


function bounds(::Type{T}, block::Block, α, X) where T
	i1,i2,j1,j2 = block.i1,block.i2,block.j1,block.j2

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
		Bounds(s,s) # exact bounds
	elseif i1==j1 # block on the diagonal
		@assert i2==j2
		npoints = div((i2-i1+1)*(i2-i1),2)
		# lb,ub = npoints.*(-0.7399861849949221, 1.1968268412042982) # TODO: improve this
		lb,ub = npoints.*ϕ4extrema(0.0, (X[i2]-X[i1])/α) # TODO: improve this?
		Bounds(lb,ub)
	else
		npoints = (i2-i1+1)*(j2-j1+1)
		lb,ub = npoints.*ϕ4extrema((X[j1]-X[i2])/α, (X[j2]-X[i1])/α)
		Bounds(lb,ub)
	end
end



# Let s = ∑ᵢ∑ⱼϕⁱᵛ(α⁻¹(Xᵢ-Xⱼ)).
# returns sign(s-C)
function ssign(::Type{T},α,X,C)::Int where T
	n = length(X)
	C -= n*1.1968268412042982 # n*ϕ4(0) - get rid of the diagonal entries
	C/=2 # because of symmetry, we can sum over j>i (upper triangular part), effectively halving the sum

	root = Block(1,n,1,n)
	queue = PriorityQueue{Block,Bounds}()
	rootBounds = bounds(T, root, α, X)
	enqueue!(queue,root,rootBounds)

	lb,ub = rootBounds.lb,rootBounds.ub

	while true
		lb > C && return 1
		ub < C && return -1
		isempty(queue) && break
		block,blockBounds = dequeue_pair!(queue)

		# remove existing bounds
		lb -= blockBounds.lb
		ub -= blockBounds.ub

		# create child blocks
		midI = div(block.i1+block.i2,2)
		midJ = div(block.j1+block.j2,2)

		child11 = Block(block.i1, midI, block.j1, midJ)
		childBounds11 = bounds(T, child11, α, X)
		lb += childBounds11.lb
		ub += childBounds11.ub
		isexact(childBounds11) || enqueue!(queue, child11, childBounds11)

		child12 = Block(block.i1, midI, midJ+1, block.j2)
		childBounds12 = bounds(T, child12, α, X)
		lb += childBounds12.lb
		ub += childBounds12.ub
		isexact(childBounds12) || enqueue!(queue, child12, childBounds12)

		if block.j1>block.i1 # only care about upper triangular part
			child21 = Block(midI+1, block.i2, block.j1, midJ)
			childBounds21 = bounds(T, child21, α, X)
			lb += childBounds21.lb
			ub += childBounds21.ub
			isexact(childBounds21) || enqueue!(queue, child21, childBounds21)
		end

		child22 = Block(midI+1, block.i2, midJ+1, block.j2)
		childBounds22 = bounds(T, child22, α, X)
		lb += childBounds22.lb
		ub += childBounds22.ub
		isexact(childBounds22) || enqueue!(queue, child22, childBounds22)
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
	lower = λ*1e-9
	upper = λ*1e9

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
function SD(α, X::AbstractArray{T}) where T
	# @show α
	n = length(X)
	s = 2*SDrec(promote_type(T,Float64),α,X,1,n,1,n)
	s += n*ϕ4(0) # diagonal

	# @info s
	@assert s>0 "s>0 failed: $s"

	s/(n*(n-1)*α^5)
end



# unoptimized version
function SDreference(α, X::AbstractArray{T}) where T
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

