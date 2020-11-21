struct SmoothBlock
	depth::Int
	i::Int
	lb::Float64
	ub::Float64
end
Base.isless(b1::SmoothBlock,b2::SmoothBlock) = isless(b1.ub-b1.lb, b2.ub-b2.lb)


adjustbounds(lb,ub,::WeightTree,::Any,::Any) = lb,ub
function adjustbounds(lb,ub,tree::MinMaxTree,depth,i)
	m = tree.mins[depth][i]
	M = tree.maxes[depth][i]
	min(m*lb,m*ub), max(M*lb,M*ub) # assumes that 0≤lb≤ub
end


function smoothexact(::Type{T},i,C,D,X,Y,xeval,leafSize) where T
	n = length(X)
	k1 = (i-1)*leafSize+1
	k2 = min(i*leafSize, n)

	# TODO: optimize (avoid allocations etc)
	# TODO: use dispatch instead
	if isnothing(Y)
		sum(ϕ.(C.*(X[k1:k2].-xeval),D))
	else
		sum(Y[k1:k2].*ϕ.(C.*(X[k1:k2].-xeval),D))
	end
end


function smoothbounds(::Type{T}, depth, i, C, D, X, tree, xeval) where T
	n = length(X)
	intervalSize = tree.leafSize * 2^(length(tree.intervalSums)-depth)

	k1 = (i-1)*intervalSize+1
	k2 = min(i*intervalSize, n)
	npoints = k2-k1+1

	a = C*(X[k1]-xeval)
	b = C*(X[k2]-xeval)

	xMean = C*(tree.intervalSums[depth][i]/npoints-xeval)

	if a<0 && b>0
		# special case, interval covers xeval
		# TODO: improve bounds? (We can use Jensen's inequality for this case too, if -1≤a≤b≤1, but ϕbounds() implementation needs to be updated)
		lb = ϕ(max(-a,b),D)
		ub = 1/√(2π)
		# return npoints.*(lb,ub)
		return adjustbounds((npoints.*(lb,ub))..., tree, depth, i)
	elseif a<0
		a,b,xMean = -b,-a,-xMean # mirror
	end

	# npoints.*ϕbounds(a,b,xMean)
	adjustbounds((npoints.*ϕbounds(a,b,xMean,D))..., tree, depth, i)
end

function smoothstep!(::Type{T},heap,C,D,X,Y,tree,xeval) where T
	block = pop!(heap)

	# remove existing bounds
	lb = -block.lb
	ub = -block.ub

	if block.depth >= depth(tree) # Fallback to exact sum
		# s = 0.0 # TODO: compute exact sum
		s = smoothexact(T,block.i,C,D,X,Y,xeval,tree.leafSize)
		lb += s
		ub += s
	else # create child blocks
		d = block.depth+1
		i = 2block.i-1 # first child
		nbrBlocks = length(tree.intervalSums[d])

		lb1,ub1 = smoothbounds(T,d,i,C,D,X,tree,xeval)
		lb += lb1
		ub += ub1
		lb1!=ub1 && push!(heap, SmoothBlock(d, i, lb1, ub1))

		if i<nbrBlocks
			lb2,ub2 = smoothbounds(T,d,i+1,C,D,X,tree,xeval)
			lb += lb2
			ub += ub2
			lb2!=ub2 && push!(heap, SmoothBlock(d, i+1, lb2, ub2))
		end
	end
	lb,ub
end


"""
	GaussianKernelSmoother(x,y; leafSize=10)

Create a callable `GaussianKernelSmoother` object of a set of data points with coordinates `x` and values `y`.

See also `smooth`.
"""
struct GaussianKernelSmoother{T,TX,TY}
	x::TX
	y::TY
	minMaxTree::MinMaxTree{T}
	weightTree::WeightTree{T}
end
function GaussianKernelSmoother(x, y; leafSize=10)
	@assert length(x)==length(y)
	if !issorted(x)
		perm = sortperm(x)
		x = x[perm]
		y = y[perm]
	end
	minMaxTree = MinMaxTree(x,y,leafSize)
	weightTree = WeightTree(minMaxTree)
	GaussianKernelSmoother(x,y,minMaxTree,weightTree)
end



function smoothapprox(gks::GaussianKernelSmoother{T},C,D,xeval;atol,rtol) where T
	heapn = BinaryMaxHeap{SmoothBlock}()
	heapd = BinaryMaxHeap{SmoothBlock}()

	lbn,ubn = smoothbounds(T,1,1,C,D,gks.x,gks.minMaxTree,xeval)
	push!(heapn,SmoothBlock(1,1,lbn,ubn))
	lbd,ubd = smoothbounds(T,1,1,C,D,gks.x,gks.weightTree,xeval)
	push!(heapd,SmoothBlock(1,1,lbd,ubd))

	while true
		lb = min(lbn/lbd,lbn/ubd)
		ub = max(ubn/lbd,ubn/ubd)

		isempty(heapn) && isempty(heapd) && return lb,ub
		isapprox(lb,ub;atol=atol,rtol=rtol) && return lb,ub

		if !isempty(heapn)
			lbn,ubn = (lbn,ubn) .+ smoothstep!(T,heapn,C,D,gks.x,gks.y,  gks.minMaxTree,xeval)
		end
		if !isempty(heapd)
			lbd,ubd = (lbd,ubd) .+ smoothstep!(T,heapd,C,D,gks.x,nothing,gks.weightTree,xeval)
		end
	end
end


"""
	(::GaussianKernelSmoother)(bandwidth, xeval; rtol=atol>0 ? 0 : 1e-3, atol=0)

Calling a GaussianKernelSmoother object evaluates the smoothed function for the given `bandwidth` and `xeval`.
The value of the smoothed function `f` at `x₀` is given by `f(x₀) := ∑ᵢyᵢwᵢ / ∑ᵢwᵢ`, where `wᵢ := exp(-(x₀-xᵢ)²/2bandwidth²)`.

The accuracy of the result is controlled by `rtol` and `atol`. Lower and upper bounds `lb≤f(x₀)≤ub` are gradually improved until `isapprox(lb,ub;rtol=rtol,atol=atol)`.

# Examples
```julia-repl
julia> g = GaussianKernelSmoother([0.0,1.0],[2.0,4.0]);
julia> g(1, 0.9)
3.197375320224904
```
"""
function (gks::GaussianKernelSmoother)(bandwidth::Real, xeval::Real; atol=0, rtol=atol>0 ? 0 : 1e-3)
	x,y = gks.x,gks.y
	C = 1 / bandwidth
	# compute constant D such that the weight of the closest point is normalized to 1, which is critical for numerical precision.
	r = searchsorted(x,xeval)
	D1 = (x[min(first(r),length(x))]-xeval)
	D2 = (x[max(last(r),1)]-xeval)
	D = C * (abs(D1)<abs(D2) ? D1 : D2)

	lb,ub = smoothapprox(gks,C,D,xeval; atol=atol, rtol=rtol)
	(lb+ub)/2
end

"""
	smooth(x, y, bandwidth, xeval; leafSize=10, rtol=atol>0 ? 0 : 1e-3, atol=0)

Evaluate the Gaussian Kernel Smoother of a set of data points with coordinates `x` and values `y`, with the specified `bandwidth`, at the coordinates in `xeval`.
The value of the smoothed function `f` at `x₀` is given by `f(x₀) := ∑ᵢyᵢwᵢ / ∑ᵢwᵢ`, where `wᵢ := exp(-(x₀-xᵢ)²/2bandwidth²)`.

The accuracy of the result is controlled by `rtol` and `atol`. Lower and upper bounds `lb≤f(x₀)≤ub` are gradually improved until `isapprox(lb,ub;rtol=rtol,atol=atol)`.

It is much more efficient to call `smooth` once with vector/matrix arguments for `xeval` and/or `bandwidth` than to call `smooth` multiple times.

# Examples
```julia-repl
julia> smooth([0.0,1.0], [2.0,4.0], 1.0, [0.1, 0.5, 0.9])
3-element Array{Float64,1}:
 2.802624679775096
 3.0
 3.197375320224904
```
See also `GaussianKernelSmoother`.
"""
function smooth(x::AbstractVector{T}, y::AbstractVector{T}, bandwidth, xeval; leafSize=10, kwargs...) where T
	gks = GaussianKernelSmoother(x,y; leafSize=leafSize)
	gks.(bandwidth, xeval; kwargs...)
end
