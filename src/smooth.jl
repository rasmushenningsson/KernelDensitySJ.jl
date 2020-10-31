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


function smoothapprox(::Type{T},C,D,X,Y,tree,xeval;rtol) where T
	heap = BinaryMaxHeap{SmoothBlock}()

	lb,ub = smoothbounds(T,1,1,C,D,X,tree,xeval)

	push!(heap,SmoothBlock(1,1,lb,ub))
	while !isempty(heap) && !isapprox(lb,ub;rtol=rtol)
		lb,ub = (lb,ub) .+ smoothstep!(T,heap,C,D,X,Y,tree,xeval)
	end
	lb,ub
end



function _gaussiansmoothing(weightTree::WeightTree{T}, minMaxTree::MinMaxTree{T}, x, y, C, xeval, rtol) where T
	# compute constant D such that the weight of the closest point is normalized to 1, which is critical for numerical precision.
	r = searchsorted(x,xeval)
	D1 = (x[min(first(r),length(x))]-xeval)
	D2 = (x[max(last(r),1)]-xeval)
	D = C * (abs(D1)<abs(D2) ? D1 : D2)

	# TODO: tighten rtol for numerator & denominator to get the desired accuracy in the final result
	denom = smoothapprox(T, C, D, x, nothing, weightTree, xeval; rtol=rtol)
	numer = smoothapprox(T, C, D, x, y,       minMaxTree, xeval; rtol=rtol)

	# use that the denominator is positive
	lb = min((numer[1]./denom)...)
	ub = max((numer[2]./denom)...)
	lb,ub
	(lb+ub)/2
end


function gaussiansmoothing(x::AbstractVector{T}, y::AbstractVector{T}, bandwidth, xeval; leafSize=10, rtol=1e-3) where T
	@assert length(x)==length(y)
	c = 1/(2*bandwidth^2)

	if !issorted(x)
		perm = sortperm(x)
		x = x[perm]
		y = y[perm]
	end

	# ϕ(x) = 1/√(2π) * exp(-x^2/2)
	# c = 1/(2*bandwidth^2)
	# w .= exp.(-c*(x.-xeval[i]).^2)
	# exp(-c*(x-xeval).^2)
	# exp(-((x-xeval)/bandwidth).^2 / 2)
	C = 1/bandwidth

	# weightTree = WeightTree(x,leafSize)
	minMaxTree = MinMaxTree(x,y,leafSize)
	weightTree = WeightTree(minMaxTree)

	_gaussiansmoothing.(Ref(weightTree),Ref(minMaxTree),Ref(x),Ref(y),C,xeval,rtol)
end
