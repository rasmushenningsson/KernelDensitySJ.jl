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

