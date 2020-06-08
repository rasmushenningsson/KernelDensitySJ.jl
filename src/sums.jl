function leafsums(::Type{T}, X, leafSize::Int) where T
	n = length(X)
	nbrRanges = div(n-1,leafSize)+1
	intervalSums = Vector{T}(undef, nbrRanges)
	diagonalSums = Vector{T}(undef, nbrRanges)

	for k in 1:nbrRanges
		offs = (k-1)*leafSize
		nk = k<nbrRanges ? leafSize : mod1(n,leafSize)

		s = zero(T)
		sDiag = zero(T)
		for i in 1:nk
			x = X[i+offs]
			s     += x
			sDiag += (2i-nk-1)*x
		end
		intervalSums[k] = s
		diagonalSums[k] = sDiag
	end
	intervalSums,diagonalSums
end
leafsums(X::AbstractArray{T}, leafSize::Int) where T = leafsums(promote_type(T,Float64),X,leafSize)

function aggregatesums(intervalSums::Vector{T},diagonalSums::Vector{T},n::Int,intervalSize::Int) where T
	nbrRanges = length(intervalSums)
	@assert nbrRanges==length(diagonalSums)
	nbrRanges2 = div(nbrRanges-1,2)+1
	intervalSums2 = Vector{T}(undef, nbrRanges2)
	diagonalSums2 = Vector{T}(undef, nbrRanges2)

	for k in 1:nbrRanges2-1
		intervalSums2[k] = intervalSums[2k-1] + intervalSums[2k]
		diagonalSums2[k] = diagonalSums[2k-1] + diagonalSums[2k] + intervalSize*(intervalSums[2k]-intervalSums[2k-1])
	end
	# special case for final one or two (possibly partial) subintervals
	k = nbrRanges2
	if isodd(nbrRanges)
		intervalSums2[k] = intervalSums[2k-1]
		diagonalSums2[k] = diagonalSums[2k-1]
	else
		intervalSums2[k] = intervalSums[2k-1] + intervalSums[2k]
		diagonalSums2[k] = diagonalSums[2k-1] + diagonalSums[2k] + intervalSize*intervalSums[2k] - mod1(n,intervalSize)*intervalSums[2k-1]
	end

	intervalSums2, diagonalSums2
end


struct SumTree{T}
	leafSize::Int
	intervalSums::Vector{Vector{T}}
	diagonalSums::Vector{Vector{T}}
end
function SumTree(X, leafSize)
	intervalSums,diagonalSums = leafsums(X,leafSize)
	allIntervalSums,allDiagonalSums = [intervalSums],[diagonalSums]

	intervalSize = leafSize
	while length(intervalSums)>1
		intervalSums,diagonalSums = aggregatesums(intervalSums,diagonalSums,length(X),intervalSize)
		push!(allIntervalSums,intervalSums)
		push!(allDiagonalSums,diagonalSums)
		intervalSize *= 2
	end

	SumTree(leafSize, reverse(allIntervalSums), reverse(allDiagonalSums))
end
depth(t::SumTree) = length(t.intervalSums)

