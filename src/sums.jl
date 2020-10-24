function binary_reduce(op::F, ::Type{T}, x, final) where {F,T}
	N = length(x)
	out = Vector{T}(undef, div(N+1,2))
	for i=1:div(N,2)
		out[i] = op(x[2i-1], x[2i])
	end
	isodd(N) && (out[end] = final)
	out
end
binary_reduce(op, x::AbstractVector{T}, final) where T = binary_reduce(op, T, x, final)

diagonalsum(x) = (N=length(x); sum(z->z[1]*z[2], zip(x,range(-N+1;length=N,step=2))))

# Function barrier to help compiler
diagonal_reduction(::Type{T},N,diagonalSums,intervalSums,intervalSize) where T =
	binary_reduce((i1,i2)->diagonalSums[i1]+diagonalSums[i2]+intervalSize*intervalSums[i2] - min(N-i1*intervalSize,intervalSize)*intervalSums[i1], T, 1:length(diagonalSums), diagonalSums[end])

struct SumTree{T}
	leafSize::Int
	intervalSums::Vector{Vector{T}}
	diagonalSums::Vector{Vector{T}}
end
function SumTree(X::AbstractVector{T}, leafSize) where T
	N = length(X)

	intervalSums = map(sum, Iterators.partition(X,leafSize))
	diagonalSums = map(diagonalsum, Iterators.partition(X,leafSize))
	allIntervalSums,allDiagonalSums = [intervalSums],[diagonalSums]

	intervalSize = leafSize
	while length(intervalSums)>1
		diagonalSums = diagonal_reduction(T, N, diagonalSums, intervalSums, intervalSize)
		intervalSums = binary_reduce(+,intervalSums,intervalSums[end])

		push!(allIntervalSums,intervalSums)
		push!(allDiagonalSums,diagonalSums)
		intervalSize *= 2
	end

	SumTree(leafSize, reverse(allIntervalSums), reverse(allDiagonalSums))
end
depth(t::SumTree) = length(t.intervalSums)
