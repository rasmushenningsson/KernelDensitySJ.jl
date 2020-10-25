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

diagonalsum(x) = (N=length(x); sum(z->z[1]*z[2], zip(x,-N+1:2:N-1)))


# TESTING
mutable struct DiagonalContext{T}
	allIntervalSums::T
	intervalSize::Int
	N::Int
	level::Int
end
DiagonalContext(allIntervalSums, leafSize, N) =
	DiagonalContext(allIntervalSums, leafSize, N, 1)

# Function barrier to help compiler
diagonal_reduction(::Type{T},N,diagonalSums,intervalSums,intervalSize) where T =
	binary_reduce((i1,i2)->diagonalSums[i1]+diagonalSums[i2]+intervalSize*intervalSums[i2] - min(N-i1*intervalSize,intervalSize)*intervalSums[i1], T, 1:length(diagonalSums), diagonalSums[end])
function diagonal_reduction(diagonalSums::AbstractVector{T}, context::DiagonalContext) where T
	diagonalSums = diagonal_reduction(T, context.N, diagonalSums, context.allIntervalSums[context.level], context.intervalSize)
	context.level += 1
	context.intervalSize *= 2
	diagonalSums
end

function pyramid(f, r)
	v = [r]
	while length(r)>1
		r = f(r)
		push!(v, r)
	end
	v
end


struct SumTree{T}
	leafSize::Int
	intervalSums::Vector{Vector{T}}
	diagonalSums::Vector{Vector{T}}
end
function SumTree(X::AbstractVector{T}, leafSize) where T
	N = length(X)

	# TESTING
	intervalSums = map(sum, Iterators.partition(X,leafSize))
	allIntervalSums = pyramid(z->binary_reduce(+,z,z[end]), intervalSums)

	context = DiagonalContext(allIntervalSums, leafSize, length(X))
	diagonalSums = map(diagonalsum, Iterators.partition(X,leafSize))
	allDiagonalSums = pyramid(z->diagonal_reduction(z,context), diagonalSums)


	SumTree(leafSize, reverse(allIntervalSums), reverse(allDiagonalSums))
end
depth(t::SumTree) = length(t.intervalSums)
