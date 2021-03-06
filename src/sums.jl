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

function pyramid(f, x, nbrLevels)
	v = [x]
	for i=2:nbrLevels
		push!(v, f(v[end]))
	end
	v
end

diagonalsum(x) = (N=length(x); sum(z->z[1]*z[2], zip(x,-N+1:2:N-1)))

# Function barrier to help compiler
diagonal_reduction(::Type{T},N,diagonalSums,intervalSums,intervalSize) where T =
	binary_reduce((i1,i2)->diagonalSums[i1]+diagonalSums[i2]+intervalSize*intervalSums[i2] - min(N-i1*intervalSize,intervalSize)*intervalSums[i1], T, 1:length(diagonalSums), diagonalSums[end])
function diagonal_reduction(t::Tuple{AbstractVector{T},Int}, allIntervalSums, N, leafSize) where T
	diagonalSums,level = t[1],t[2]
	diagonalSums = diagonal_reduction(T, N, diagonalSums, allIntervalSums[level], leafSize*(2^(level-1)))
	diagonalSums,level+1
end

nbrlevels(N,leafSize) = Int(log2(nextpow(2, div(N+leafSize-1,leafSize))))+1

abstract type Tree end
depth(t::Tree) = length(t.intervalSums)


struct SumTree{T} <: Tree
	leafSize::Int
	intervalSums::Vector{Vector{T}}
	diagonalSums::Vector{Vector{T}}
end
function SumTree(X::AbstractVector{T}, leafSize) where T
	N = length(X)
	nbrLevels = nbrlevels(N,leafSize)

	intervalSums = map(sum, Iterators.partition(X,leafSize))
	allIntervalSums = pyramid(z->binary_reduce(+,z,z[end]), intervalSums, nbrLevels)

	diagonalSums = map(diagonalsum, Iterators.partition(X,leafSize))
	allDiagonalSums = first.(pyramid(z->diagonal_reduction(z,allIntervalSums,N,leafSize), (diagonalSums,1), nbrLevels))

	SumTree(leafSize, reverse!(allIntervalSums), reverse!(allDiagonalSums))
end

struct MinMaxTree{T} <: Tree
	leafSize::Int
	intervalSums::Vector{Vector{T}}
	mins::Vector{Vector{T}}
	maxes::Vector{Vector{T}}
end
function MinMaxTree(X, Y, leafSize)
	@assert length(X)==length(Y)
	nbrLevels = nbrlevels(length(X),leafSize)
	intervalSums = map(sum, Iterators.partition(X,leafSize))
	allIntervalSums = pyramid(z->binary_reduce(+,z,z[end]), intervalSums, nbrLevels)
	mins = map(minimum, Iterators.partition(Y,leafSize))
	allMins = pyramid(z->binary_reduce(min,z,z[end]), mins, nbrLevels)
	maxes = map(maximum, Iterators.partition(Y,leafSize))
	allMaxes = pyramid(z->binary_reduce(max,z,z[end]), maxes, nbrLevels)
	MinMaxTree(leafSize,reverse!(allIntervalSums),reverse!(allMins),reverse!(allMaxes))
end

struct WeightTree{T} <: Tree
	leafSize::Int
	intervalSums::Vector{Vector{T}}
end
function WeightTree(X, leafSize)
	nbrLevels = nbrlevels(length(X),leafSize)
	intervalSums = map(sum, Iterators.partition(X,leafSize))
	allIntervalSums = pyramid(z->binary_reduce(+,z,z[end]), intervalSums, nbrLevels)
	WeightTree(leafSize,reverse!(allIntervalSums))
end
WeightTree(tree::MinMaxTree) = WeightTree(tree.leafSize, tree.intervalSums)
