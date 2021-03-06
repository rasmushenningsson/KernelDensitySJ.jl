using Roots

# Unoptimized reference implementation
# NB: ϕ is scaled by √2π compared to standard definition
function SD_reference(α, X::AbstractArray{T}) where T
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

# NB: ϕ is scaled by √2π compared to standard definition
function TD_reference(b, X::AbstractArray{T}) where T
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


function objective_reference(h, X, α2Constant)
	# [ R(K) / (nσ⁴ₖSD(α₂(h))) ]^(1/5) - h = 0
	# R(K) = 1/(2√π)
	# σ⁴ₖ = 1
	n = length(X)
	α2 = α2Constant*h^(5/7) # 1.357[SD(a)/TD(b)]^(1/7) * h^(5/7)
	1/(√2*n*h^5) - SD_reference(α2,X)
end


function _bwsj_reference(X; rtol=0.1)
	n = length(X)
	q25,q75 = quantile(X,(0.25,0.75))
	λ = min(q75-q25, std(X)*1.349) # Givivng the same scale estimate as e.g. the KernSmooth R package.
	a = 0.920λ*n^(-1/7)
	b = 0.912λ*n^(-1/9)
	α2Constant = 1.357*(SD_reference(a,X)/TD_reference(b,X))^(1/7)

	# Figure out lower and upper bounds for search.
	# The Bandwidth scales roughly with n^(-1/5). TODO: Use this fact!

	lower = 0.1*λ
	for i=1:11
		i==11 && error("Failed to find lower bound")
		objective_reference(lower,X,α2Constant)>0 && break
		lower/=10
	end
	upper = 1*λ
	for i=1:11
		i==11 && error("Failed to find upper bound")
		objective_reference(upper,X,α2Constant)<0 && break
		upper*=10
	end

	find_zero(h->objective_reference(h,X,α2Constant), (lower,upper), Bisection(), xrtol=rtol)
end


bwsj_reference(X; kwargs...) = _bwsj_reference(issorted(X) ? X : sort(X); kwargs...)


function _smooth_reference(x,y,c,xeval,w)
	w .= -c.*(x.-xeval).^2
	w .= exp.(w .- maximum(w)) # multiply all weights with a factor for improved precision
	y'w / sum(w)
end

function smooth_reference(x::AbstractVector{<:T}, y::AbstractVector{<:S}, bandwidth, xeval) where {T<:Real,S<:Real}
	@assert length(x)==length(y)
	c = 1 ./(2 .* bandwidth.^2)
	w = similar(x, promote_type(T,eltype(bandwidth),eltype(xeval)))
	_smooth_reference.(Ref(x),Ref(y),c,xeval,Ref(w))
end

function _density_reference(x,bandwidth,c,xeval,w)
	w .= exp.(-c.*(x.-xeval).^2)
	sum(w) / (sqrt(2π)*bandwidth*length(x))
end

function density_reference(x::AbstractVector{<:T}, bandwidth, xeval) where {T<:Real}
	c = 1 ./(2 .* bandwidth.^2)
	w = similar(x, promote_type(T,eltype(bandwidth),eltype(xeval)))
	_density_reference.(Ref(x),bandwidth,c,xeval,Ref(w))
end
