# initial tests for implementation of the "bw.sj" method, from 
# Sheather and Jones, "A Reliable Data-based Bandwidth Selection Method for Kernel Density Estimation"



ϕ(x) = 1/√(2π) * exp(-x^2/2)

ϕ4(x) = (x^4 - 6x^2 + 3) * ϕ(x)
ϕ6(x) = (x^6 - 15x^4 + 45x^2 - 15) * ϕ(x)


# unoptimized version
function S(α, X)
	# (n(n-1))⁻¹α⁻⁵∑ᵢ∑ⱼϕⁱᵛ(α⁻¹(Xᵢ-Xⱼ))
	n = length(X)
	C = 1/(n*(n-1)*α^5)
	Y = X./α
	C * sum(ϕ4.(Y.-Y')) # TODO: fix inefficient summation. No need to allocate N^2 matrix.
end

function T(b, X)
	# -(n(n-1))⁻¹b⁻⁷∑ᵢ∑ⱼϕᵛⁱ(b⁻¹(Xᵢ-Xⱼ))
	n = length(X)
	C = -1/(n*(n-1)*b^7)
	Y = X./b
	C * sum(ϕ6.(Y.-Y')) # TODO: fix inefficient summation. No need to allocate N^2 matrix.
end


function objective(h, X, α2Constant)
	# [ R(K) / (nσ⁴ₖS(α₂(h))) ]^(1/5) - h = 0
	# R(K) = 1/(2√π)
	# σ⁴ₖ = 1
	n = length(X)
	α2 = α2Constant*h^(5/7) # 1.357[S(a)/T(b)]^(1/7) * h^(5/7)
	denom = (n*2*√π) * S(α2,X)
	denom^(-1/5) - h
end


function bwsj(X)
	n = length(X)
	q25,q75 = quantile(X,(0.25,0.75))
	λ = min(q75-q25, std(X)*1.349) # Givivng the same scale estimate as e.g. the KernSmooth R package.
	a = 0.920λ*n^(-1/7)
	b = 0.912λ*n^(-1/9)
	α2Constant = 1.357*(S(a,X)/T(b,X))^(1/7)


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

	find_zero(h->objective(h,X,α2Constant), (lower,upper), Bisection())
end


