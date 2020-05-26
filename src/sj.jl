# initial tests for implementation of the "bw.sj" method, from 
# Sheather and Jones, "A Reliable Data-based Bandwidth Selection Method for Kernel Density Estimation"
# Notation below follows roughly the notation in the paper.



ϕ(x) = 1/√(2π) * exp(-x^2/2)

ϕ4(x) = (x^4 - 6x^2 + 3) * ϕ(x)
ϕ6(x) = (x^6 - 15x^4 + 45x^2 - 15) * ϕ(x)


function ϕ4extrema(a,b)
	# (x-√(5+√10))(x+√(5+√10))(x-√(5-√10))(x+√(5-√10))x
	# -√(5+√BigFloat(10)), √(5-√BigFloat(10)), 0, √(5-√BigFloat(10)), √(5+√BigFloat(10))
	# xExt = (-2.8569700138728056, 1.355626179974266, 0.0, 1.355626179974266, 2.8569700138728056)
	# yExt = (0.13912184973492348, -0.7399861849949221, 1.1968268412042982, -0.7399861849949221, 0.13912184973492348)
	# we only need to consider x>0 since the function is even
	xExt = (0.0, 1.355626179974266, 2.8569700138728056)
	yExt = (1.1968268412042982, -0.7399861849949221, 0.13912184973492348)

	# TODO: clean up code

	if a < 0 < b # covers 0
		M = yExt[1] # global maximum
		c = max(-a,b) # move to x>0
		m = c>xExt[2] ? yExt[2] : ϕ4(c) # covers global minimum or not?
		return m,M
	end

	if a<0
		a,b = -b,-a # move to x>0
	end

	a>xExt[3] && return ϕ4(b),ϕ4(a)

	if a>xExt[2]
		yb = ϕ4(b)
		M = b>xExt[3] ? yExt[3] : yb
		m = min(ϕ4(a),yb)
		return m,M
	end

	ya = ϕ4(a)
	b>xExt[3] && return yExt[2], max(ya,yExt[3])
	yb = ϕ4(b)
	b>xExt[2] && return yExt[2], max(ya,yb)
	yb,ya
end



function SDrec(::Type{T},α,X,i1,i2,j1,j2)::T where T
	j1<i1 && return zero(T) # early out, skip lower triangular part


	if (i2-i1)<=10 || (j2-j1)<=10 # fallback - loop over the data points
		s = zero(T)
		# @inbounds for i in i1:i2, j in j1:j2
		# 	i>=j && continue # UGLY way to do upper triangular part. TODO: do it properly.
		# 	s += ϕ4((X[i]-X[j])/α)
		# end
		@inbounds for i in i1:i2
			for j in max(i+1,j1):j2 # upper triangular part
				s += ϕ4((X[i]-X[j])/α)
			end
		end
		return s
	end


	if i1!=j1# && (i2-i1)<20
		npoints = (i2-i1+1)*(j2-j1+1)
		# NB: these are not proper bounds since ϕ4 is not monotone (even for x>0)
		# lb = npoints*ϕ4((X[j2]-X[i1])/α)
		# ub = npoints*ϕ4((X[j1]-X[i2])/α)

		lb,ub = npoints.*ϕ4extrema((X[j1]-X[i2])/α, (X[j2]-X[i1])/α) # proper bounds

		# @show i1:i2,j1:j2,lb,ub
		# if abs(ub-lb) < 0.01 # TESTING. (abs() because not monotone which potentially creates weird issues)
		if ub-lb < 1 # 20 # 0.1 # 0.01 # what is a good bound?
			# @show i1:i2,j1:j2,lb,ub
			#return (i2-i1+1)*(j2-j1+1)*(lb+ub)/2 # approximation
			approx = (lb+ub)/2 # approximation
			return approx

			# # DEBUG
			# s = zero(T)
			# @inbounds for i in i1:i2
			# 	for j in max(i+1,j1):j2 # upper triangular part
			# 		s += ϕ4((X[i]-X[j])/α)
			# 	end
			# end

			# d = (abs(approx-s))
			# # d>1 && @show (i1,i2,j1,j2,d)
			# return approx
			# # return s
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
	# for yi in Y, yj in Y # TODO: fix inefficient summation.
	# 	s += ϕ4(yi-yj)
	# end
	@inbounds for i=1:n
		for j=i:n
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
	for yi in Y, yj in Y # TODO: fix inefficient summation.
		s += ϕ6(yi-yj)
	end
	C * s
end


function objective(h, X, α2Constant)
	# [ R(K) / (nσ⁴ₖSD(α₂(h))) ]^(1/5) - h = 0
	# R(K) = 1/(2√π)
	# σ⁴ₖ = 1
	n = length(X)
	α2 = α2Constant*h^(5/7) # 1.357[SD(a)/TD(b)]^(1/7) * h^(5/7)
	denom = (n*2*√π) * SD(α2,X)
	denom^(-1/5) - h
end


function _bwsj(X)
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

	find_zero(h->objective(h,X,α2Constant), (lower,upper), Bisection())
end


bwsj(X) = _bwsj(issorted(X) ? X : sort(X))

