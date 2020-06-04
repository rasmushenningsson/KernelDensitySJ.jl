

ϕ5(x) = (x2=x*x; x4=x2*x2; -x*(x4 - 10x2 + 15) * ϕ(x))


"""
Compute linear functions that bound ϕ⁽⁴⁾(x) from below and above on the given interval.
Returns kₗ,mₗ,kᵤ,mᵤ.
"""
function ϕ4affinebounds(a,b)::Tuple{Float64,Float64,Float64,Float64}
	# TODO: Use a cache instead of computing every time?
	@assert 0<=a<=b "$a, $b"

	x = (0.6167065901925941, 1.8891758777537109, 3.3242574335521193)

	if b-a < 1e-9 # to avoid div by zero
		y = ϕ4((a+b)/2)
		return 0.,y,0.,y
	end

	ai = (a>x[1]) + (a>x[2]) + (a>x[3])
	bi = (b>x[1]) + (b>x[2]) + (b>x[3])

	if ai!=bi # not in an interval where the function is convex/concave.
		# TODO: improve fallback.
		m,M = ϕ4extrema(a,b)
		return 0.,m,0.,M
	end

	ya,yb = ϕ4(a),ϕ4(b)

	# line above (convex) or below (concave) ϕ4 in the interval
	k1 = (yb-ya) / (b-a)
	m1 = ya - k1*a


	# tangent line at midpoint, i.e.
	# line below (convex) or above (concave) ϕ4 in the interval
	mid = (a+b)/2
	k2 = ϕ5(mid)
	m2 = ϕ4(mid) - k2*mid

	if ai&1 == 0 # concave
		k1,m1,k2,m2
	else # convex
		k2,m2,k1,m1
	end
end



# test
function testaffine(N)
	k = rand(1:4)
	x = (0.,0.6167065901925941, 1.8891758777537109, 3.3242574335521193, 10.)
	m = x[k]
	M = x[k+1]
	a,b = extrema(m.+(M-m)*rand(2))

	kl,ml,ku,mu = ϕ4affinebounds(a,b)

	Z = a .+ (b-a)*rand(N)

	Yl = kl*Z .+ ml
	Y  = ϕ4.(Z)
	Yu = ku*Z .+ mu

	@assert all(Yl .<= Y .+ 1e-9) "$a, $b"
	@assert all(Y .<= Yu .+ 1e-9) "$a, $b"
end
