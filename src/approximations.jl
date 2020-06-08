function ϕ4bounds(a,b,x,npoints)::Tuple{Float64,Float64}
	# TODO: Use a cache instead of computing affine approximation every time?
	@assert 0<=a<=x<=b "$a, $x, $b"

	breakpoints = (0.6167065901925941, 1.8891758777537109, 3.3242574335521193)

	if b-a < 1e-9 # to avoid div by zero
		y = ϕ4((a+b)/2)
		return y,y
	end

	ai = (a>breakpoints[1]) + (a>breakpoints[2]) + (a>breakpoints[3])
	bi = (b>breakpoints[1]) + (b>breakpoints[2]) + (b>breakpoints[3])

	if ai!=bi # not in an interval where the function is convex/concave.
		# TODO: improve fallback.
		m,M = npoints.*ϕ4extrema(a,b)
		return m,M
	end

	ya,yb = ϕ4(a),ϕ4(b)

	# line above (convex) or below (concave) ϕ4 in the interval
	k = (yb-ya) / (b-a)
	m = ya - k*a
	y2 = npoints*(k*x+m)

	# lower bound (convex) or upper bound (concave) using Jensen's inqueality.
	y1 = ϕ4(x)*npoints


	if ai&1 == 0 # concave
		y2, y1
	else # convex
		y1, y2
	end
end



function ϕ6bounds(a,b,x,npoints)::Tuple{Float64,Float64}
	# TODO: Use a cache instead of computing affine approximation every time?
	@assert 0<=a<=x<=b "$a, $x, $b"

	breakpoints = (0.5390798113513751, 1.636519042435108, 2.8024858612875416, 4.144547186125894)
	if b-a < 1e-9 # to avoid div by zero
		y = ϕ6((a+b)/2)
		return y,y
	end

	ai = (a>breakpoints[1]) + (a>breakpoints[2]) + (a>breakpoints[3]) + (a>breakpoints[4])
	bi = (b>breakpoints[1]) + (b>breakpoints[2]) + (b>breakpoints[3]) + (b>breakpoints[4])

	if ai!=bi # not in an interval where the function is convex/concave.
		# TODO: improve fallback.
		m,M = npoints.*ϕ6extrema(a,b)
		return m,M
	end

	ya,yb = ϕ6(a),ϕ6(b)

	# line above (convex) or below (concave) ϕ6 in the interval
	k = (yb-ya) / (b-a)
	m = ya - k*a
	y2 = npoints*(k*x+m)

	# lower bound (convex) or upper bound (concave) using Jensen's inqueality.
	y1 = ϕ6(x)*npoints


	if ai&1 == 1 # concave
		y2, y1
	else # convex
		y1, y2
	end
end


ϕbounds(::Val{4},a,b,x,npoints) = ϕ4bounds(a,b,x,npoints)
ϕbounds(::Val{6},a,b,x,npoints) = ϕ6bounds(a,b,x,npoints)

