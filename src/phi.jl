ϕ(x) = 1/√(2π) * exp(-x^2/2)

ϕ4(x) = (x2=x*x; (x2*x2 - 6x2 + 3) * ϕ(x))
# ϕ5(x) = (x2=x*x; x4=x2*x2; -x*(x4 - 10x2 + 15) * ϕ(x)) # For reference: used to find local extrema of ϕ4
ϕ6(x) = (x2=x*x; x4=x2*x2; (x4*x2 - 15x4 + 45x2 - 15) * ϕ(x))
# ϕ7(x) = (x2=x*x; x4=x2*x2; x6=x4*x2; x*(x6 - 21x4 + 105x2 - 105)*ϕ(x))                 # For reference: used to find local extrema of ϕ6
# ϕ8(x) = (x2=x*x; x4=x2*x2; x6=x4*x2; x8=x4*x4; (x8 - 28x6 + 210x4 - 420x2 + 105)*ϕ(x)) # For reference: used to find intervals where ϕ6 is convex or concave.

ϕ(::Val{4},x) = ϕ4(x)
ϕ(::Val{6},x) = ϕ6(x)


function ϕ4bounds(a,b,x)::Tuple{Float64,Float64}
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
		return ϕ4extrema(a,b)
	end

	ya,yb = ϕ4(a),ϕ4(b)

	# line above (convex) or below (concave) ϕ4 in the interval
	k = (yb-ya) / (b-a)
	m = ya - k*a
	y2 = (k*x+m)

	# lower bound (convex) or upper bound (concave) using Jensen's inqueality.
	y1 = ϕ4(x)


	if ai&1 == 0 # concave
		y2, y1
	else # convex
		y1, y2
	end
end

function ϕ6bounds(a,b,x)::Tuple{Float64,Float64}
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
		return ϕ6extrema(a,b)
	end

	ya,yb = ϕ6(a),ϕ6(b)

	# line above (convex) or below (concave) ϕ6 in the interval
	k = (yb-ya) / (b-a)
	m = ya - k*a
	y2 = k*x+m

	# lower bound (convex) or upper bound (concave) using Jensen's inqueality.
	y1 = ϕ6(x)


	if ai&1 == 1 # concave
		y2, y1
	else # convex
		y1, y2
	end
end

ϕbounds(::Val{4},a,b,x) = ϕ4bounds(a,b,x)
ϕbounds(::Val{6},a,b,x) = ϕ6bounds(a,b,x)




function ϕ4extrema(a,b)
	# (x-√(5+√10))(x+√(5+√10))(x-√(5-√10))(x+√(5-√10))x
	# -√(5+√BigFloat(10)), √(5-√BigFloat(10)), 0, √(5-√BigFloat(10)), √(5+√BigFloat(10))
	# xExt = (-2.8569700138728056, 1.355626179974266, 0.0, 1.355626179974266, 2.8569700138728056)
	# yExt = (0.13912184973492348, -0.7399861849949221, 1.1968268412042982, -0.7399861849949221, 0.13912184973492348)
	# we only need to consider x>0 since the function is even
	xExt = (1.355626179974266, 2.8569700138728056)
	yExt = (-0.7399861849949221, 0.13912184973492348)

	# the implementation assumes 0<=a<=b
	if b<=xExt[1] || a>=xExt[2]
		ϕ4(b),ϕ4(a) # decreasing function in these intervals
	elseif a>=xExt[1]
		if b<=xExt[2]
			ϕ4(a),ϕ4(b) # increasing function in this interval
		else
			min(ϕ4(a),ϕ4(b)),yExt[2] # the function covers a local maximum
		end
	elseif b<=xExt[2]
		yExt[1],max(ϕ4(a),ϕ4(b)) # the function covers the global minimum
	else
		yExt[1], max(ϕ4(a),yExt[2]) # the function covers the global minimum and a local maximum
	end
end

function ϕ6extrema(a,b)
	# ϕ6 has global minimum at 0, global maximum at x≈1.15, two more local extrema and ϕ6(x)→0 as x→∞.
	xExt = (1.154405394739968, 2.366759410734541, 3.750439717725742)
	yExt = (4.2406124820663305, -1.4018233232694841, 0.15246889807941977)
	# the implementation assumes 0<=a<=b


	if xExt[3]<=a || (xExt[1]<=a && b<=xExt[2])
		ϕ6(b),ϕ6(a) # decreasing function in these intervals
	elseif b<=xExt[1] || (xExt[2]<=a && b<=xExt[3])
		ϕ6(a),ϕ6(b) # increasing function in these intervals
	elseif b<=xExt[2]
		min(ϕ6(a),ϕ6(b)),yExt[1] # covers global maximum
	elseif xExt[1]<=a && b<=xExt[3]
		yExt[2],max(ϕ6(a),ϕ6(b)) # covers local minimum
	elseif xExt[2]<=a
		min(ϕ6(a),ϕ6(b)),yExt[3] # covers second local maximum
	elseif xExt[1]<=a
		yExt[2],max(ϕ6(a),yExt[3])  # covers the local minimum and the second local maximum
	# elseif b<=xExt[3] # same as last case
	# 	min(ϕ6(a),yExt[2]),yExt[1] # covers global maximum and local minimum
	else
		min(ϕ6(a),yExt[2]),yExt[1] # covers global maximum and local minimum
	end
end

ϕextrema(::Val{4},a,b) = ϕ4extrema(a,b)
ϕextrema(::Val{6},a,b) = ϕ6extrema(a,b)
