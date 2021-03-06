# NB: ϕ is scaled by √2π compared to standard definition
ϕ(x) = exp(-x^2/2)
ϕ(x,D) = exp((D-x)*(D+x)/2) # rescaling that can be used for improved numerical precision

# ϕ2(x) = (x*x-1)*ϕ(x) # For reference: used to find intervals where ϕ is convex or concave.

ϕ4(x) = (x2=x*x; (x2*x2 - 6x2 + 3) * ϕ(x))
# ϕ5(x) = (x2=x*x; x4=x2*x2; -x*(x4 - 10x2 + 15) * ϕ(x)) # For reference: used to find local extrema of ϕ4
ϕ6(x) = (x2=x*x; x4=x2*x2; (x4*x2 - 15x4 + 45x2 - 15) * ϕ(x))
# ϕ7(x) = (x2=x*x; x4=x2*x2; x6=x4*x2; x*(x6 - 21x4 + 105x2 - 105)*ϕ(x))                 # For reference: used to find local extrema of ϕ6
# ϕ8(x) = (x2=x*x; x4=x2*x2; x6=x4*x2; x8=x4*x4; (x8 - 28x6 + 210x4 - 420x2 + 105)*ϕ(x)) # For reference: used to find intervals where ϕ6 is convex or concave.

ϕ(::Val{4},x) = ϕ4(x)
ϕ(::Val{6},x) = ϕ6(x)



"""
	ϕbounds(a,b,x,D)

Let `x` be the mean of a set of points {xᵢ} such that xᵢ∈[`a`,`b`] ∀i.
ϕbounds computes lower and upper bounds for ∑ᵢϕ(xᵢ)/n, where n is the number of points.
`D` is a rescaling that can be used for better numerical precision. ϕ(x) := 1/√(2π) * exp(D^2/2 - x^2/2)
"""
function ϕbounds(a,b,x,D)::Tuple{Float64,Float64}
	if b-a < 1e-9 # to avoid div by zero
		(abs(x-a)<1e-9 && abs(b-x)<1e-9) || throw(ArgumentError("x=$x must be in interval [a=$a, b=$b]"))
		y = ϕ((a+b)/2,D)
		return y,y
	end
	0<=a<=x<=b || throw(ArgumentError("x=$x must be in interval [a=$a, b=$b]"))

	ya,yb = ϕ(a,D),ϕ(b,D)
	if a<1 && b>1 # not in an interval where the function is convex/concave.
		# TODO: improve fallback
		return yb,ya
	end

	# bound by line above (convex) or below (concave) ϕ in the interval
	y2 = (yb-ya)*((x-a)/(b-a)) + ya

	# lower bound (convex) or upper bound (concave) using Jensen's inqueality.
	y1 = ϕ(x,D)

	if b<=1 # concave
		y2, y1
	else # convex
		y1, y2
	end
end


"""
	ϕ4bounds(a,b,x)

Let `x` be the mean of a set of points {xᵢ} such that xᵢ∈[`a`,`b`] ∀i.
ϕ4bounds computes lower and upper bounds for ∑ᵢϕ⁽⁴⁾(xᵢ)/n, where n is the number of points.
"""
function ϕ4bounds(a,b,x)::Tuple{Float64,Float64}
	breakpoints = (0.6167065901925941, 1.8891758777537109, 3.3242574335521193)
	if b-a < 1e-9 # to avoid div by zero
		(abs(x-a)<1e-9 && abs(b-x)<1e-9) || throw(ArgumentError("x=$x must be in interval [a=$a, b=$b]"))
		y = ϕ4((a+b)/2)
		return y,y
	end
	0<=a<=x<=b || throw(ArgumentError("x=$x must be in interval [a=$a, b=$b]"))

	ai = (a>breakpoints[1]) + (a>breakpoints[2]) + (a>breakpoints[3])
	bi = (b>breakpoints[1]) + (b>breakpoints[2]) + (b>breakpoints[3])

	if ai!=bi # not in an interval where the function is convex/concave.
		# TODO: improve fallback.
		return ϕ4extrema(a,b)
	end

	ya,yb = ϕ4(a),ϕ4(b)

	# bound by line above (convex) or below (concave) ϕ4 in the interval
	y2 = (yb-ya)*((x-a)/(b-a)) + ya


	# lower bound (convex) or upper bound (concave) using Jensen's inqueality.
	y1 = ϕ4(x)


	if ai&1 == 0 # concave
		y2, y1
	else # convex
		y1, y2
	end
end

"""
	ϕ6bounds(a,b,x)

Let `x` be the mean of a set of points {xᵢ} such that xᵢ∈[`a`,`b`] ∀i.
ϕ6bounds computes lower and upper bounds for ∑ᵢϕ⁽⁶⁾(xᵢ)/n.
"""
function ϕ6bounds(a,b,x)::Tuple{Float64,Float64}
	breakpoints = (0.5390798113513751, 1.636519042435108, 2.8024858612875416, 4.144547186125894)
	if b-a < 1e-9 # to avoid div by zero
		(abs(x-a)<1e-9 && abs(b-x)<1e-9) || throw(ArgumentError("x=$x must be in interval [a=$a, b=$b]"))
		y = ϕ6((a+b)/2)
		return y,y
	end
	0<=a<=x<=b || throw(ArgumentError("x=$x must be in interval [a=$a, b=$b]"))

	ai = (a>breakpoints[1]) + (a>breakpoints[2]) + (a>breakpoints[3]) + (a>breakpoints[4])
	bi = (b>breakpoints[1]) + (b>breakpoints[2]) + (b>breakpoints[3]) + (b>breakpoints[4])

	if ai!=bi # not in an interval where the function is convex/concave.
		# TODO: improve fallback.
		return ϕ6extrema(a,b)
	end

	ya,yb = ϕ6(a),ϕ6(b)

	# bound by line above (convex) or below (concave) ϕ6 in the interval
	y2 = (yb-ya)*((x-a)/(b-a)) + ya

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
	# x at extrema: -√(5+√10), √(5-√10), 0, √(5-√10), √(5+√10)
	# we only need to consider x>0 since the function is even
	xExt = (1.355626179974266, 2.8569700138728056)
	yExt = (-1.8548702941445978, 0.34872676216452453)

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
	# we only need to consider x>0 since the function is even
	xExt = (1.154405394739968, 2.366759410734541, 3.750439717725742)
	yExt = (10.62963914930061, -3.513849978144481, 0.3821828509277058)
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
