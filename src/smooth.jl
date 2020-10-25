# Work in progress. Only for one point so far.
function gaussiansmoothing(x, y, bandwidth, xeval::Number; leafSize=10, rtol=1e-3)
	@assert length(x)==length(y)
	c = 1/(2*bandwidth^2)

	if !issorted(x)
		perm = sortperm(x)
		x = x[perm]
		y = y[perm]
	end

	# 


end
