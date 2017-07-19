# A truncated version of the conjugate gradient method.
# CG(A,b,Δ,atol,rtol,verbose) solves the linear system 'A * x = b' within a region of radius Δ.

"""A truncated version of the conjugate gradient method to solve the symmetric linear system Ax=b.
A can be positive definite or not.
"""
function CG(A, b, Δ::Float64=10., atol::Float64=1e-8, rtol::Float64=1e-6, verbose::Bool=false)

  x = zeros(length(b)); # initial estimation x = 0
  xNorm = 0.0
  xNorms = [xNorm] # Values of ‖x‖
	r = -b # initial residual r = Ax-b = -b
	d = b # first descent direction
	Δ² = Δ^2
	rNorm = norm(r, 2)
	ϵ = atol + rtol * rNorm
  q = 0.0
	qvalues = [q] # values of the quadratic model

  iter = 0
	verbose && @printf("%5s %5s %14s %12s\n", "Iter", "‖x‖", "‖r‖", "q")
  verbose && @printf("    %d  %8.1e    %8.1e    %8.1e", iter, xNorm, rNorm, q)

	while rNorm > ϵ # stopping criterion : ‖r‖ < = ϵ
    iter += 1
		Ad = A * d
		dAd = dot(d, Ad)
		# solving ‖x+t1d‖²-Δ² = 0 with t1 >= 0
		c = norm(d, 2)^2
		a = dot(x, d)
		f = norm(x, 2)^2 - Δ²
		t = sqrt(a^2 - c * f)

		if a < 0.0
			t1 = (-a + t) / c
		else
			t1 = f / (-a - t)
		end

		# convexity test along the direction d : if the model is not convexe,
		#d is followed until the edge of the trust region
		if dAd <= 0.0
			x = x + t1 * d
			xNorm = norm(x,2)
      xNorms = push!(xNorms, xNorm)
			q = dot(-b, x) + 1/2 * dot(x, A * x)
			push!(qvalues, q)
			rNorm = norm(A * x - b)
			verbose && @printf("%d %8.1e %8.1e %8.1e\n", iter, xNorm, rNorm, q)
			return x
		end

		α = rNorm^2 / dAd # step for x estimation
		# if x is out of the trust region, d is followed until the edge of the
		#trust region
		if α >= t1
			x = x + t1 * d
			xNorm = norm(x, 2)
	    push!(xNorms, xNorm)
			q = dot(-b, x) + 1/2 * dot(x, A * x)
			push!(qvalues, q)
			rNorm = norm(A * x - b)
			verbose && @printf("%d %8.1e %8.1e %8.1e\n", iter, xNorm, rNorm, q)
			return x
		end

		x = x + α * d # new estimation
		xNorm = norm(x, 2)
    push!(xNorms, xNorm)
		q = dot(-b, x) + 0.5 * dot(x, A * x)
		push!(qvalues, q)
		roldNorm = rNorm
		r = r + α * Ad # new residual
		rNorm = norm(r, 2)

		β = rNorm^2 / roldNorm^2 # step for the next descent direction
		d = -r + β * d # new descent direction
		verbose && @printf("%d %8.1e %8.1e %8.1e\n", iter, xNorm, rNorm, q)

	end
	return x
end
