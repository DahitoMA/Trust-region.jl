import Krylov

# A truncated version of the conjugate gradient method.
# CG(A,b,Δ,atol,rtol,verbose) solves the linear system 'A * x = b' within a region of radius Δ.

"""A truncated version of the conjugate gradient method to solve the symmetric linear system Ax=b.
A can be positive definite or not.
"""
function CG(A, b, Δ::Float64=10., atol::Float64=1e-8, rtol::Float64=1e-6, itmax::Int=0, verbose::Bool=false)
    n = size(b, 1) # size of the problem
    (size(A, 1) == n & size(A, 2) == n) || error("Inconsistent problem size")
    verbose && @printf("CG: system of %d equations in %d variables\n", n, n)

    x = zeros(n) # initial estimation x = 0
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
    itmax == 0 && (itmax = 2 * n)
    verbose && @printf("%5s %6s %10s %10s %10s %10s\n", "Iter", "‖x‖", "‖r‖", "q", "α", "t1")
    verbose && @printf("    %d  %8.1e    %8.1e    %8.1e", iter, xNorm, rNorm, q)

    solved = rNorm <= ϵ
    tired = iter >= itmax
    on_boundary = false

    while ! (solved || tired)
        iter += 1
        Ad = A * d
        dAd = dot(d, Ad)
        α = rNorm^2 / dAd # step for x estimation
        verbose && @printf("  %7.1e", α)

        if Δ > 0.0
            # solving ‖x+t1d‖²-Δ² = 0 with t1 >= 0
            xNorm² = xNorm^2
            t1 = maximum(Krylov.to_boundary(x, d, Δ; flip = false, xNorm2 = xNorm²))

            verbose && @printf("   %7.1e\n", t1);

            # if the model is not convexe or x is out of the trust region,
            # d is followed until the edge of the trust region
            if (dAd <= 0.0) | (α >= t1)
                α = t1
                on_boundary = true
            end

        else
            verbose && @printf("\n")

        end

        x = x + α * d # new estimation
        xNorm = norm(x, 2)
        push!(xNorms, xNorm)
        q = dot(-b, x) + 0.5 * dot(x, A * x)
        push!(qvalues, q)
        roldNorm = rNorm
        r = r + α * Ad # new residual
        rNorm = norm(r, 2)

        verbose && @printf("    %d  %8.1e    %8.1e    %8.1e", iter, xNorm, rNorm, q)
        solved = (rNorm <= ϵ) | on_boundary
        tired = iter >= itmax

        (solved || tired) && continue
        β = rNorm^2 / roldNorm^2 # step for the next descent direction
        d = -r + β * d # new descent direction

        end
        verbose && @printf("\n")
        return x
    end
