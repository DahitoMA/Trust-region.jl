import Krylov

# A truncated version of the conjugate gradient method.
# CG(A,b,Δ,itmax; quad) solves the linear system 'A * x = b' within a region of radius Δ.

"""A truncated version of the conjugate gradient method to solve the symmetric linear system Ax=b.
If quad = true, the values of the quadratic model are computed.
A can be positive definite or not.
"""
function CG(A, b, Δ::Float64=10., ϵa::Float64=1e-8, ϵr::Float64=1e-6, itmax::Int=0, ε::Float64=1e-6; quad::Bool=false)
    n = size(b, 1) # size of the problem
    (size(A, 1) == n & size(A, 2) == n) || error("Inconsistent problem size")
    @info(loggerCG, @sprintf("CG: system of %d equations in %d variables", n, n))

    x = zeros(n) # initial estimation x = 0
    xNorm = 0.0
    xNorms = [xNorm] # Values of ‖x‖
    r = -b # initial residual r = Ax-b = -b
    d = b # first descent direction
    rNorm = norm(r, 2)
    ϵ = ϵa + ϵr * rNorm # tolerance
    if quad
        q = 0.0
        qvalues = [q] # values of the quadratic model
        qstr = @sprintf("%8.1e", q)
    else
        qstr = ""
    end


    iter = 0
    itmax == 0 && (itmax = 2 * n)
    @info(loggerCG, @sprintf("%5s %6s %10s %10s %10s %10s\n", "Iter", "‖x‖", "‖r‖", "q", "α", "t1"))
    @info(loggerCG, @sprintf("    %d  %8.1e    %8.1e    %5s", iter, xNorm, rNorm, qstr))

    solved = rNorm ≤ ϵ
    tired = iter ≥ itmax
    on_boundary = false

    while ! (solved || tired)
        iter += 1
        Ad = A * d
        dAd = dot(d, Ad)

        if dAd ≤ 0 && Δ == 0
            @critical(loggerCG, "indefinite system and no trust region")
            return x, d
        end

        if Δ > 0
          # find t1 > 0 such that ‖x + t1 * d‖² = Δ²
            xNorm² = xNorm^2
            @debug(loggerCG, @sprintf("Δ = %8.1e > 0 and ‖x‖² = %8.1e", Δ, xNorm²))
            t1 = maximum(Krylov.to_boundary(x, d, Δ; flip = false, xNorm2 = xNorm²))
            @assert t1 > 0

            # if the model is not convexe or x is out of the trust region,
            # d is followed until the edge of the trust region
            if dAd ≤ ε * norm(d) * norm(Ad)
                @debug(loggerCG, @sprintf("non positive curvature dAd = %8.1e", dAd))
                α = t1
                on_boundary = true
            else
                α = rNorm^2 / dAd # step for x estimation
                if α ≥ t1
                    @debug(loggerCG, @sprintf("α = %8.1e ≥ t1 = %8.1e", α, t1))
                    α = t1
                    on_boundary = true
                end
            end

        end

        x = x + α * d # new estimation
        xNorm = norm(x, 2)
        push!(xNorms, xNorm)

        if quad
            q = -dot(b, x) + 0.5 * dot(x, A * x)
            push!(qvalues, q)
            qsrt = @sprintf("%8.1e", q)
        end
        roldNorm = rNorm
        r = r + α * Ad # new residual
        rNorm = norm(r, 2)

        @info(loggerCG, @sprintf("    %d  %8.1e    %8.1e    %5s  %8.1e  %8.1e", iter, xNorm, rNorm, qstr, α, t1))

        solved = (rNorm ≤ ϵ) | on_boundary
        tired = iter ≥ itmax

        (solved || tired) && continue
        β = rNorm^2 / roldNorm^2 # step for the next descent direction
        d = -r + β * d # new descent direction

        end
        return x
    end
