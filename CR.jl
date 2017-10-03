import Krylov

# A truncated version of Stiefel’s Conjugate Residual method
# CR(A, b, Δ, atol, rtol, itmax, verbose) solves the linear system 'A * x = b' or the least-squares problem :
# 'min ‖b - A * x‖²' within a region of fixed radius Δ.

"""A truncated version of Stiefel’s Conjugate Residual method to solve the symmetric linear system Ax=b.
"""
function CR(A, b, Δ::Float64=10., atol::Float64=1.0e-8, rtol::Float64=1.0e-6, itmax::Int=0) #, verbose::Bool=true)
    n = size(b, 1) # size of the problem
    (size(A, 1) == n & size(A, 2) == n) || error("Inconsistent problem size")
    @info(logger, @sprintf("CR: system of %d equations in %d variables", n, n))

    x = zeros(n) # initial estimation x = 0
    xNorm = 0.0
    xNorms = [xNorm] # Values of ‖x‖
    r = b # initial residual r = b - Ax = b
    rNorm = norm(r, 2) # ‖r‖
    rNorm² = rNorm * rNorm
    s = A * r
    ρ = dot(r, s)
    absρ = abs(ρ)
    p = r
    q = s
    Δ² = Δ^2
    m = 0.0
    mvalues = [m] # values of the quadratic model
    ϵ = atol + rtol * rNorm
    pr = rNorm²
    abspr = pr
    pAp = ρ
    abspAp = absρ

    iter = 0
    itmax == 0 && (itmax = 2 * n)
    @info(logger, @sprintf("%5s %7s %7s %8s %8s %8s %8s %8s", "Iter", "‖x‖", "‖r‖", "q", "p'r", "α", "t1", "t2"))

    descent = pr > 0  # p'r > 0 means p is a descent direction
    solved = rNorm <= ϵ
    tired = iter >= itmax
    on_boundary = false

    while ! (solved || tired)
        info_line = @sprintf("%5d %7.1e %7.1e %8.1e %8.1e", iter, xNorm, rNorm, m, pr)
        iter += 1
        α = ρ / dot(q, q) # step

        if pAp ≤ 0 && Δ == 0
            @critical(logger, "indefinite system and no trust region")
            return x, p
        end

        if Δ > 0

            # find t1 > 0 and t2 < 0 such that ‖x + ti * p‖² = Δ²  (i = 1, 2)
            xNorm² = xNorm^2
            t = Krylov.to_boundary(x, p, Δ; flip = false, xNorm2 = xNorm²)
            t1 = maximum(t)
            t2 = minimum(t)

            # pour debugger
            @assert t1 > 0
            @assert t2 < 0

            if pAp ≤ 0
                @debug(logger, @sprintf("nonpositive curvature: pAp = %8.1e", pAp))

                # according to Fong and Saunders, p'r = 0 can only happen if pAp ≤ 0
                if abspr <= eps() * norm(p) * rNorm
                    @debug(logger, @sprintf("p'r = %8.1e ≃ 0", pr))

                    p = r # - ∇q(x)
                    pAp = dot(p, q)
                    abspAp = abs(pAp)
                    pr = abspr = rNorm²
                    descent = true

                    # TODO: find a way to not recompute t1 and t2
                    t = Krylov.to_boundary(x, p, Δ; flip = false, xNorm2 = xNorm²)
                    t1 = maximum(t)
                    t2 = minimum(t)

                    @assert t1 > 0
                    @assert t2 < 0

                    # TODO: continue normally if α < t1?
                    if α ≥ t1
                        α = t1
                        on_boundary = true
                    end

                else

                    α = descent ? t1 : t2
                    on_boundary = true

                end

            elseif α ≥ t1
                # at this point, it is not possible that α < 0 because pAp > 0
                # (cf. Fong and Saunders)
                @debug(logger, @sprintf("pAp = %8.1e > 0 but α = %8.1e ≥ t1 = %8.1e", pAp, α, t1))

                α = t1  # > 0
                on_boundary = true

            end

            info_line *= @sprintf(" %8.1e %8.1e %8.1e", α, t1, t2)

        end

        x = x + α * p
        xNorm = norm(x, 2)
        push!(xNorms, xNorm)
        Ax = A * x
        m = -dot(b, x) + 0.5 * dot(x, Ax)
        push!(mvalues, m)
        r = r - α * q  # residual
        rNorm = norm(r)

        @info(logger, info_line)

        solved = (rNorm <= ϵ) | on_boundary
        @debug(logger, @sprintf("solved = %s", solved))
        @debug(logger,@sprintf("on_boundary = %s", on_boundary))
        tired = iter >= itmax
        (solved || tired) && continue

        oldpAp = pAp
        s = A * r
        ρbar = ρ
        ρ = dot(r, s)
        absρ = abs(ρ)
        β = ρ / ρbar # step for the direction calculus
        p = r + β * p # search direction
        q = s + β * q

        pAp = dot(p, q)
        abspAp = abs(pAp)
        rNorm² = rNorm * rNorm
        pr = rNorm² + β * pr - β * α * oldpAp # p'r
        abspr = abs(pr)
        descent = pr > 0.0

    end
    # verbose && @printf("\n")

    return x
end
