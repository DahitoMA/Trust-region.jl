import Krylov

# A truncated version of Stiefel’s Conjugate Residual method
# CR(A, b, Δ, atol, rtol, itmax, verbose) solves the linear system 'A * x = b' or the least-squares problem :
# 'min ‖b - A * x‖²' within a region of fixed radius Δ.

"""A truncated version of Stiefel’s Conjugate Residual method to solve the symmetric linear system Ax=b.
"""
function CR(A, b, Δ::Float64=10., atol::Float64=1.0e-8, rtol::Float64=1.0e-6, itmax::Int=0) #, verbose::Bool=true)
    n = size(b, 1) # size of the problem
    (size(A, 1) == n & size(A, 2) == n) || error("Inconsistent problem size")
    @info(loggerCR, @sprintf("CR: system of %d equations in %d variables", n, n))

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
    @info(loggerCR, @sprintf("%5s %7s %7s %8s %8s %8s %8s %8s", "Iter", "‖x‖", "‖r‖", "q", "p'r", "α", "t1", "t2"))
    @info(loggerCR, @sprintf("%5d %7.1e %7.1e %8.1e %8.1e", iter, xNorm, rNorm, m, pr))

    descent = pr > 0  # p'r > 0 means p is a descent direction
    solved = rNorm ≤ ϵ
    tired = iter ≥ itmax
    on_boundary = false

    while ! (solved || tired)
        iter += 1
        α = ρ / dot(q, q) # step

        if pAp ≤ 0 && Δ == 0
            @critical(loggerCR, "indefinite system and no trust region")
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

            if abspAp ≤ eps() * norm(p) * norm(q)
                @debug(loggerCR, @sprintf("p'Ap = %8.1e ≃ 0", pAp))
                # according to Fong and Saunders, p'r = 0 can only happen if pAp ≤ 0
                if abspr ≤ eps() * norm(p) * rNorm
                    @debug(loggerCR, @sprintf("p'r = %8.1e ≃ 0, redefining p := r", pr))

                    p = r # - ∇q(x)
                    # q(x + αr) = q(x) - α ‖r‖² + ½ α² r'Ar
                    # 1) if rAr > 0, q decreases from α = 0 to α = ‖r‖² / r'Ar
                    # 2) if rAr ≤ 0, q decreases to -∞ in the direction r

                    q = s # = Ar = Ap
                    pAp = ρ # = dot(p, q) = pAp = rAr
                    abspAp = abs(pAp)
                    pr = abspr = rNorm²
                    descent = true

                    # TODO: find a way to not recompute t1 and t2
                    t = Krylov.to_boundary(x, p, Δ; flip = false, xNorm2 = xNorm²)
                    t1 = maximum(t)
                    # t2 = minimum(t)

                    # @assert t1 > 0
                    # @assert t2 < 0

                    if ρ > 0  # case 1
                        @debug(loggerCR,
                               @sprintf("quadratic is convex in direction r, curv = %7.1e", ρ))

                        α = rNorm² / ρ

                        if α ≥ t1
                            α = t1
                            on_boundary = true
                        end
                        solved = true

                    else  # case 2
                        @debug(loggerCR,
                               @sprintf("r is a direction of nonpositive curvature: %8.1e", ρ))

                        α = t1
                        on_boundary = true
                    end
                elseif descent
                    α = t1
                    on_boundary = true
                else
                    α = t2
                    on_boundary = true
                end

            elseif pAp > 0 && ρ > 0
                @debug(loggerCR, @sprintf("positive curvatures along p and r. p'Ap = %8.1e and r'Ar = %8.1e", pAp, ρ))
                if α ≥ t1
                    α = t1
                    on_boundary = true
                end

            elseif pAp > 0 && ρ < 0
                @debug(loggerCR, @sprintf("p'Ap = %8.1e > 0 and r'Ar = %8.1e < 0", pAp, ρ))
                p = r
                q = s # = Ar = Ap
                pAp = ρ # = dot(p, q) = pAp = rAr
                abspAp = abs(pAp)
                pr = abspr = rNorm²
                descent = true

                t = Krylov.to_boundary(x, p, Δ; flip = false, xNorm2 = xNorm²)
                α = t1 = maximum(t)
                on_boundary = true

            elseif pAp < 0 && ρ > 0
                @debug(loggerCR, @sprintf("p'Ap = %8.1e < 0 and r'Ar = %8.1e > 0", pAp, ρ))
                α = descent ? t1 : t2
                on_boundary = true

            elseif pAp < 0 && ρ < 0
                @debug(loggerCR, @sprintf("negative curvatures along p and r. p'Ap = %8.1e and r'Ar = %8.1e ", pAp, ρ))
                # q_p = q(x + ti * p) - q(x) = -ti * r'p + ½ (ti)² * p'Ap, i = 1, 2
                      # i = 1 if p is a descent direction and t2 otherwise
                # q_r = q(x + t1 * r) - q(x) = -tr * ‖r‖² + ½ (tr)² * r'Ar
                # dif = q_p - q_r
                # if dif > 0, r is followed until the edge of the trust-region
                # else p is followed until the edge of the trust-region

                t = Krylov.to_boundary(x, r, Δ; flip = false, xNorm2 = xNorm²)
                tr = maximum(t)

                if descent
                    α = t1
                else
                    α = t2
                end

                dif = -α * pr + tr * rNorm² + 0.5 * (α^2 * pAp - (tr)^2 * ρ)

                if dif > 0
                    @debug(loggerCR, @sprintf("direction r engenders a bigger decrease. q_p - q_r = %8.1e > 0", dif))
                    @debug(loggerCR, "redefining p := r")
                    p = r
                    q = s # = Ar = Ap
                    pAp = ρ # = dot(p, q) = pAp = rAr
                    abspAp = abs(pAp)
                    pr = abspr = rNorm²
                    descent = true
                    α = tr
                else @debug(loggerCR, @sprintf("direction p engenders an equal or a bigger decrease. q_p - q_r = %8.1e ≤ 0", dif))
                end
                on_boundary = true

            end

        end

        x = x + α * p
        xNorm = norm(x, 2)
        push!(xNorms, xNorm)
        Ax = A * x
        m = -dot(b, x) + 0.5 * dot(x, Ax)
        push!(mvalues, m)
        r = r - α * q  # residual
        rNorm = norm(r)

        @info(loggerCR, @sprintf("%5d %7.1e %7.1e %8.1e %8.1e %8.1e %8.1e %8.1e", iter, xNorm, rNorm, m, pr, α, t1, t2))

        solved = (rNorm <= ϵ) | on_boundary
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
        descent = pr > 0

    end

    return x
end
