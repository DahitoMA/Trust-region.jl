import Krylov

# A truncated version of Stiefel’s Conjugate Residual method
# CR(A, b, Δ, itmax; quad) solves the linear system 'A * x = b' or the least-squares problem :
# 'min ‖b - A * x‖²' within a region of fixed radius Δ.

"""A truncated version of Stiefel’s Conjugate Residual method to solve the symmetric linear system Ax=b.
If quad = true, the values of the quadratic model are computed.
"""
function CR(A, b, Δ::Float64=10., ϵa::Float64=1e-8, ϵr::Float64=1e-6, itmax::Int=0, ε::Float64=1e-6; quad::Bool=false)
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
    p = r
    q = s
    if quad
        m = 0.0
        mvalues = [m] # values of the quadratic model
        mstr = @sprintf("%8.1e", m)
    else
        mstr = ""
    end
    ϵ = ϵa + ϵr * rNorm # tolerance
    pr = rNorm²
    abspr = pr
    pAp = ρ # = dot(p, q) = dot(r, s)
    abspAp = abs(pAp)

    iter = 0
    itmax == 0 && (itmax = 2 * n)
    @info(loggerCR, @sprintf("%5s %7s %7s %8s %8s %8s %8s %8s", "Iter", "‖x‖", "‖r‖", "q", "p'r", "α", "t1", "t2"))
    @info(loggerCR, @sprintf("%5d %7.1e %7.1e %9s %8.1e", iter, xNorm, rNorm, mstr, pr))

    descent = pr > 0  # p'r > 0 means p is a descent direction
    solved = rNorm ≤ ϵ
    tired = iter ≥ itmax
    on_boundary = false
    npcurv = false

    while ! (solved || tired)
        iter += 1

        if pAp ≤ 0 && Δ == 0
            @critical(loggerCR, "indefinite system and no trust region")
            return x, p
        end

        if Δ > 0

            # find t1 > 0 and t2 < 0 such that ‖x + ti * p‖² = Δ²  (i = 1, 2)
            xNorm² = xNorm^2
            @debug(loggerCR, @sprintf("Δ = %8.1e > 0 and ‖x‖² = %8.1e", Δ, xNorm²))
            pNorm = norm(p)
            t = Krylov.to_boundary(x, p, Δ; flip = false, xNorm2 = xNorm², dNorm = pNorm)
            t1 = maximum(t)
            t2 = minimum(t)
            tr = maximum(Krylov.to_boundary(x, r, Δ; flip = false, xNorm2 = xNorm², dNorm = rNorm))
            @debug(loggerCR, @sprintf("t1 = %8.1e, t2 = %8.1e and tr = %8.1e", t1, t2, tr))

            # pour debugger
            @assert t1 > 0
            @assert t2 < 0

            if abspAp ≤ ε * pNorm * norm(q) # p'Ap ≃ 0
                @debug(loggerCR, @sprintf("p'Ap = %8.1e ≃ 0", pAp))
                # according to Fong and Saunders, p'r = 0 can only happen if pAp ≤ 0
                if abspr ≤ ε * pNorm * rNorm # p'r ≃ 0
                    @debug(loggerCR, @sprintf("p'r = %8.1e ≃ 0, redefining p := r", pr))

                    p = r # - ∇q(x)
                    # q(x + αr) = q(x) - α ‖r‖² + ½ α² rᵀAr
                    # 1) if rᵀAr > 0, the quadratic decreases from α = 0 to α = ‖r‖² / rᵀAr
                    # 2) if rᵀAr ≤ 0, the quadratic decreases to -∞ in the direction r

                    q = s # = Ar = Ap
                    pr = abspr = rNorm²
                    descent = true

                    if ρ > 0  # case 1
                        @debug(loggerCR,
                                @sprintf("quadratic is convex in direction r, curv = %8.1e", ρ))

                        α = rNorm² / ρ

                        if α ≥ tr
                            α = tr
                            on_boundary = true
                        end
                        solved = true

                    else  # case 2
                        @debug(loggerCR,
                                @sprintf("r is a direction of nonpositive curvature: %8.1e", ρ))

                        α = tr
                        on_boundary = true
                    end
                else
                    # q_p = q(x + α_p * p) - q(x) = -α_p * r'p + ½ (α_p)² * p'Ap
                    # q_r = q(x + α_r * r) - q(x) = -α_r * ‖r‖² + ½ (α_r)² * r'Ar
                    # dif = q_p - q_r
                    # if dif > 0, r is followed in the limits of the trust-region
                    # else p is followed in the limits of the trust-region

                    α = descent ?  t1 : t2
                    ρ > 0 && (tr = min(tr, rNorm² / ρ))
                    dif = -α * pr + tr * rNorm² - 0.5 * (tr)^2 * ρ # as pAp = 0
                    if dif > 0
                        @debug(loggerCR, @sprintf("direction r engenders a bigger decrease. q_p - q_r = %8.1e > 0", dif))
                        @debug(loggerCR, "redefining p <- r")
                        p = r
                        q = s # = Ar = Ap
                        descent = true
                        α = tr
                    else
                        @debug(loggerCR, @sprintf("direction p engenders an equal or a bigger decrease. q_p - q_r = %8.1e ≤ 0", dif))
                    end
                    on_boundary = true
                end

            elseif pAp > 0 && ρ > 0
                @debug(loggerCR, @sprintf("positive curvatures along p and r. p'Ap = %8.1e and r'Ar = %8.1e", pAp, ρ))
                α = ρ / dot(q, q) # step
                if α ≥ t1
                    α = t1
                    on_boundary = true
                end

            elseif pAp > 0 && ρ < 0
                @debug(loggerCR, @sprintf("p'Ap = %8.1e > 0 and r'Ar = %8.1e < 0", pAp, ρ))
                # q_p = q(x + α_p * p) - q(x) = -α_p * r'p + ½ (α_p)² * p'Ap
                # q_r = q(x + α_r * r) - q(x) = -α_r * ‖r‖² + ½ (α_r)² * r'Ar
                # q_p is minimal for α_p = rᵀp / pᵀAp
                # if pᵀr>0, the condition on α_p is: 0 < α_p ≤ t1
                # else t2 ≦ α_p < 0
                # q_r is minimal for α_r = ‖r²‖ / rᵀAr
                # the condition on α_r is: α_r ≦ t1
                # dif = q_p - q_r
                # if dif > 0, r is followed in the limits of the trust-region
                # else p is followed in the limits of the trust-region
                α = descent ?  min(t1, pr / pAp) : max(t2, pr / pAp)
                dif = -α * pr + tr * rNorm² + 0.5 * (α^2 * pAp - (tr)^2 * ρ)

                if dif > 0
                    @debug(loggerCR, @sprintf("direction r engenders a bigger decrease. q_p - q_r = %8.1e > 0", dif))
                    @debug(loggerCR, "redefining p := r")
                    p = r
                    q = s # = Ar = Ap
                    descent = true
                    α = tr
                else
                    @debug(loggerCR, @sprintf("direction p engenders an equal or a bigger decrease. q_p - q_r = %8.1e ≤ 0", dif))
                end
                on_boundary = true

            elseif pAp < 0 && ρ > 0
                @debug(loggerCR, @sprintf("p'Ap = %8.1e < 0 and r'Ar = %8.1e > 0", pAp, ρ))
                # pAp < 0 => the minimum of the quadratic along p is on the boundary
                α = descent ? t1 : t2
                tr = min(tr, rNorm² / ρ)
                dif = -α * pr + tr * rNorm² + 0.5 * (α^2 * pAp - (tr)^2 * ρ)

                if dif > 0
                    @debug(loggerCR, @sprintf("direction r engenders a bigger decrease. q_p - q_r = %8.1e > 0", dif))
                    @debug(loggerCR, "redefining p := r")
                    p = r
                    q = s # = Ar = Ap
                    descent = true
                    α = tr
                else
                    @debug(loggerCR, @sprintf("direction p engenders an equal or a bigger decrease. q_p - q_r = %8.1e ≤ 0", dif))
                end
                on_boundary = true

            elseif pAp < 0 && ρ < 0
                @debug(loggerCR, @sprintf("negative curvatures along p and r. p'Ap = %8.1e and r'Ar = %8.1e ", pAp, ρ))
                α = descent ?  t1 : t2
                dif = -α * pr + tr * rNorm² + 0.5 * (α^2 * pAp - (tr)^2 * ρ)

                if dif > 0
                    @debug(loggerCR, @sprintf("direction r engenders a bigger decrease. q_p - q_r = %8.1e > 0", dif))
                    @debug(loggerCR, "redefining p := r")
                    p = r
                    q = s # = Ar = Ap
                    descent = true
                    α = tr
                else
                    @debug(loggerCR, @sprintf("direction p engenders an equal or a bigger decrease. q_p - q_r = %8.1e ≤ 0", dif))
                end
                on_boundary = true

            end

        elseif Δ == 0
            α = ρ / dot(q, q)
        end

        x = x + α * p
        xNorm = norm(x, 2)
        xNorm ≈ Δ && (on_boundary = true)
        push!(xNorms, xNorm)
        if quad
            m = -dot(b, x) + 0.5 * dot(x, A * x)
            push!(mvalues, m)
            mstr = @sprintf("%8.1e", m)
        end
        r = r - α * q  # residual
        rNorm² = abs(rNorm² - α * ρ)
        rNorm = sqrt(rNorm²)

        @info(loggerCR, @sprintf("%5d %7.1e %7.1e %9s %8.1e %8.1e %8.1e %8.1e", iter, xNorm, rNorm, mstr, pr, α, t1, t2))

        solved = (rNorm <= ϵ) || on_boundary || npcurv
        tired = iter >= itmax
        (solved || tired) && continue

        s = A * r
        ρbar = ρ
        ρ = dot(r, s)
        β = ρ / ρbar # step for the direction calculus
        p = r + β * p # search direction
        q = s + β * q

        pr = rNorm² + β * pr - β * α * pAp # p'r
        pAp = ρ + β^2 * pAp # dot(p, q)
        abspAp = abs(pAp)
        abspr = abs(pr)
        descent = pr > 0

    end

    return x
end
