import Krylov

# A truncated version of MINRES
# minres(A,b,Δ,M,λ,atol,rtol,etol,window,itmax,conlim) solves the linear
# system '(A + λ*M)x = b' or the least-squares problem : 'min ‖b - (A + λ*M)x‖²'
# within a region of radius Δ.

"""A truncated version of MINRES.
"""
function minres{T <: Number}(A :: AbstractLinearOperator, b :: Vector{T}, Δ::Float64=10.,
                             M :: AbstractLinearOperator=opEye(size(A,1)),
                             λ :: Float64=0.0,
                             atol :: Float64=1.0e-8, rtol :: Float64=1.0e-6,
                             etol :: Float64=0.0, window :: Int=5,
                             itmax :: Int=0, conlim :: Float64=1.0e+8)


    (m, n) = size(A)
    m == n || error("System must be square")
    size(b, 1) == m || error("Inconsistent problem size")
    @info(loggerminres, @sprintf("MINRES: system of size %d", n))
    ϵM = eps(T)
    x = zeros(T, n)
    ctol = conlim > 0.0 ? 1./conlim : 0.0;

    q = 0.0
    qvalues = [q] # values of the quadratic model

    # Initialize Lanczos process.
    # β₁ M v₁ = b.
    r1 = copy(b)
    y = M * r1
    β₁ = dot(r1, y)
    β₁ < 0.0 && error("Preconditioner is not positive definite")
    β₁ == 0.0 && return (x, "x = 0 is a zero-residual solution")
    β₁ = sqrt(β₁)
    β = β₁

    oldβ = 0.0
    δbar = 0.0
    ϵ = 0.0
    rNorm = β₁
    rNorms = [β₁]
    ϕbar = β₁
    rhs1 = β₁
    rhs2 = 0.0
    γmax = 0.0
    γmin = Inf
    cs = -1.0
    sn = 0.0
    v = zeros(T, n)
    w = zeros(T, n)
    w1 = zeros(T, n)
    w2 = zeros(T, n)
    r2 = copy(r1)

    ANorm² = 0.0
    ANorm = 0.0
    Acond = 0.0
    ArNorm = 0.0
    ArNorms = [ArNorm]
    xNorm = 0.0

    xENorm² = 0.0
    err_lbnd = 0.0
    err_vec = zeros(T, window)
    @debug(loggerminres, @sprintf("%5s  %7s  %7s  %7s  %8s  %8s  %7s  %7s",
                       "Aprod", "‖r‖", "‖A'r‖", "β", "cos", "sin", "‖A‖", "κ(A)"))
    @debug(loggerminres, @sprintf("%5d  %7.1e  %7.1e  %7.1e  %8.1e  %8.1e  %7.1e  %7.1e",
                       0, rNorm, ArNorm, β, cs, sn, ANorm, Acond))

    iter = 0
    itmax == 0 && (itmax = n)

    @info(loggerminres, @sprintf("%5s %8s %5s %5s", "Iter", "‖x‖", "‖r‖", "q"))
    @info(loggerminres, @sprintf("%5d %8.1e %8.1e %8.1e", iter, xNorm, rNorm, q))

    tol = atol + rtol * β₁
    status = "unknown"
    solved = solved_mach = solved_lim = (rNorm ≤ tol)
    tired  = iter ≥ itmax
    ill_cond = ill_cond_mach = ill_cond_lim = false
    zero_resid = zero_resid_mach = zero_resid_lim = (rNorm ≤ tol)
    fwd_err = false

    while ! (rNorm ≤ tol || tired) # ! (solved || tired || ill_cond)
        iter = iter + 1

        # Generate next Lanczos vector.
        v = 1./β * y
        y = A * v
        λ != 0.0 && (y = y - λ * v)
        iter ≥ 2 && (y = y - (β / oldβ) * r1)

        α = dot(v, y)
        y = y - (α / β) * r2
        r1 = copy(r2)
        r2 = copy(y)
        y = M * r2
        oldβ = β
        β = dot(r2, y)
        β < 0.0 && error("Preconditioner is not positive definite")
        β = sqrt(β)
        ANorm² = ANorm² + α * α + oldβ * oldβ + β * β

        # Apply rotation to obtain
        #  [ δₖ   ϵₖ₊₁    ] = [ cs  sn ] [ δbarₖ  0    ]
        #  [ γbar  δbarₖ₊₁ ]   [ sn -cs ] [ αₖ     βₖ₊₁ ]
        oldϵ = ϵ
        δ = cs * δbar + sn * α
        γbar = sn * δbar - cs * α
        ϵ = sn * β
        δbar = -cs * β
        root = sqrt(γbar * γbar + δbar * δbar)
        ArNorm = ϕbar * root  # = ‖A'rₖ₋₁‖
        push!(ArNorms, ArNorm)

        # Compute the next plane rotation.
        γ = sqrt(γbar * γbar + β * β)
        γ = max(γ, ϵM)
        cs = γbar / γ
        sn = β / γ
        ϕ = cs * ϕbar # step
        ϕbar = sn * ϕbar

        # Update x.
        w1 = copy(w2)
        w2 = copy(w)
        w = (v - oldϵ * w1 - δ * w2) / γ # descent direction

        if Δ > 0
            # find t1 > 0 and t2 < 0 such that ‖x + ti * w‖² = Δ² (i=1,2)
            xNorm² = xNorm^2
            t = Krylov.to_boundary(x, w, Δ; flip = false, xNorm2 = xNorm²)
            t1 = maximum(t)
            t2 = minimum(t)

            # if x is out of the trust region, w is followed until the edge of the
            #trust region
            if ϕ >= t1
                x = x + t1 * w
                q = -dot(b, x) + 0.5 * dot(x, A * x)
                qvalues = push!(qvalues, q)
                xNorm = norm(x, 2)
                rNorm = norm(A*x - b)
                @info(loggerminres, @sprintf("%5d %8.1e %8.1e %8.1e", iter, xNorm, rNorm, q))
                return x
            elseif ϕ <= t2
                x = x + t2 * w
                q = -dot(b, x) + 0.5 * dot(x, A * x)
                qvalues = push!(qvalues, q)
                xNorm = norm(x, 2)
                rNorm = norm(A*x - b)
                @info(loggerminres, @sprintf("%5d %8.1e %8.1e %8.1e", iter, xNorm, rNorm, q))
                return x
            end
        end

        x = x + ϕ * w
        xNorm = norm(x, 2)
        q = -dot(b, x) + 0.5 * dot(x, A * x)
        qvalues = push!(qvalues, q)

        xENorm² = xENorm² + ϕ * ϕ

        # Compute lower bound on forward error.
        err_vec[mod(iter, window) + 1] = ϕ
        iter ≥ window && (err_lbnd = norm(err_vec))

        γmax = max(γmax, γ)
        γmin = min(γmin, γ)
        ζ = rhs1 / γ
        rhs1 = rhs2 - δ * ζ
        rhs2 = -ϵ * ζ

        # Estimate various norms.
        ANorm = sqrt(ANorm²)

        ϵA = ANorm * ϵM
        ϵx = ANorm * xNorm * ϵM
        ϵr = ANorm * xNorm * rtol
        d = γbar
        d == 0.0 && (d = ϵA)

        rNorm = ϕbar

        test1 = rNorm / (ANorm * xNorm)
        test2 = root / ANorm
        push!(rNorms, rNorm)

        Acond = γmax / γmin

        @debug(loggerminres, @sprintf("%5d  %7.1e  %7.1e  %7.1e  %8.1e  %8.1e  %7.1e  %7.1e",
                         iter, test1, test2, β, cs, sn, ANorm, Acond))

        @info(loggerminres, @sprintf("%5d %8.1e %8.1e %8.1e", iter, xNorm, rNorm, q))

        # if iter == 1
        #     # A'b = 0 so x = 0 is a minimum least-squares solution
        #     β / β₁ ≤ 10 * ϵM && return x, "x = 0 is a minimum least-squares solution"
        # end

        # Stopping conditions that do not depend on user input.
        # This is to guard against tolerances that are unreasonably small.
        ill_cond_mach = (1.0 + 1.0 / Acond ≤ 1.0)
        solved_mach = (1.0 + test2 ≤ 1.0)
        zero_resid_mach = (1.0 + test1 ≤ 1.0)
        # solved_mach = (ϵx ≥ β₁)

        # Stopping conditions based on user-provided tolerances.
        tired = iter ≥ itmax
        ill_cond_lim = (1.0 / Acond ≤ ctol)
        solved_lim = (test2 ≤ tol)
        zero_resid_lim = (test1 ≤ tol)
        iter ≥ window && (fwd_err = err_lbnd ≤ etol * sqrt(xENorm²))

        ill_cond = ill_cond_mach | ill_cond_lim
        solved = solved_mach | solved_lim |zero_resid_mach |zero_resid_lim |fwd_err

    end

    tired         && (status = "maximum number of iterations exceeded")
    ill_cond_mach && (status = "condition number seems too large for this machine")
    ill_cond_lim  && (status = "condition number exceeds tolerance")
    solved        && (status = "found approximate minimum least-squares solution")
    zero_resid    && (status = "found approximate zero-residual solution")
    fwd_err       && (status = "truncated forward error small enough")

    return x
end
