# A truncated version of Stiefel’s Conjugate Residual method
# CR(A, b, Δ, atol, rtol, itmax, verbose) solves the linear system 'A * x = b' or the least-squares problem :
# 'min ‖b - A * x‖²' within a region of fixed radius Δ.

"""A truncated version of Stiefel’s Conjugate Residual method to solve the symmetric linear system Ax=b.
"""
function CR(A, b, Δ::Float64=10., atol::Float64=1.0e-8, rtol::Float64=1.0e-6, itmax::Int=0, verbose::Bool=true)
  n = size(b, 1) # size of the problem
  (size(A, 1) == n & size(A, 2) == n) || error("Inconsistent problem size")
  verbose && @printf("CR: system of %d equations in %d variables\n", n, n)

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
  verbose && @printf("%5s %6s %10s %10s %10s %10s %10s\n", "Iter", "‖x‖", "‖r‖", "q", "α", "t1", "t2")
  verbose && @printf("    %d  %8.1e    %8.1e    %8.1e", iter, xNorm, rNorm, m)

  descent = pr > 0.0 # p'r > 0 means p is a descent direction
  solved = rNorm <= ϵ
  tired = iter >= itmax
  on_boundary = false

  while ! (solved || tired)
    iter += 1
    α = ρ / dot(q, q) # step
    verbose && @printf("  %7.1e", α)

    if Δ > 0.0
      # solving ‖x+ti*p‖²-Δ² = 0 with i=1,2
      xNorm² = xNorm^2
      t1 = to_boundary(x, p, Δ; flip = false, xNorm2 = xNorm²)
      t2 = - to_boundary(x, p, Δ; flip = true, xNorm2 = xNorm²)

      verbose && @printf("   %7.1e   %7.1e\n", t1, t2)

      if abspr <= ϵ # pr = 0
        p = r # - ∇q(x)
        pAp = dot(p, q)
        abspAp = abs(pAp)
        pr = rNorm²
        abspr = pr
        descent = true
      end

      if abspAp <= ϵ # p'q = 0

        if (absρ <= ϵ) | (descent) # ρ = 0 or descent = true
          α = t1 # > 0
          on_boundary = true
        elseif !descent # descent = false
          α = t2 # < 0
          on_boundary = true
        end

      elseif pAp < 0.0 # negative curvature
          if descent
              α = t1 # > 0
          else α = t2 # < 0
          end
          on_boundary = true
      elseif (!descent) & (α > 0)
        p = - p
        pr = - pr
        abspr = pr
        descent = true

        if α >= t1
          α = t1 # > 0
          on_boundary = true
        end

      elseif (!descent) | (α <= t2)
        α = t2 # < 0
        on_boundary = true
      elseif α >= t1
        α = t1 # > 0
        on_boundary = true
      end

    else
      verbose && @printf("\n")

    end

    x = x + α * p # new estimation
    xNorm = norm(x, 2)
    push!(xNorms, xNorm)
    Ax = A * x
    m = - dot(b, x) + 0.5 * dot(x, Ax)
    push!(mvalues, m)
    r = r - α * q # residual
    rNorm = norm(r, 2) # ‖r‖

    verbose && @printf("    %d  %8.1e    %8.1e    %8.1e", iter, xNorm, rNorm, m)
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
    descent = pr > 0.0

  end
  verbose && @printf("\n")

  return x
end
