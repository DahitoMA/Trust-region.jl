# A truncated version of Stiefel’s Conjugate Residual method
# cr(A, b, Δ, atol, rtol, verbose) solves the linear system 'A * x = b' or the least-squares problem :
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
  s = A * r
  ρ = dot(r, s)
  p = r
  P = [zeros(n), p] # [0, p_0]
  q = s
  Q = [zeros(n), q] # [0, q_0]
  Δ² = Δ^2
  m = 0.0
  mvalues = [m] # values of the quadratic model
  ϵ = atol + rtol * rNorm
  pr = rNorm * rNorm
  pAp = ρ

  iter = 0
  itmax == 0 && (itmax = 2* n)
  verbose && @printf("%5s %6s %10s %10s %10s %10s %10s\n", "Iter", "‖x‖", "‖r‖", "q", "α", "t1", "t2")
  verbose && @printf("    %d  %8.1e    %8.1e    %8.1e", iter, xNorm, rNorm, m)

  descent = pr > 0.0 # p'r > 0 means p is a descent direction
  solved = rNorm <= ϵ
  tired = iter >= itmax
  on_boundary = false

  while ! (solved || tired)
    iter += 1
    qNorm² = dot(q, q)
    α = ρ / qNorm² # step
    verbose && @printf("  %7.1e", α)

    if iter >= 2
      Q = [Q[2], q] # [q_{k-2}, q_{k-1}]
      P = [P[2], p] # [p_{k-2}, p_{k-1}]
    end

    if Δ > 0.0
      # solving ‖x+ti*p‖²-Δ² = 0 with i=1,2
      c = dot(p, p)
      a = dot(x, p)
      f = xNorm^2 - Δ²
      t = sqrt(a^2 - c * f)

      if a < 0.0
        t1 = (-a + t) / c
        t2 = f / (-a + t)
      else
        t1 = f / (-a - t)
        t2 = (-a - t) / c
      end

      if α <= eps(Float64) # 1.e-15
        γ = dot(s, A *  q) / qNorm²
        δ = dot(s, A * Q[1]) / dot(Q[1], Q[1])
        p = s - γ * p - δ * P[1]
        q = A * p
        α = dot(r, q) / dot(q, q)
        pAp = dot(p, q)
        pr = dot(p, r - α * q)
        descent = pr > 0.0
      end

      if abs(pr) < eps(Float64) # pr == 0.0 # p'r = 0
        solved = true
        continue
      elseif (!descent) & (pAp > 0.0) # p rise direction of positive curvature
        p = - p
        descent = true
      end

      verbose && @printf("   %7.1e   %7.1e\n", t1, t2)

      if pAp == 0.0 # p'q = 0
        if (ρ == 0.0) | (descent)
          α = t1 # > 0
          on_boundary = true
        elseif !descent
          α = t2 # < 0
          on_boundary = true
        end
      elseif (!descent) | (α < t2)
        α = t2 # < 0
        on_boundary = true
      elseif (pAp < 0.0) | (α > t1)
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
    m = dot(-b, x) + 0.5 * dot(x, Ax)
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
    β = ρ / ρbar # step for the direction calculus
    p = r + β * p # descent direction
    q = s + β * q

    pAp = dot(p, q)
    pr = rNorm * rNorm + β * pr - β * α * oldpAp # p'r
    descent = pr > 0.0

  end
  verbose && @printf("\n")

  return x
end
