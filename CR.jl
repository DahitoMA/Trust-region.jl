# A truncated version of Stiefel’s Conjugate Residual method
# cr(A, b, Δ, atol, rtol, verbose) solves the linear system 'A * x = b' or the least-squares problem :
# 'min ‖b - A * x‖²' within a region of fixed radius Δ.

"""A truncated version of Stiefel’s Conjugate Residual method to solve the symmetric linear system Ax=b.
"""
function CR(A, b, Δ::Float64=10., atol::Float64=1.0e-8, rtol::Float64=1.0e-6, itmax::Int=0, verbose::Bool=false)
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
  q = s
  Δ² = Δ^2
  m = 0.0
  mvalues = [m] # values of the quadratic model
  ϵ = atol + rtol * rNorm
  γ = rNorm * rNorm
  pAp = ρ

  iter = 0
  itmax == 0 && (itmax = 2* n)
  verbose && @printf("%5s %6s %10s %10s %10s %10s\n", "Iter", "‖x‖", "‖r‖", "q", "α", "t1")
  verbose && @printf("    %d  %8.1e    %8.1e    %8.1e", iter, xNorm, rNorm, m)

  solved = rNorm <= ϵ
  tired = iter >= itmax
  on_boundary = false

  while ! (solved || tired)
    iter += 1
    α = ρ / dot(q, q) # step
    verbose && @printf("  %7.1e", α)

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

      verbose && @printf("   %7.1e   %7.1e\n", t1, t2)

      if pAp == 0.0
        if (ρ == 0.0) | (γ > 0)
          α = t1
          on_boundary = true
        elseif γ < 0.0
          α = t2
          on_boundary = true
        end
      elseif (γ < 0.0) | (α < 0.0)
        α = t2
        on_boundary = true
      elseif (pAp < 0.0) | (α > t1)
        α = t1
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
    γ = rNorm * rNorm + β * γ - β * α * oldpAp # p'r

    if (γ < 0.0) & (pAp > 0.0)
      solved = true
    elseif γ == 0.0
      solved = true
    end

  end
  verbose && @printf("\n")

  return x
end
