# An implemantation of the trust-region method for least squares problems
# TrustRegionLS(lsmodel,algo;filename,Δ,ϵa,ϵr,itemax) solves a continuous optimization problem 'lsmodel' within
# a trust region of initial radius Δ and with tolerances ϵa and ϵr.
# Steps are calculated using the argument 'algo', a truncated optimization method.

"""An implementation of the trust-region method for least squares problems
"""
function TrustRegionLS(lsmodel, algo; filename::String=string("result", string(algo)), Δ::Float64=10., ϵa::Float64=1e-6, ϵr::Float64=1e-6, itemax::Int=10000)#, args...)

    @info(loggerTRLS, @sprintf("TrustRegionLS: resolution of %s using %s, initial radius = %8.1e", lsmodel.meta.name, string(algo), Δ))

    n = lsmodel.meta.nvar # size of the problem
    x = lsmodel.meta.x0 # initial estimation from the lsmodel
    r = residual(lsmodel, x) # F(x_0)
    A = jac_op_residual(lsmodel, x) # J(x_0)
    g =  A' * r # ∇f(x_0) where f(x) = ½ ‖F(x)‖²

    normg0 = norm(g,2) # ‖∇f(x_0)‖
	normg = normg0 # ‖∇f(x_k)‖
	ϵ = ϵa + ϵr * normg # tolerance

	fx0 = 0.5 * norm(r)^2 # f(x_0)
    fx = fx0 # f(x_k)

	fvalues = [fx] # improving values of the objective
	k = 0 # number of iterations
	ite = [] # number of "successfull" iterations
	vs_ite = [] # very successfull iterations
    γ = A * x + r
    γNorm² = norm(γ)^2
    m = 0.5 * γNorm²
	mvalues = [m] # values of the quadratic lsmodel

    @debug(loggerTRLS, @sprintf("%4s  %9s  %7s  %7s", "Iter", "f", "‖∇f‖", "Radius"))
    @debug(loggerTRLS, @sprintf("%4d  %9.2e  %7.1e  %7.1e", k, fx, normg, Δ))

    jldname = string(filename, ".jld")
	txtname = string(filename, ".txt")

	while normg > ϵ && k < itemax # stopping criterion : ‖∇f(x_k)‖ <= ϵ or k >= itemax
		k = k + 1
		s = algo(A, -r, atol = 0., rtol = min(0.1, sqrt(normg)), radius = Δ , verbose = false) # step
        s = s[1]
		xtrial = x + s # x_k + s
        rtrial = residual(lsmodel, xtrial) # F(x_k + s)
		fxtrial =  0.5 * norm(rtrial)^2 # f(x_k + s)
        @debug(loggerTRLS, @sprintf("fxtrial = %8.1e", fxtrial))
		Δf = fxtrial - fx # f(x_k + s) - f(x_k)
        @debug(loggerTRLS, @sprintf("Δf = %8.1e", Δf))

		As = A * s
        Δm = 0.5 * norm(As)^2 + dot(s, g) # m_k(x_k + s) - m_k(x_k)
        @debug(loggerTRLS, @sprintf("Δm = %8.1e", Δm))
        abs(Δm) < eps(Float64) && @critical(loggerTRLS, "Δm ≃ 0")

		# ratio ρ
		ρ = Δf / Δm
        @debug(loggerTRLS, @sprintf("ratio Δf/Δq = %8.1e", ρ))

		# adaptation of the trust region
		if ρ < 1.e-4
            Δ = Δ / 3
            @debug(loggerTRLS, @sprintf("Δf/Δq < 1e-4, new radius = %8.1e", Δ))
		else
            @debug(loggerTRLS, "Δf/Δq ≥ 1e-4")
			x = xtrial
			fx = fxtrial
            r = rtrial
            A = jac_op_residual(lsmodel, x)
            g = A' * r # ∇f(x_k)
			normg = norm(g)
			fvalues = push!(fvalues, fx)
			mvalues = push!(mvalues, Δm + 0.5 * γNorm²)
            γ = A * x + r
            γNorm² = norm(γ)^2
			ite = push!(ite, k)
			if ρ >= 0.99 && ρ >= 1.e-4
				Δ = 3 * Δ
                @debug(loggerTRLS, @sprintf("Δf/Δq ≥ 0.99, new radius = %8.1e", Δ))
				vs_ite = push!(vs_ite, k)
			end
            ρ < 0.99 && @debug(loggerTRLS, @sprintf("1e-4 ≤ Δf/Δq < 0.99, new radius = %8.1e", Δ))
		end

        @debug(loggerTRLS, @sprintf("%4s  %9s  %7s  %7s  %8s", "Iter", "f", "‖∇f‖", "Radius", "Ratio"))
		@debug(loggerTRLS, @sprintf("%4d  %9.2e  %7.1e  %7.1e  %7.1e", k, fx, normg, Δ, ρ))

	end

	# save(jldname,"x",x,"ite",ite,"fvalues",fvalues,"mvalues",mvalues)
	# X = load(jldname)
	# writedlm(txtname,X)

  	@info(loggerTRLS, @sprintf("%30s %s %9s %9s %9s %9s %3s %3s %3s %4s %4s %4s %4s","name", "nvar", "f(x*)", "/ f(x0)", "‖∇f(x*)‖", "/ ‖∇f(x0)‖", "#r", "#Av", "#A'v", "#it", "s_it", "vs_it", "reject_it"))
  	@info(loggerTRLS, @sprintf("%30s %d  %8.1e  %8.1e    %7.1e  %7.1e     %d   %d   %d    %d    %d   %d     %d", lsmodel.meta.name, n, fx, fx0, normg, normg0, neval_residual(lsmodel), neval_jprod_residual(lsmodel), neval_jtprod_residual(lsmodel), k, length(ite)-length(vs_ite),length(vs_ite), k-length(ite)))

    optimal = normg ≤ ϵ
    tired = k ≥ itemax
    status = optimal ? "optimal" : "tired"

    # return x, fx, normg, k, optimal, tired, status # for use of two_solvers()
    # return [@sprintf("%5s %5s %4d %8.1e %8.1e %7.1e %7.1e %4d %4d %4d %4d %4d %4d %4d", lsmodel.meta.name, string(algo)[8:end], n, fx, fx0, normg, normg0, neval_residual(lsmodel), neval_jprod_residual(lsmodel), neval_jtprod_residual(lsmodel), k, length(ite)-length(vs_ite), length(vs_ite), k-length(ite))]
    return [neval_residual(lsmodel) neval_jprod_residual(lsmodel) neval_jtprod_residual(lsmodel)] #for NLStwo_solvers()
end
