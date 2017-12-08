# An implemantation of the trust-region method
# TrustRegion(model,algo;filename,Δ,ϵa,ϵr,itemax) solves a continuous optimization problem 'model' within
# a trust region of initial radius Δ and with tolerances ϵa and ϵr.
# Steps are calculated using the argument 'algo', a truncated optimization method.

"""An implementation of the trust-region method
"""
function TrustRegion(model, algo; filename::String=string("result", string(algo)), Δ::Float64=10., ϵa::Float64=1e-6, ϵr::Float64=1e-6, itemax::Int=10000)

    @info(loggerTR, @sprintf("TrustRegion: resolution of %s using %s, initial radius = %8.1e", model.meta.name, string(algo), Δ))
    x = model.meta.x0 # initial estimation from the model
	n = model.meta.nvar # size of the problem
	g = grad(model, x) # ∇f(x_0)
	H = hess_op(model, x) # ∇²f(x_0)

	normg0 = norm(g,2) # ‖∇f(x_0)‖
	normg = normg0 # ‖∇f(x_k)‖
	ϵ = ϵa + ϵr * normg # tolerance

	fx0 = obj(model, x) # f(x_0)
	fx = fx0 # f(x_k)

	fvalues = [fx] # improving values of the objective
	k = 0 # number of iterations
	ite = [] # number of "successfull" iterations
	vs_ite = [] # very successfull iterations
	mvalues = [fx] # values of the quadratic model

    @debug(loggerTR, @sprintf("%4s  %9s  %7s  %7s", "Iter", "f", "‖∇f‖", "Radius"))
    @debug(loggerTR, @sprintf("%4d  %9.2e  %7.1e  %7.1e", k, fx, normg, Δ))

    jldname = string(filename, ".jld")
	txtname = string(filename, ".txt")

	while normg > ϵ && k < itemax # stopping criterion : ‖∇f(x_k)‖ <= ϵ or k >= itemax
		k = k + 1
		s = algo(H, -g, Δ) # step
		xtrial = x + s # x_k + s
		fxtrial = obj(model, xtrial) # f(x_k + s)
        @debug(loggerTR, @sprintf("fxtrial = %8.1e", fxtrial))
		Δf = fxtrial - fx # f(x_k + s) - f(x_k)
        @debug(loggerTR, @sprintf("fx = %8.1e", fx))

		Hs = H * s
		Δm = dot(g,s) + 0.5 * dot(s,Hs) # m_k(x_k + s) - m_k(x_k)
        @debug(loggerTR, @sprintf("Δm = %8.1e", Δm))
        abs(Δm) < eps(Float64) && @critical(loggerTR, "Δm ≃ 0")

		# ratio ρ
		ρ = Δf / Δm
        @debug(loggerTR, @sprintf("ratio Δf/Δq = %8.1e", ρ))

		# adaptation of the trust region
		if ρ < 1.e-4
            if Δ > eps()
                Δ = Δ / 3
            end
            @debug(loggerTR, @sprintf("Δf/Δq < 1e-4, new radius = %8.1e", Δ))
		else
            @debug(loggerTR, "Δf/Δq ≥ 1e-4")
			x = xtrial
			fx = fxtrial
			g = grad(model,x) # ∇f(x_k)
			normg = norm(g,2)
			H = hess_op(model, x) # ∇²f(x_k)
			fvalues = push!(fvalues, fx)
			mvalues = push!(mvalues, Δm + fx)
			ite = push!(ite, k)
			if ρ >= 0.99 && ρ >= 1.e-4
				Δ = 3 * Δ
                @debug(loggerTR, @sprintf("Δf/Δq ≥ 0.99, new radius = %8.1e", Δ))
				vs_ite = push!(vs_ite, k)
			end
            ρ < 0.99 && @debug(loggerTR, @sprintf("1e-4 ≤ Δf/Δq < 0.99, new radius = %8.1e", Δ))
		end

        @debug(loggerTR, @sprintf("%4s  %9s  %7s  %7s  %8s", "Iter", "f", "‖∇f‖", "Radius", "Ratio"))
		@debug(loggerTR, @sprintf("%4d  %9.2e  %7.1e  %7.1e  %7.1e", k, fx, normg, Δ, ρ))

	end

	# save(jldname,"x",x,"ite",ite,"fvalues",fvalues,"mvalues",mvalues)
	# X = load(jldname)
	# writedlm(txtname,X)

  	@info(loggerTR, @sprintf("%30s %s %9s %9s %9s %9s %3s %3s %4s %4s %4s %4s %4s","name", "nvar", "f(x*)", "/ f(x0)", "‖∇f(x*)‖", "/ ‖∇f(x0)‖", "#f", "#g", "#Hv", "#it", "s_it", "vs_it", "reject_it"))
  	@info(loggerTR, @sprintf("%30s %d  %8.1e  %8.1e    %7.1e  %7.1e     %d   %d   %d    %d    %d   %d     %d", model.meta.name, n, fx, fx0, normg, normg0, neval_obj(model), neval_grad(model), neval_hprod(model), k, length(ite)-length(vs_ite),length(vs_ite), k-length(ite)))

    optimal = normg ≤ ϵ
    tired = k ≥ itemax
    status = optimal ? :optimal : :tired

    # return x, fx, normg, k, optimal, tired, status # for use of two_solvers()
    # return fx, fx0, normg, normg0, k, status, neval_obj(model), neval_grad(model), neval_hprod(model), length(ite)-length(vs_ite), length(vs_ite), k-length(ite) # for CUTEstpbs
    # return [@sprintf("%5s %5s %4d %8.1e %8.1e %7.1e %7.1e %4d %4d %4d %4d %4d %4d %4d", model.meta.name[22:end], string(algo), n, fx, fx0, normg, normg0, neval_obj(model), neval_grad(model), neval_hprod(model), k, length(ite)-length(vs_ite), length(vs_ite), k-length(ite))]
    return [@sprintf("%5s %5s %4d %8.1e %8.1e %7.1e %7.1e %4d %4d %4d %4d %4d %4d %4d", model.meta.name, string(algo), n, fx, fx0, normg, normg0, neval_obj(model), neval_grad(model), neval_hprod(model), k, length(ite)-length(vs_ite), length(vs_ite), k-length(ite))]
end
