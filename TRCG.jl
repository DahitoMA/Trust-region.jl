# The trust-region method using a truncated version of the conjugate gradient method to find the steps.
# TRCG(model,Δ,ϵa,ϵr,itemax,verbose) solves a continous optimisation problem 'model' within
# a trust region and with tolerances ϵa and ϵr

"""The trust-region method using truncated conjugate gradient method
"""
function TRCG(model, Δ::Float64=10., ϵa::Float64=1e-6, ϵr::Float64=1e-6, itemax::Int=10000, verbose::Bool=false)

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

	jldname = "resultTRCG.jld"
	txtname = "resultTRCG.txt"

	while normg > ϵ && k < itemax # stopping criterion : ‖∇f(x_k)‖ <= ϵ or k >= itemax
		k = k + 1
		s = CG(H,-g,Δ) # step
		xtrial = x + s # x_k + s
		fxtrial = obj(model, xtrial) # f(x_k + s)
		Δf = fxtrial - fx # f(x_k + s) - f(x_k)

    Hs = H * s
    Δm = dot(g,s) + 1/2 * dot(s,Hs) # m_k(x_k + s) - m_k(x_k)

		# ratio ρ
		ρ = Δf / Δm

		# adaptation of the trust region
		if ρ < 1.e-4
			Δ = Δ / 3
		else
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
				vs_ite = push!(vs_ite, k)
			end
		end

		if verbose == true
			@printf("%4s  %9s  %7s  %7s  %8s\n", "Iter", "f", "‖∇f‖", "Radius", "Ratio")
			@printf("%4d  %9.2e  %7.1e  %7.1e  %7.1e\n", k, fx, normg, Δ, ρ)
		end

	end

	save(jldname,"x",x,"ite",ite,"fvalues",fvalues,"mvalues",mvalues)
	X = load(jldname)
	writedlm(txtname,X)
  @printf("%30s %s %9s %9s %9s %9s %3s %3s %4s %4s %4s %4s %4s\n","name", "nvar", "f(x*)", "/ f(x0)", "‖∇f(x*)‖", "/ ‖∇f(x0)‖", "#f", "#g", "#Hv", "#it", "s_it", "vs_it", "reject_it")
  @printf("%30s %d  %8.1e  %8.1e    %7.1e  %7.1e     %d   %d   %d    %d    %d   %d     %d\n", model.meta.name, n, fx, fx0, normg, normg0, neval_obj(model), neval_grad(model), neval_hprod(model), k, length(ite)-length(vs_ite),length(vs_ite), k-length(ite))
	return x

end
