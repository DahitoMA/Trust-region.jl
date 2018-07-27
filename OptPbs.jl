import Optimize
using Krylov
using OptimizationProblems

# Nonconvex problems var ≥ 10 (50 pbs)
Problems = [broydn7d, brybnd, chainwoo, chnrosnb_mod, cosine, cragglvy, curly,
dixmaane, dixmaanf, dixmaang, dixmaanh, dixmaani, dixmaanj, dixmaank, dixmaanl,
dixmaanm, dixmaann, dixmaano, dixmaanp, edensch, eg2, errinros_mod, extrosnb,
fletcbv2, fletcbv3_mod, fletchcr, fminsrf2, freuroth, genhumps, genrose,
genrose_nash, indef_mod, liarwhd, morebv, ncb20, ncb20b, noncvxu2, noncvxun,
nondia, penalty2, powellsg, schmvett, scosine, sinquad, sparsine, sparsqur,
srosenbr, tointgss, tquartic, woods]

# Convex problems var ≥ 10 (17 pbs)
# Problems = [arglina, arglinb, arglinc, arwhead, bdqrtic, clplatea, clplateb, clplatec,
# dixon3dq, dqdrtic, dqrtic, engval1, nondquar, power, quartc, tridia, vardim]

n = 100

function OptPbs()
    solvers = Dict{Symbol, Function}(:CR => TRCR,:CG => TRCG)
    bmark_args = Dict{Symbol, Any}(:skipif => model -> !unconstrained(model))

    Optimize.bmark_solvers(solvers,
                (MathProgNLPModel(p(n), name=string(p)) for p in Problems);
                bmark_args...)

end


stats = OptPbs()
tf = font(10) # titlefont
f = font(10)

pb_type = "nonconvex problems"
p = performance_profile(hcat([p[:, 1] for p in values(stats)]...),
                            collect(String, [string(s) for s in keys(stats)]),
                            title=string("Performance profile : #f in the solution of ", pb_type),
                            titlefont = tf, legendfont = f, guidefont = f, legend=:bottomright) # Profile for #f

Plots.savefig(p, string("profil_f_", pb_type, ".pdf"))

p = performance_profile(hcat([p[:, 2] for p in values(stats)]...),
                            collect(String, [string(s) for s in keys(stats)]),
                            title=string("Performance profile : #g in the solution of ", pb_type),
                            titlefont = tf, legendfont = f, guidefont = f, legend=:bottomright) # Profile for #g

Plots.savefig(p, string("profil_g_", pb_type, ".pdf"))


p = performance_profile(hcat([p[:, 3] for p in values(stats)]...),
                            collect(String, [string(s) for s in keys(stats)]),
                            title=string("Performance profile : #Hv in the solution of ", pb_type),
                            titlefont = tf, legendfont = f, guidefont = f, legend=:bottomright) # Profile for #Hv

Plots.savefig(p, string("profil_Hv_", pb_type, ".pdf"))

p = performance_profile(hcat([sum(p, 2) for p in values(stats)]...),
                            collect(String, [string(s) for s in keys(stats)]),
                            title=string("Performance profile : #f + #g  + #Hv in the solution of ", pb_type),
                            titlefont = tf, legendfont = f, guidefont = f, legend=:bottomright) # Profile for #f + #g + #Hv

Plots.savefig(p, string("profil_f+g+Hv_", pb_type, ".pdf"))
