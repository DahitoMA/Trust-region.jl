using Optimize
using CUTEst
using Krylov
using OptimizationProblems

# CUTEst problems
Problems = open(readlines, "Problems.txt")
# Problems = deleteat!(Problems,56)

n = 100

function two_solversCUTEst()
    solvers = Dict{Symbol, Function}(:CG => TRCG_CUTEst,:CR => TRCR_CUTEst)
    bmark_args = Dict{Symbol, Any}(:skipif => model -> !unconstrained(model))

    bmark_solvers(solvers,
                (CUTEstModel(p) for p in Problems);
                bmark_args...) # for CUTEst problems

end


stats = two_solversCUTEst()
tf = font(13) # titlefont
f = font(12)
pb_type = "CUTEst problems"
algo_used = "CG and CR"
p = performance_profile(hcat([p[:, 1] for p in values(stats)]...),
                            collect(String, [string(s) for s in keys(stats)]),
                            title=string("Performance profile : #f in the resolution of ", pb_type, " with ", algo_used),
                            titlefont = tf, legendfont = f, guidefont = f, size = (800,800), legend=:bottomright) # Profile for #f

Plots.savefig(p, string("profil_f_", pb_type, ".pdf"))

p = performance_profile(hcat([p[:, 2] for p in values(stats)]...),
                            collect(String, [string(s) for s in keys(stats)]),
                            title=string("Performance profile : #g in the resolution of ", pb_type, " with ", algo_used),
                            titlefont = tf, legendfont = f, guidefont = f, size = (800,800), legend=:bottomright) # Profile for #g

Plots.savefig(p, string("profil_g_", pb_type, ".pdf"))


p = performance_profile(hcat([p[:, 3] for p in values(stats)]...),
                            collect(String, [string(s) for s in keys(stats)]),
                            title=string("Performance profile : #Hv in the resolution of ", pb_type, " with ", algo_used),
                            titlefont = tf, legendfont = f, guidefont = f, size = (800,800), legend=:bottomright) # Profile for #Hv

Plots.savefig(p, string("profil_Hv_", pb_type, ".pdf"))

p = performance_profile(hcat([sum(p, 2) for p in values(stats)]...),
                            collect(String, [string(s) for s in keys(stats)]),
                            title=string("Performance profile : #f + #g  + #Hv in the resolution of ", pb_type, " with ", algo_used),
                            titlefont = tf, legendfont = f, guidefont = f, size = (800,800), legend=:bottomright) # Profile for #f + #g + #Hv

Plots.savefig(p, string("profil_f+g+Hv_", pb_type, ".pdf"))
