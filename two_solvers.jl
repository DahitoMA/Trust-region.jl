using Optimize
using CUTEst
using Krylov

# non convex problems
# Problems = [AMPGO02, AMPGO03, AMPGO04, AMPGO05, AMPGO06, AMPGO08, AMPGO09,
# AMPGO10, AMPGO11, AMPGO12, AMPGO13, AMPGO14, AMPGO15, AMPGO18,
# AMPGO20, AMPGO21, AMPGO22,
# beale, brownbs, broydn7d, brybnd,
# chainwoo, chnrosnb_mod, cosine, cragglvy, curly, dixmaane, dixmaanf, dixmaang,
# dixmaanh, dixmaani, dixmaanj, dixmaank, dixmaanl, dixmaanm, dixmaann, dixmaano, dixmaanp,
# Dus2_1, Dus2_3, Dus2_9, Duscube,
# edensch, eg2, errinros_mod, extrosnb, fletcbv2, fletcbv3_mod, fletchcr, fminsrf2, freuroth, genhumps, genrose, genrose_nash,
# indef_mod, liarwhd, meyer3, morebv, ncb20, ncb20b, noncvxu2, noncvxun, nondia, penalty2, powellsg, schmvett,
# scosine, Shpak1, Shpak2, Shpak3, Shpak4, Shpak5, Shpak6, sinquad, sparsine, sparsqur, srosenbr, tointgss, tquartic, woods]

# convex problems
# Problems = [arglina, arglinb, arglinc, arwhead, bdqrtic, brownden, cliff, clplatea, clplateb, clplatec, dixon3dq, dqdrtic, dqrtic, engval1, nasty, nondquar,
# palmer1c, palmer1d, palmer2c, palmer3c, palmer4c, palmer5c, palmer5d, palmer6c, palmer7c, palmer8c, power, quartc, tridia, vardim]

# unconstrained OptimizationProblems
# Problems = [AMPGO02, AMPGO03, AMPGO04, AMPGO05, AMPGO06, AMPGO08, AMPGO09,
# AMPGO10, AMPGO11, AMPGO12, AMPGO13, AMPGO14, AMPGO15, AMPGO18,
# AMPGO20, AMPGO21, AMPGO22,
# beale, brownbs, broydn7d, brybnd,
# chainwoo, chnrosnb_mod, cosine, cragglvy, curly, dixmaane, dixmaanf, dixmaang,
# dixmaanh, dixmaani, dixmaanj, dixmaank, dixmaanl, dixmaanm, dixmaann, dixmaano, dixmaanp,
# Dus2_1, Dus2_3, Dus2_9, Duscube,
# edensch, eg2, errinros_mod, extrosnb, fletcbv2, fletcbv3_mod, fletchcr, fminsrf2, freuroth, genhumps, genrose, genrose_nash,
# indef_mod, liarwhd, meyer3, morebv, ncb20, ncb20b, noncvxu2, noncvxun, nondia, penalty2, powellsg, schmvett,
# scosine, Shpak1, Shpak2, Shpak3, Shpak4, Shpak5, Shpak6, sinquad, sparsine, sparsqur, srosenbr, tointgss, tquartic, woods,
# arglina, arglinb, arglinc, arwhead, bdqrtic, brownden, cliff, clplatea, clplateb, clplatec, dixon3dq, dqdrtic, dqrtic, engval1, nasty, nondquar,
# palmer1c, palmer1d, palmer2c, palmer3c, palmer4c, palmer5c, palmer5d, palmer6c, palmer7c, palmer8c, power, quartc, tridia, vardim]

# CUTEst problems
# Problems = open(readlines, "Problems.txt")
# Problems = deleteat!(Problems,56)

# NLSProblems without mgh11
Problems = [mgh01, mgh02, mgh03,mgh04, mgh05, mgh06, mgh07, mgh08, mgh09,
mgh10, mgh12, mgh13, mgh14, mgh15, mgh16, mgh17, mgh18, mgh19,
mgh20, mgh21, mgh22, mgh23, mgh24, mgh25, mgh26, mgh27, mgh28, mgh29,
mgh30, mgh31, mgh32, mgh33, mgh34, mgh35]

n = 100

function two_solvers()
    # solvers = Dict{Symbol, Function}(:CG => TRCG,:CR => TRCR)
    # solvers = Dict{Symbol, Function}(:CG => TRCG,:MINRES => TRminres)
    # solvers = Dict{Symbol, Function}(:CR => TRCR,:MINRES => TRminres)
    solvers = Dict{Symbol, Function}(:cgls => TRCGLS,:crls => TRCRLS)
    bmark_args = Dict{Symbol, Any}(:skipif => model -> !unconstrained(model))

    # bmark_solvers(solvers,
    #             (CUTEstModel(p) for p in Problems);
    #             bmark_args...) # for CUTEst problems

    # bmark_solvers(solvers,
    #             (MathProgNLPModel(p(n), name=string(p)) for p in Problems);
    #             bmark_args...) # for OptimizationProblems

end


stats = two_solvers()
tf = font(15) # titlefont
f = font(12)
# pb_type = "CUTEst problems"
pb_type = "NLSProblems"
algo_used = "cgls and crls"
p = performance_profile(hcat([p[:, 1] for p in values(stats)]...),
                            collect(String, [string(s) for s in keys(stats)]),
                            title=string("Performance profile : #f in the resolution of ", pb_type, " with ", algo_used),
                            titlefont = tf, legendfont = f, guidefont = f) # Profile for #f

Plots.savefig(p, string("profil_f_", pb_type, ".pdf"))

# p = performance_profile(hcat([p[:, 2] for p in values(stats)]...),
#                             collect(String, [string(s) for s in keys(stats)]),
#                             title="Performance profile: #g in the resolution of nonconvex problems with TRCG and TRCR")
#                             # titlefont = tf, legendfont = f, guidefont = f) # Profile for #g
#
# Plots.savefig(p, string("profil_g_", pb_type, ".pdf"))
#
#
# p = performance_profile(hcat([p[:, 3] for p in values(stats)]...),
#                             collect(String, [string(s) for s in keys(stats)]),
#                             title="Performance profile: #Hv in the resolution of nonconvex problems with TRCG and TRCR")
#                             # titlefont = tf, legendfont = f, guidefont = f) # Profile for #Hv
#
# Plots.savefig(p, string("profil_Hv_", pb_type, ".pdf"))
#
# p = performance_profile(hcat([sum(p, 2) for p in values(stats)]...),
#                             collect(String, [string(s) for s in keys(stats)]),
#                             title="Performance profile: #f + #g + #Hv in the resolution of nonconvex problems with TRCG and TRCR")
#                             # titlefont = tf, legendfont = f, guidefont = f) # Profile for #f + #g + #Hv
#
# Plots.savefig(p, string("profil_f+g+Hv_", pb_type, ".pdf"))
