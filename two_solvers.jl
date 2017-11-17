using Optimize
using CUTEst
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

# All unconstrained problems
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

#CUTEst problems
Problems = open(readlines, "Problems.txt")

n = 100

function two_solvers()
    solvers = [TRCG, TRCR]
    # solvers = [TRCG, TRminres]
    # solvers = [TRCR, TRminres]
    bmark_args = Dict{Symbol, Any}(:skipif => model -> !unconstrained(model))
    profile_args = Dict{Symbol, Any}(:title => "f+g+hprod")
    bmark_and_profile(solvers,
                    (MathProgNLPModel(p(n), name=string(p)) for p in Problems),
                    bmark_args=bmark_args, profile_args=profile_args)
    # bmark_and_profile(solvers,
    #                 (CUTEstModel(p(n), name=string(p)) for p in Problems),
    #                 bmark_args=bmark_args, profile_args=profile_args)
end

stats, profile = two_solvers()
tf = font(15) # titlefont
f = font(12)
performance_profile(hcat([p[:, 1] for p in values(stats)]...),
                            collect(String, [string(s) for s in keys(stats)]),
                            title="Performance profile : #f in the resolution of convex and nonconvex problems with TRCG and TRCR",
                            titlefont = tf, legendfont = f, guidefont = f) # Profile for #f

# performance_profile(hcat([p[:, 2] for p in values(stats)]...),
#                             collect(String, [string(s) for s in keys(stats)]),
#                             title="Performance profile: #g in the resolution of convex and nonconvex problems with TRCG and TRCR",
#                             titlefont = tf, legendfont = f, guidefont = f) # Profile for #g
#
# performance_profile(hcat([p[:, 3] for p in values(stats)]...),
#                             collect(String, [string(s) for s in keys(stats)]),
#                             title="Performance profile: #Hv in the resolution of nonconvex problems with TRCG and TRCR",
#                             titlefont = tf, legendfont = f, guidefont = f) # Profile for #Hv
#
# performance_profile(hcat([sum(p, 2) for p in values(stats)]...),
#                             collect(String, [string(s) for s in keys(stats)]),
#                             title="Performance profile: #f + #g + #Hv in the resolution of nonconvex problems with TRCG and TRCR",
#                             titlefont = tf, legendfont = f, guidefont = f) # Profile for #f + #g + #Hv
