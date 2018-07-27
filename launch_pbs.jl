using OptimizationProblems

## All pbs dimension > 1 and not in CUTEst
# Problems = [brownden, cliff, clplatea, clplateb, clplatec, nasty,
# palmer1c, palmer1d, palmer2c, palmer3c, palmer4c, palmer5c, palmer5d,
# palmer6c, palmer7c, palmer8c,
# beale, brownbs,
# chnrosnb_mod,
# errinros_mod, fletcbv3_mod, genrose_nash,
# meyer3, scosine]

# All pbs var ≥ 10 and not in CUTEst
Problems = [clplatea, clplateb, clplatec, chnrosnb_mod, errinros_mod,
fletcbv3_mod, genrose_nash, scosine]

algo = CG
# algo = CR
# csvname = string("noncvxTR",algo,".csv")

D = ["model     algo nvar   f(x)    f(x0)   ‖g(x)‖  ‖g(x0)‖   #f  #g  #Hv  #it #s_it #vs_it #rej_it optimal"]

for problem in Problems
    model = MathProgNLPModel(problem(), name=string(problem))
    T = TrustRegion(model, algo)
    S = @sprintf("%5s %5s %4d %5s %5s %5s %5s %4d %4d %4d %4d %4d %4d %4d %5s", T[1], T[2], T[3], T[4], T[5], T[6], T[7], T[8], T[9], T[10], T[11], T[12], T[13], T[14], T[15])
    D = vcat(D, S)
    reset!(model)
end
writedlm("OptPbsCG.txt", D)
