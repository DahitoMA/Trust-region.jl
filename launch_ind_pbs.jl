using OptimizationProblems
import CSV
import DataFrames

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
## Dimension > 1
# Problems = [beale, brownbs, broydn7d, brybnd,
# chainwoo, chnrosnb_mod, cosine, cragglvy, curly, dixmaane, dixmaanf, dixmaang,
# dixmaanh, dixmaani, dixmaanj, dixmaank, dixmaanl, dixmaanm, dixmaann, dixmaano, dixmaanp,
# edensch, eg2, errinros_mod, extrosnb, fletcbv2, fletcbv3_mod, fletchcr, fminsrf2, freuroth, genhumps, genrose, genrose_nash,
# indef_mod, liarwhd, meyer3, morebv, ncb20, ncb20b, noncvxu2, noncvxun, nondia, penalty2, powellsg, schmvett,
# scosine, sinquad, sparsine, sparsqur, srosenbr, tointgss, tquartic, woods]

# Nonconvex problems var ≥ 10 (50 pbs)
Problems = [broydn7d, brybnd, chainwoo, chnrosnb_mod, cosine, cragglvy, curly,
dixmaane, dixmaanf, dixmaang, dixmaanh, dixmaani, dixmaanj, dixmaank, dixmaanl,
dixmaanm, dixmaann, dixmaano, dixmaanp, edensch, eg2, errinros_mod, extrosnb,
fletcbv2, fletcbv3_mod, fletchcr, fminsrf2, freuroth, genhumps, genrose,
genrose_nash, indef_mod, liarwhd, morebv, ncb20, ncb20b, noncvxu2, noncvxun,
nondia, penalty2, powellsg, schmvett, scosine, sinquad, sparsine, sparsqur,
srosenbr, tointgss, tquartic, woods]

# algo = CG
algo = CR
csvname = string("noncvxTR",algo,".csv")

D = ["model     algo nvar   f(x)    f(x0)   ‖g(x)‖  ‖g(x0)‖   #f  #g  #Hv  #it #s_it #vs_it #rej_it optimal"]

MiniLogging.@info(loggerind, D[1])
n=100
for problem in Problems
    model = MathProgNLPModel(problem(n), name=string(problem))
    T = TrustRegion(model, algo)
    S = @sprintf("%5s %5s %4d %5s %5s %5s %5s %4d %4d %4d %4d %4d %4d %4d %5s", T[1], T[2], T[3], T[4], T[5], T[6], T[7], T[8], T[9], T[10], T[11], T[12], T[13], T[14], T[15])
    MiniLogging.@info(loggerind, S)
    D = vcat(D, S)
    reset!(model)
end
# DF = DataFrames.DataFrame(D)
# CSV.write(csvname, DF)
writedlm("nonconvpbsCR.txt", D)
