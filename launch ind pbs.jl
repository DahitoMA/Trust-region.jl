Problems = [AMPGO02, AMPGO03, AMPGO04, AMPGO05, AMPGO06, AMPGO08, AMPGO09,
AMPGO10, AMPGO11, AMPGO12, AMPGO13, AMPGO14, AMPGO15, AMPGO18,
AMPGO20, AMPGO21, AMPGO22,
beale, brownbs, broydn7d, brybnd,
chainwoo, chnrosnb_mod, cosine, cragglvy, curly, dixmaane, dixmaanf, dixmaang,
dixmaanh, dixmaani, dixmaanj, dixmaank, dixmaanl, dixmaanm, dixmaann, dixmaano, dixmaanp,
Dus2_1, Dus2_3, Dus2_9, Duscube,
edensch, eg2, errinros_mod, extrosnb, fletcbv2, fletcbv3_mod, fletchcr, fminsrf2, freuroth, genhumps, genrose, genrose_nash,
indef_mod, liarwhd, meyer3, morebv, ncb20, ncb20b, noncvxu2, noncvxun, nondia, penalty2, powellsg, schmvett,
scosine, Shpak1, Shpak2, Shpak3, Shpak4, Shpak5, Shpak6, sinquad, sparsine, sparsqur, srosenbr, tointgss, tquartic, woods]

# penalty3, sbrybnd

Algos = [CG, CR]

loggerInd = MiniLogging.Logger("loggerInd", MiniLogging.INFO)
push!(loggerInd.handlers, MiniLogging.Handler(open("indpbs.txt", "w"), "%Y-%m-%d %H:%M:%S"))

for problem in Problems
    model = MathProgNLPModel(problem(), name=string(problem))
    for algo in Algos
        @info(loggerInd, TrustRegion(model, algo))
        reset!(model)
  end
end
