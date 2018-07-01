using Krylov
# using NLSProblems

# NLSProblems without mgh11
# Problems = [LVcon501, LVcon502, LVcon503, LVcon504, LVcon511, LVcon512,
# LVcon513, LVcon514, LVcon515, LVcon516, LVcon517, LVcon518, NZF1,
# hs01, hs02, hs06, hs13, hs14, hs16, hs17, hs18, hs20, hs21, hs22,
# hs23, hs25, hs26, hs27, hs28, hs30, hs31, hs32, hs42, hs43, hs46, hs48, hs49,
# hs50, hs51, hs52, hs53, hs57, hs60, hs61, hs65, hs77, hs79,
# mgh01, mgh02, mgh03,mgh04, mgh05, mgh06, mgh07, mgh08, mgh09,
# mgh10, mgh12, mgh13, mgh14, mgh15, mgh16, mgh17, mgh18, mgh19,
# mgh20, mgh21, mgh22, mgh23, mgh24, mgh25, mgh26, mgh27, mgh28, mgh29,
# mgh30, mgh31, mgh32, mgh33, mgh34, mgh35]

ALSQR = readdlm("LSPbs/AllPbsLSQR.txt")
ALSMR = readdlm("LSPbs/AllPbsLSMR.txt")
Problems = ALSQR[2:end,1]

function NLStwo_solvers() #for NLSProblems
    # solvers = Dict{String, Function}("cgls" => TRCGLS,"crls" => TRCRLS)
    # stats = Dict{String, Any}("cgls" => ["#r" "#Av" "#A'v"], "crls" => ["#r" "#Av" "#A'v"])
    stats = Dict{String, Any}("LSQR" => ["#r" "#Av" "#A'v"], "LSMR" => ["#r" "#Av" "#A'v"])
    # for (name, solver) in solvers
    #     for p in Problems
    #         model = p()
    #         T = solver(model)
    #         stat = T[15] ? T[8:10] : -T[8:10] # T[15] = optimal
    #         stats[name] = vcat(values(stats[name]), stat')
    #         reset!(model)
    #     end
    # end
    k = 1
    for p in Problems
        k += 1
        statLSMR = ALSMR[k,15] ? ALSMR[k,8:10] : -ALSMR[k,8:10]
        stats["LSMR"] = vcat(values(stats["LSMR"]), statLSMR')
        statLSQR = ALSQR[k,15] ? ALSQR[k,8:10] : -ALSQR[k,8:10]
        stats["LSQR"] = vcat(values(stats["LSQR"]), statLSQR')
    end
    return stats
end

stats = NLStwo_solvers()
tf = font(10) # titlefont
f = font(10)
pb_type = "allpbs"
algo_used = "LSQR and LSMR"

p = performance_profile(hcat([p[2:end, 1] for p in values(stats)]...),
                            collect(String, [string(s) for s in keys(stats)]),
                            title=string("Performance profile : #F"),
                            titlefont = tf, legendfont = f, guidefont = f, legend=:bottomright) # Profile for #r

savefig(p, string("profil_r_", pb_type, ".pdf"))

p = performance_profile(hcat([p[2:end, 2]+p[2:end, 3] for p in values(stats)]...),
                            collect(String, [string(s) for s in keys(stats)]),
                            title=string("Performance profile : #Jv + #Jᵀv"),
                            titlefont = tf, legendfont = f, guidefont = f, legend=:bottomright) # Profile for #Av + #A'v

savefig(p, string("profil_Av+A'v_", pb_type, ".pdf"))

p = performance_profile(hcat([p[2:end, 1]+p[2:end, 2]+p[2:end, 3] for p in values(stats)]...),
                            collect(String, [string(s) for s in keys(stats)]),
                            title=string("Performance profile : #F + #Jv + #Jᵀv"),
                            titlefont = tf, legendfont = f, guidefont = f, legend=:bottomright) # Profile for #r + #Av + #A'v

savefig(p, string("profil_r+Av+A'v_", pb_type, ".pdf"))
