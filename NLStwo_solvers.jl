using Optimize
using Krylov
using NLSProblems

# NLSProblems without mgh11
Problems = [mgh01, mgh02, mgh03,mgh04, mgh05, mgh06, mgh07, mgh08, mgh09,
mgh10, mgh12, mgh13, mgh14, mgh15, mgh16, mgh17, mgh18, mgh19,
mgh20, mgh21, mgh22, mgh23, mgh24, mgh25, mgh26, mgh27, mgh28, mgh29,
mgh30, mgh31, mgh32, mgh33, mgh34, mgh35]

function NLStwo_solvers() #for NLSProblems
    solvers = Dict{String, Function}("cgls" => TRCGLS,"crls" => TRCRLS)
    stats = Dict{String, Any}("cgls" => ["#r" "#Av" "A'v"], "crls" => ["#r" "#Av" "A'v"])
    # solvers = Dict{String, Function}("lsqr" => TRLSQR,"lsmr" => TRLSMR)
    # stats = Dict{String, Any}("lsqr" => ["#r" "#Av" "A'v"], "lsmr" => ["#r" "#Av" "A'v"])
    for (name, solver) in solvers
        for p in Problems
            model = ADNLSModel(p())
            stat = solver(model)
            stats[name] = vcat(values(stats[name]), stat)
            reset!(model)
        end
    end
    return stats
end

stats = NLStwo_solvers()
tf = font(17) # titlefont
f = font(15)
pb_type = "NLSProblems"
algo_used = "cgls and crls"
# algo_used = "lsqr and lsmr"

p = performance_profile(hcat([p[2:end, 1] for p in values(stats)]...),
                            collect(String, [string(s) for s in keys(stats)]),
                            title=string("Performance profile : #r in the resolution of ", pb_type, " with ", algo_used),
                            titlefont = tf, legendfont = f, guidefont = f, size = (800,800)) # Profile for #r

# savefig(p, string("profil_r_", pb_type, ".pdf"))
#
# p = performance_profile(hcat([p[2:end, 2] for p in values(stats)]...),
#                             collect(String, [string(s) for s in keys(stats)]),
#                             title=string("Performance profile : #Av in the resolution of ", pb_type, " with ", algo_used),
#                             titlefont = tf, legendfont = f, guidefont = f, size = (800,800)) # Profile for #Av
#
# savefig(p, string("profil_Av_", pb_type, ".pdf"))
#
#
# p = performance_profile(hcat([p[2:end, 3] for p in values(stats)]...),
#                             collect(String, [string(s) for s in keys(stats)]),
#                             title=string("Performance profile : #A'v in the resolution of ", pb_type, " with ", algo_used),
#                             titlefont = tf, legendfont = f, guidefont = f, size = (800,800)) # Profile for #A'v
#
# Plots.savefig(p, string("profil_A'v_", pb_type, ".pdf"))
#
# p = performance_profile(hcat([sum(p[2:end,:], 2) for p in values(stats)]...),
#                             collect(String, [string(s) for s in keys(stats)]),
#                             title=string("Performance profile : #r + #Av + A'v in the resolution of ", pb_type, " with ", algo_used),
#                             titlefont = tf, legendfont = f, guidefont = f, size = (800,800)) # Profile for #r + #Av + #A'v
#
# savefig(p, string("profil_r+Av+A'v_", pb_type, ".pdf"))
