using Krylov

# convex
# ACG = readdlm("NewResults_TR/convpbsCG.txt")
# ACR = readdlm("NewResults_TR/convpbsCR.txt")

# nonconvex
ACG = readdlm("NewResults_TR/nonconvpbsCG.txt")
ACR = readdlm("NewResults_TR/nonconvpbsCR.txt")

Problems = ACG[2:end,1]

function perfprofileOpt()
    stats = Dict{String, Any}("CR" => ["#f" "#g" "#Hv"], "CG" => ["#f" "#g" "#Hv"])
    k = 1
    for p in Problems
        k += 1
        statCR = ACR[k,15] ? ACR[k,8:10] : -ACR[k,8:10]
        stats["CR"] = vcat(values(stats["CR"]), statCR')
        statCG = ACG[k,15] ? ACG[k,8:10] : -ACG[k,8:10]
        stats["CG"] = vcat(values(stats["CG"]), statCG')
    end
    return stats
end

stats = perfprofileOpt()
tf = font(17) # titlefont
f = font(17)

pb_type = "nonconvex problems"

p = performance_profile(hcat([p[2:end, 1] for p in values(stats)]...),
                            collect(String, [string(s) for s in keys(stats)]),
                            title=string("#f (", pb_type, ")"),
                            titlefont = tf, legendfont = f, guidefont = f,
                            legend=:bottomright, palette=:blues, w=2, style=:auto) # Profile for #f
# Plots.xlabel!("Within this factor of the best (log₂ scale)")
savefig(p, string("profil_f_", pb_type, ".pdf"))

p = performance_profile(hcat([p[2:end, 2] for p in values(stats)]...),
                            collect(String, [string(s) for s in keys(stats)]),
                            title=string("#g (", pb_type, ")"),
                            titlefont = tf, legendfont = f, guidefont = f,
                            legend=:bottomright, palette=:blues, w=2, style=:auto) # Profile for #g
# Plots.xlabel!("Within this factor of the best (log₂ scale)")
savefig(p, string("profil_g_", pb_type, ".pdf"))


p = performance_profile(hcat([p[2:end, 3] for p in values(stats)]...),
                            collect(String, [string(s) for s in keys(stats)]),
                            title=string("#Hv (", pb_type, ")"),
                            titlefont = tf, legendfont = f, guidefont = f,
                            legend=:bottomright, palette=:blues, w=2, style=:auto) # Profile for #Hv
# Plots.xlabel!("Within this factor of the best (log₂ scale)")
Plots.savefig(p, string("profil_Hv_", pb_type, ".pdf"))

p = performance_profile(hcat([p[2:end, 1]+p[2:end, 2]+p[2:end, 3] for p in values(stats)]...),
                            collect(String, [string(s) for s in keys(stats)]),
                            title=string("#f + #g  + #Hv (", pb_type, ")"),
                            titlefont = tf, legendfont = f, guidefont = f,
                            legend=:bottomright, palette=:blues, w=2, style=:auto) # Profile for #f + #g + #Hv
# Plots.xlabel!("Within this factor of the best (log₂ scale)")
savefig(p, string("profil_f+g+Hv_", pb_type, ".pdf"))
