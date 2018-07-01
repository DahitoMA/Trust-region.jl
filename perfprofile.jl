using Krylov

ACG = readdlm("TRPbs/AllPbsCG.txt")
ACR = readdlm("TRPbs/AllPbsCR.txt")
Problems = ACG[2:end,1]

function perfprofile()
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

stats = perfprofile()
tf = font(10) # titlefont
f = font(10)
pb_type = "allpbs"

p = performance_profile(hcat([p[2:end, 1] for p in values(stats)]...),
                            collect(String, [string(s) for s in keys(stats)]),
                            title="Performance profile : #f",
                            titlefont = tf, legendfont = f, guidefont = f, legend=:bottomright) # Profile for #f

savefig(p, string("profil_f_", pb_type, ".pdf"))

p = performance_profile(hcat([p[2:end, 2] for p in values(stats)]...),
                            collect(String, [string(s) for s in keys(stats)]),
                            title="Performance profile : #g",
                            titlefont = tf, legendfont = f, guidefont = f, legend=:bottomright) # Profile for #g

savefig(p, string("profil_g_", pb_type, ".pdf"))


p = performance_profile(hcat([p[2:end, 3] for p in values(stats)]...),
                            collect(String, [string(s) for s in keys(stats)]),
                            title="Performance profile : #Hv",
                            titlefont = tf, legendfont = f, guidefont = f, legend=:bottomright) # Profile for #Hv

Plots.savefig(p, string("profil_Hv_", pb_type, ".pdf"))

p = performance_profile(hcat([p[2:end, 1]+p[2:end, 2]+p[2:end, 3] for p in values(stats)]...),
                            collect(String, [string(s) for s in keys(stats)]),
                            title="Performance profile : #f + #g + #Hv",
                            titlefont = tf, legendfont = f, guidefont = f, legend=:bottomright) # Profile for #f + #g + #Hv

savefig(p, string("profil_f+g+Hv_", pb_type, ".pdf"))
