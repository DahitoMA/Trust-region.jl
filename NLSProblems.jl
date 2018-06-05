using NLSProblems
using Krylov

#Problems without mgh11
Problems = [LVcon501, LVcon502, LVcon503, LVcon504, LVcon511, LVcon512,
LVcon513, LVcon514, LVcon515, LVcon516, LVcon517, LVcon518, NZF1,
hs01, hs02, hs06, hs13, hs14, hs16, hs17, hs18, hs20, hs21, hs22,
hs23, hs25, hs26, hs27, hs28, hs30, hs31, hs32, hs42, hs43, hs46, hs48, hs49,
hs50, hs51, hs52, hs53, hs57, hs60, hs61, hs65, hs77, hs79,
mgh01, mgh02, mgh03,mgh04, mgh05, mgh06, mgh07, mgh08, mgh09,
mgh10, mgh12, mgh13, mgh14, mgh15, mgh16, mgh17, mgh18, mgh19,
mgh20, mgh21, mgh22, mgh23, mgh24, mgh25, mgh26, mgh27, mgh28, mgh29,
mgh30, mgh31, mgh32, mgh33, mgh34, mgh35]

Algos = [lsqr, lsmr] # [cgls, lsqr, crls, lsmr]

D = ["model     algo nvar   f(x)    f(x0)   ‖g(x)‖  ‖g(x0)‖   #r  #Av  #A'v  #it #s_it #vs_it #rej_it"]
@info(loggerNLS, D[1])

for problem in Problems
    model = problem() #MathProgNLSModel(problem())
    for algo in Algos
        T = TrustRegionLS(model, algo)
        S = @sprintf("%5s %5s %4d %8.1e %8.1e %7.1e %7.1e %4d %4d %4d %4d %4d %4d %4d", T[1], T[2], T[3], T[4], T[5], T[6], T[7], T[8], T[9], T[10], T[11], T[12], T[13], T[14])
        @info(loggerNLS, S)
        D = vcat(D, S)
        reset!(model)
    end
end
txtname = "NLSProblems.txt"
writedlm(txtname, D)

# lines = open(readlines, txtname)
# V = []
# i = 0
# for l in lines
#     i += 1
#     if i < 35
#         l = replace(l, "Generic", string("hs", Int(ceil((i-1)/4))))
#     else
#         l = replace(l, "Generic", string("mgh", Int(ceil((i-1)/4))))
#     end
#     V = vcat(V, l)
# end
# writedlm(txtname, V)
